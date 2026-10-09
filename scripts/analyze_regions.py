#!/usr/bin/env python3
"""Optional genomic windows and MACS3 hotspots from deduplicated BLISS ends."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.regions import consensus_regions, load_sites, union_regions, window_regions, write_matrix
from helixbusters.reporting import write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for option in ('samples', 'counts', 'headers', 'molecules'):
        p.add_argument(f'--{option}', nargs='+', required=True)
    p.add_argument('--design-file', required=True)
    p.add_argument('--windows', default='')
    p.add_argument('--peaks', action='store_true')
    p.add_argument('--min-reps-consensus', type=int, default=2)
    p.add_argument('--peak-width', type=int, default=100)
    p.add_argument('--peak-qvalue', type=float, default=0.01)
    p.add_argument('--effective-genome-size', type=int, default=2913022398)
    args = p.parse_args()
    if len(set(args.samples)) != len(args.samples) or not all(len(items) == len(args.samples) for items in (args.counts, args.headers, args.molecules)):
        p.error('Samples and input lists must be unique and have matching lengths')
    widths = sorted(set(int(value) for value in args.windows.split(',') if value))
    if any(width < 1 for width in widths) or args.peak_width < 2 or args.peak_width % 2 or not 0 < args.peak_qvalue < 1 or args.effective_genome_size < 1:
        p.error('Widths must be positive; peak-width must be even; qvalue must be between 0 and 1')
    design = json.loads(Path(args.design_file).read_text())
    metadata = {row['sample']: row for row in design['sample_metadata']}
    if set(metadata) != set(args.samples):
        p.error('Design metadata and analysis samples differ')
    groups = defaultdict(list)
    for sample in args.samples:
        groups[metadata[sample]['group']].append(sample)
    if args.peaks:
        if shutil.which('macs3') is None:
            p.error('MACS3 is required only when peak calling is enabled; activate an environment containing macs3')
        if any(not 1 <= args.min_reps_consensus <= len(samples) for samples in groups.values()):
            p.error('min-reps-consensus must not exceed the replicate count of any condition')
    header, sites = load_sites(args.counts, args.headers)
    report = {'design': design['design'], 'model_formula': design['model_formula'],
              'status': 'exploratory region discovery and counting; no differential model fitted',
              'windows': {},
              'samples': {sample: {'molecules': sum(count for rows in sample_sites.values() for _, count in rows),
                                   'occupied_sites': sum(len(rows) for rows in sample_sites.values())}
                          for sample, sample_sites in zip(args.samples, sites, strict=True)},
              'conditions': {}}
    with Path('analysis.samples.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=('sample', 'group', 'replicate', 'donor'), delimiter='\t')
        writer.writeheader()
        writer.writerows(metadata[sample] for sample in args.samples)
    for width in widths:
        regions = window_regions(sites, header, width)
        write_matrix(f'windows_{width}', regions, args.samples, sites)
        report['windows'][str(width)] = {'regions': len(regions), 'unobserved_bins': 'omitted', 'coordinate_system': '0-based half-open'}
    all_peaks = {}
    version = None
    if args.peaks:
        version = subprocess.check_output(['macs3', '--version'], text=True).strip()
        for sample, molecule_path, sample_sites in zip(args.samples, args.molecules, sites, strict=True):
            expected = sum(count for rows in sample_sites.values() for _, count in rows)
            observed = defaultdict(int)
            with Path(molecule_path).open() as handle:
                for line in handle:
                    fields = line.split()
                    if len(fields) != 6 or fields[5] not in ('+', '-') or int(fields[2]) != int(fields[1]) + 1:
                        raise ValueError(f'Invalid molecular BED for {sample}')
                    observed[(fields[0], int(fields[1]))] += 1
            if dict(observed) != {(chrom, start): count for chrom, rows in sample_sites.items() for start, count in rows}:
                raise ValueError(f'Molecular BED and counts differ for {sample}')
            folder = Path('SingleReplicate') / sample / 'peaks'
            folder.mkdir(parents=True)
            peakfile = folder / f'{sample}_peaks.narrowPeak'
            command = ['macs3', 'callpeak', '-t', str(molecule_path), '-f', 'BED', '-g', str(args.effective_genome_size),
                       '-n', sample, '--outdir', str(folder), '--nomodel', '--shift', str(-args.peak_width // 2),
                       '--extsize', str(args.peak_width), '--keep-dup', 'all', '-q', str(args.peak_qvalue),
                       '--min-length', str(args.peak_width), '--max-gap', str(args.peak_width)]
            if expected:
                with (folder / 'macs3.log').open('x') as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            else:
                peakfile.touch(exist_ok=False)
            peaks = []
            with peakfile.open() as handle:
                for line in handle:
                    fields = line.split()
                    peaks.append((fields[0], int(fields[1]), int(fields[2])))
            # Validate output coordinates before constructing condition consensus.
            union_regions(peaks, header)
            all_peaks[sample] = peaks
            report['samples'][sample]['peaks'] = len(peaks)
            write_json(folder / 'provenance.json', {'version': version, 'command': command, 'molecules': expected})
        common = []
        for group, samples in sorted(groups.items()):
            consensus = consensus_regions({sample: all_peaks[sample] for sample in samples}, args.min_reps_consensus)
            folder = Path('MergedReplicate') / group / 'peaks'
            folder.mkdir(parents=True)
            with (folder / f'{group}.consensus.bed').open('x') as handle:
                for i, (chrom, start, end, support) in enumerate(consensus, 1):
                    handle.write(f'{chrom}\t{start}\t{end}\t{group}_{i}\t{len(support)}\t{",".join(support)}\n')
                    common.append((chrom, start, end))
            report['conditions'][group] = {'replicates': len(samples), 'consensus_segments': len(consensus), 'minimum_replicates': args.min_reps_consensus}
        regions = union_regions(common, header)
        write_matrix('peaks_consensus', regions, args.samples, sites)
        report['peak_parameters'] = {'width': args.peak_width, 'qvalue': args.peak_qvalue, 'effective_genome_size': args.effective_genome_size, 'macs3_version': version, 'common_regions': len(regions)}
    write_json('analysis.summary.json', report)
    write_json('analysis_mqc.json', {'id': 'helixbusters_regions', 'section_name': 'Helixbusters exploratory regions',
                                   'description': report['status'], 'plot_type': 'table',
                                   'pconfig': {'id': 'helixbusters_regions_table', 'title': 'Region discovery'},
                                   'data': {**{f'windows_{width}': {'regions': info['regions']} for width, info in report['windows'].items()},
                                            **{f'consensus_{group}': {'regions': info['consensus_segments'], 'replicates': info['replicates']} for group, info in report['conditions'].items()}}})


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
