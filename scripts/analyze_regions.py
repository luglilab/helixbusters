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
from helixbusters.peaks import call_sample_peaks, read_peaks
from helixbusters.peak_overlap import write_peak_overlap


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for option in ('samples', 'counts', 'headers', 'molecules'):
        p.add_argument(f'--{option}', nargs='+', required=True)
    p.add_argument('--design-file', required=True)
    p.add_argument('--peak-files', nargs='+')
    p.add_argument('--peak-provenance', nargs='+')
    p.add_argument('--windows', default='')
    p.add_argument('--peaks', action='store_true')
    p.add_argument('--nolambda', action='store_true', help='Use the global MACS3 background instead of local lambda')
    p.add_argument('--min-reps-consensus', type=int, default=2)
    p.add_argument('--peak-width', type=int, default=100)
    p.add_argument('--peak-qvalue', type=float, default=0.01)
    p.add_argument('--effective-genome-size', type=int, default=2913022398)
    args = p.parse_args()
    if len(set(args.samples)) != len(args.samples) or not all(len(items) == len(args.samples) for items in (args.counts, args.headers, args.molecules)):
        p.error('Samples and input lists must be unique and have matching lengths')
    if args.peak_files or args.peak_provenance:
        if not args.peaks or not args.peak_files or not args.peak_provenance or len(args.peak_files) != len(args.samples) or len(args.peak_provenance) != len(args.samples):
            p.error('External peak files and provenance must match all samples with --peaks')
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
        if not args.peak_files and shutil.which('macs3') is None:
            p.error('MACS3 is required only when peak calling is enabled; activate an environment containing macs3')
        if any(not 1 <= args.min_reps_consensus <= len(samples) for samples in groups.values()):
            p.error('min-reps-consensus must not exceed the replicate count of any condition')
    header, sites = load_sites(args.counts, args.headers)
    report = {'design': design['design'], 'model_formula': design['model_formula'],
              'analysis_purpose': design.get('analysis_purpose', 'biological'),
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
        versions = set()
        for index, (sample, molecule_path, sample_sites) in enumerate(zip(args.samples, args.molecules, sites, strict=True)):
            if args.peak_files:
                peaks = read_peaks(args.peak_files[index], header)
                provenance = json.loads(Path(args.peak_provenance[index]).read_text())
                expected = report['samples'][sample]['molecules']
                if provenance.get('sample') != sample or provenance['molecules'] != expected or any(provenance.get(key) != value for key, value in
                        [('width', args.peak_width), ('qvalue', args.peak_qvalue),
                         ('effective_genome_size', args.effective_genome_size), ('nolambda', args.nolambda)]):
                    raise ValueError(f'Peak provenance differs from analysis inputs for {sample}')
            else:
                peaks, provenance = call_sample_peaks(sample, molecule_path, sample_sites, header,
                    args.peak_width, args.peak_qvalue, args.effective_genome_size, args.nolambda,
                    Path('SingleReplicate') / sample / 'peaks')
            all_peaks[sample] = peaks
            report['samples'][sample]['peaks'] = len(peaks)
            report['samples'][sample]['macs3_version'] = provenance['version']
            versions.add(provenance['version'])
        common = []
        condition_intervals = {}
        for group, samples in sorted(groups.items()):
            consensus = consensus_regions({sample: all_peaks[sample] for sample in samples}, args.min_reps_consensus)
            condition_intervals[group] = [(chrom, start, end) for chrom, start, end, _ in consensus]
            folder = Path('MergedReplicate') / group / 'peaks'
            folder.mkdir(parents=True)
            with (folder / f'{group}.consensus.bed').open('x') as handle:
                for i, (chrom, start, end, support) in enumerate(consensus, 1):
                    handle.write(f'{chrom}\t{start}\t{end}\t{group}_{i}\t{len(support)}\t{",".join(support)}\n')
                    common.append((chrom, start, end))
            report['conditions'][group] = {'replicates': len(samples), 'consensus_segments': len(consensus), 'minimum_replicates': args.min_reps_consensus}
        regions = union_regions(common, header)
        write_matrix('peaks_consensus', regions, args.samples, sites)
        report['condition_overlap'] = write_peak_overlap(condition_intervals)
        version = next(iter(versions)) if len(versions) == 1 else sorted(versions)
        report['peak_parameters'] = {'width': args.peak_width, 'qvalue': args.peak_qvalue, 'effective_genome_size': args.effective_genome_size, 'nolambda': args.nolambda, 'background': 'global' if args.nolambda else 'local', 'macs3_version': version, 'common_regions': len(regions)}
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
