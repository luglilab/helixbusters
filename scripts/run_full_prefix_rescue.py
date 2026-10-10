#!/usr/bin/env python3
"""Experimental complete-library exact20 rescue with joint UMI deduplication."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.technical import T7_REVERSE, prefix_distance, full_t7_rescue_reason
from helixbusters.design import validate_design
from helixbusters.genomes import load_reference_config
from helixbusters.mapping import map_sample, validate_index
from helixbusters.deduplication import deduplicate_bam
from inspect_fastq_layout import read_manifest, reverse_complement
from inspect_full_prefix_boundary import diagnose_bam
from pilot_aligner_comparison import accepted_reads, compare_reads


def prepare_joint(source, sample, barcode, orientation, folder):
    """One raw FASTQ scan; production subset and its exact20 superset remain separate."""
    observed = barcode if orientation == 'forward' else reverse_complement(barcode)
    counts = Counter(dict.fromkeys(('input_reads', 'barcode_mismatch_reads', 'invalid_umi_reads',
                                   'technical_prefix_reads', 'short_insert_reads', 'output_reads'), 0))
    exclusions = Counter(); rescued = 0; identifiers = set()
    opener = gzip.open if str(source).lower().endswith('.gz') else open
    end = 8 + len(barcode)
    paths = {b: folder / f'{b}.fastq.gz' for b in ('baseline', 'candidate_trim')}
    with opener(source, 'rt', encoding='ascii') as raw, gzip.open(paths['baseline'], 'xt') as baseline, \
            gzip.open(paths['candidate_trim'], 'xt') as combined, gzip.open(folder / 'read_provenance.tsv.gz', 'xt') as provenance:
        writer = csv.writer(provenance, delimiter='\t'); writer.writerow(['read', 'source_class', 'removed_insert_bases'])
        while True:
            header = raw.readline()
            if not header:
                break
            seq = raw.readline().rstrip('\r\n').upper()
            plus = raw.readline(); qual = raw.readline().rstrip('\r\n')
            counts['input_reads'] += 1
            if not header.startswith('@') or not plus.startswith('+') or not seq or len(seq) != len(qual) or set(seq) - set('ACGTN') or any(not 33 <= ord(q) <= 126 for q in qual):
                raise ValueError(f'Malformed FASTQ record {counts["input_reads"]} in {sample}')
            if counts['input_reads'] % 1000000 == 0:
                print(f'{sample}: scanned {counts["input_reads"]:,} original reads', flush=True)
            if seq[8:end] != observed:
                counts['barcode_mismatch_reads'] += 1
                continue
            umi = seq[:8]
            if len(umi) != 8 or set(umi) - set('ACGT'):
                counts['invalid_umi_reads'] += 1
                continue
            insert, quality = seq[end:], qual[end:]
            trim = 0
            if prefix_distance(insert, T7_REVERSE[:18], 2) <= 2:
                counts['technical_prefix_reads'] += 1
                reason = full_t7_rescue_reason(insert, quality)
                if reason:
                    exclusions[reason] += 1
                    continue
                trim = 20
                rescued += 1
            elif len(insert) < 20:
                counts['short_insert_reads'] += 1
                continue
            else:
                counts['output_reads'] += 1
            identifier, *description = header.rstrip('\r\n').split(maxsplit=1)
            name = identifier + '_' + umi
            if name in identifiers:
                raise ValueError(f'Duplicate prepared read identifier in {sample}')
            identifiers.add(name)
            full_header = name + (' ' + description[0] if description else '')
            if not trim:
                baseline.write(f'{full_header}\n{insert}\n+\n{quality}\n')
            combined.write(f'{full_header}\n{insert[trim:]}\n+\n{quality[trim:]}\n')
            writer.writerow([name[1:], 'exact20_rescue' if trim else 'production_retained', trim])
    if sum(v for k, v in counts.items() if k != 'input_reads') != counts['input_reads']:
        raise ValueError('Exclusive preparation counts do not conserve original reads')
    if rescued + sum(exclusions.values()) != counts['technical_prefix_reads']:
        raise ValueError('Technical-prefix exclusions do not conserve reads')
    result = {'production_counts': dict(counts), 'rescue_exclusions': dict(exclusions), 'rescued_reads': rescued,
              'baseline_input_reads': counts['output_reads'], 'combined_input_reads': counts['output_reads'] + rescued,
              'parameters': {'barcode': barcode, 'observed_barcode': observed, 'barcode_orientation': orientation,
                             'umi_length': 8, 'technical_exclusion_motif': T7_REVERSE[:18], 'technical_prefix_max_errors': 2,
                             'minimum_production_insert_length': 20, 'rescue_rule': 'exact20; tail>=40; next5 Q>=20, no N; no either-orientation T7 prefix12 within first60 tail bases'}}
    (folder / 'preparation.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--reference-config', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--barcode-orientation', choices=('forward', 'reverse_complement'), required=True)
    p.add_argument('--genome', default='hg38')
    p.add_argument('--gtf', type=Path, help='Enable matched window, annotation and PCA analysis')
    p.add_argument('--differential', action='store_true')
    p.add_argument('--contrast', nargs=2, default=['CHRONIC', 'ACUTE'])
    p.add_argument('--map-threads', type=int, default=6)
    p.add_argument('--sort-threads', type=int, default=1)
    args = p.parse_args()
    if args.outdir.exists() or args.map_threads < 1 or args.sort_threads < 0:
        p.error('Use a new directory and valid CPU limits')
    manifest = read_manifest(args.manifest.resolve())
    for row in manifest:
        if row['sample'] in {'.', '..'} or any(x in row['sample'] for x in '/\\\t\r\n'):
            p.error('Unsafe sample name')
        if not Path(row['fastq']).is_file():
            p.error(f"Missing FASTQ: {row['fastq']}")
    if len({row['fastq'] for row in manifest}) != len(manifest):
        p.error('Distinct demultiplexed FASTQs are required')
    with args.metadata.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'sample', 'group', 'replicate', 'donor'}.issubset(reader.fieldnames or []):
            p.error('Explicit sample/group/replicate/donor metadata required')
        rows = list(reader)
    by_sample = {row['sample']: row for row in rows}
    if len(by_sample) != len(rows) or set(by_sample) != {row['sample'] for row in manifest}:
        p.error('Metadata sample set differs from manifest')
    for row in manifest:
        if any(row[key] != by_sample[row['sample']][key] for key in ('group', 'replicate')):
            p.error(f"Group/replicate metadata differs for {row['sample']}")
    design = validate_design(rows, 'paired')
    reference = load_reference_config(args.reference_config, args.genome, 'bwa')
    validate_index(reference['genome_index'], 'bwa')
    if args.gtf and not args.gtf.is_file():
        p.error('GTF does not exist')
    if args.differential:
        if not args.gtf:
            p.error('--differential requires --gtf')
        subprocess.run(['Rscript', '--vanilla', '-e', 'library(DESeq2)'], check=True)
    args.outdir.mkdir(parents=True)
    report = {'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              'reference': reference, 'design': design, 'samples': {},
              'interpretation': 'Complete original FASTQs; fresh production baseline versus production plus exact20 rescue, jointly deduplicated. No barcode mismatch rescue and no 12–17-base clipping correction. These are preprocessing versions of the same biological samples, not additional biological replicates. Experimental rescue is not independent physical DSB validation.'}
    report['parameters']['all_reads'] = True
    metrics = []
    for row in manifest:
        sample = row['sample']; folder = args.outdir / sample; folder.mkdir()
        print(f'{sample}: one complete original-FASTQ scan; exact20 sensitivity preparation', flush=True)
        prep = prepare_joint(row['fastq'], sample, row['barcode'], args.barcode_orientation, folder)
        item = {'manifest': row, 'preparation': prep, 'mapping': {}}
        endpoints = {}
        for branch, expected in (('baseline', prep['baseline_input_reads']), ('candidate_trim', prep['combined_input_reads'])):
            if not expected:
                raise ValueError(f'Empty {branch} for {sample}; inspect preparation.json')
            print(f'{sample}: mapping {branch}; {expected:,} prepared reads', flush=True)
            dest = folder / branch
            mapped = map_sample(sample, folder / f'{branch}.fastq.gz', reference['genome_index'], dest,
                                aligner='bwa', min_mapq=20, threads=args.map_threads, sort_threads=args.sort_threads,
                                genome=reference['genome'], blacklist_bed=reference['blacklist_bed'], blacklist_genome=reference['blacklist_genome'])
            qc = json.loads(Path(mapped['MappingQC']).read_text())
            if qc['primary_records'] != expected:
                raise ValueError('Prepared cohort and mapping record count differ')
            dedup = deduplicate_bam(mapped['BamFilteredPath'], dest / f'{sample}.families.tsv', dest / f'{sample}.counts.bed',
                                    method='directional', umi_length=8, min_mapq=20, five_prime_policy='strict', output_qc=dest / f'{sample}.dedup.json')
            with gzip.open(dest / 'strict_boundary_reads.tsv.gz', 'xt') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['sample', 'branch', 'read', 'chrom', 'endpoint', 'strand', 'mapq', 'NM', 'prefix', *[f'mismatch_{i}' for i in range(1, 11)]])
                diagnosis = diagnose_bam(Path(mapped['BamFilteredPath']), sample, branch, writer)
            item['mapping'][branch] = {'mapping': qc, 'deduplication': dedup, 'boundary_diagnostics': diagnosis}
            endpoints[branch] = accepted_reads(mapped['BamFilteredPath'])
            metrics.append({'sample': sample, 'branch': branch, 'original_reads': prep['production_counts']['input_reads'],
                            'prepared_reads': expected, 'exact20_prepared_rescue_reads': prep['rescued_reads'],
                            'strict_reads': dedup['accepted_reads'], 'molecules': dedup['deduplicated_molecules']})
        item['endpoint_comparison'] = compare_reads(endpoints['baseline'], endpoints['candidate_trim'])
        newly_accepted = endpoints['candidate_trim'].keys() - endpoints['baseline'].keys()
        source_classes = Counter()
        with gzip.open(folder / 'read_provenance.tsv.gz', 'rt') as handle:
            for record in csv.DictReader(handle, delimiter='\t'):
                if record['read'] in newly_accepted:
                    source_classes[record['source_class']] += 1
        item['newly_accepted_source_classes'] = dict(source_classes)
        del endpoints
        report['samples'][sample] = item
        (args.outdir / 'pilot.json').write_text(json.dumps(report, indent=2) + '\n')
    with (args.outdir / 'full_rescue.metrics.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]), delimiter='\t'); writer.writeheader(); writer.writerows(metrics)
    if args.gtf:
        scripts = Path(__file__).resolve().parent
        cmd = [sys.executable, str(scripts / 'analyze_prefix_sensitivity.py'), '--sensitivity-dir', str(args.outdir.resolve()),
               '--metadata', str(args.metadata.resolve()), '--gtf', str(args.gtf.resolve()), '--reference-config', str(args.reference_config.resolve()),
               '--genome', args.genome, '--outdir', str(args.outdir.resolve() / 'Analysis'), '--contrast', *args.contrast]
        if args.differential:
            cmd += ['--differential']
        subprocess.run(cmd, check=True)
        for branch in ('baseline', 'candidate_trim'):
            folder = args.outdir.resolve() / 'Analysis' / branch
            subprocess.run([sys.executable, str(scripts / 'pca_windows.py'), '--analysis-dir', str(folder), '--outdir', str(folder / 'PCA')], check=True)
    print(f'Complete: {args.outdir}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, EOFError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
