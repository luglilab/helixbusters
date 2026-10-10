#!/usr/bin/env python3
"""Audit uniform original-read cohorts and pilot exact full-prefix removal."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.technical import prefix_distance
from helixbusters.genomes import load_reference_config
from helixbusters.mapping import map_sample, validate_index
from helixbusters.deduplication import deduplicate_bam
from inspect_fastq_layout import read_manifest, reverse_complement
from pilot_aligner_comparison import reservoir, accepted_reads, compare_reads


def audit_cohort(source, output, barcode, known_barcodes, motif, umi_length=8):
    """Keep production exclusion order; never rescue or reassign barcodes."""
    outcomes = Counter(dict.fromkeys(('barcode_mismatch', 'invalid_umi', 'technical_prefix',
                                     'short_insert', 'production_retained'), 0))
    diagnostics, prefixes, tail_prefixes, barcodes, candidates = Counter(), Counter(), Counter(), Counter(), Counter()
    barcode_quality = Counter()
    quality_sums = Counter()
    offsets = Counter()
    end = umi_length + len(barcode)
    known = sorted(set(known_barcodes))
    paths = {branch: output / f'{branch}.fastq.gz' for branch in
             ('production_retained', 'excluded_exact_untrimmed', 'excluded_exact_trim')}
    handles = {branch: gzip.open(path, 'xt') for branch, path in paths.items()}
    identifiers = set()
    try:
        with gzip.open(source, 'rt', encoding='ascii') as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline().rstrip('\r\n').upper()
                plus = handle.readline()
                quality = handle.readline().rstrip('\r\n')
                if not header.startswith('@') or not plus.startswith('+') or not sequence or len(sequence) != len(quality) or set(sequence) - set('ACGTN'):
                    raise ValueError('Malformed original FASTQ cohort')
                if any(not 33 <= ord(q) <= 126 for q in quality):
                    raise ValueError('Invalid Phred+33 quality')
                name = header.split()[0]
                if name in identifiers:
                    raise ValueError('Duplicate read identifier in cohort')
                identifiers.add(name)
                observed = sequence[umi_length:end]
                barcodes[observed] += 1
                q = [ord(x) - 33 for x in quality[umi_length:end]]
                status = 'exact' if observed == barcode else 'mismatch'
                if len(q) == len(barcode):
                    barcode_quality[status] += 1
                    quality_sums[status] += sum(q) / len(q)
                for offset in range(0, 25):
                    if sequence[offset:offset + len(barcode)] == barcode:
                        offsets[offset] += 1
                if observed != barcode:
                    outcomes['barcode_mismatch'] += 1
                    if len(observed) != len(barcode):
                        diagnostics['barcode_too_short'] += 1
                    else:
                        distance = sum(a != b for a, b in zip(observed, barcode))
                        diagnostics[f'expected_barcode_hamming_{distance}'] += 1
                        distances = {b: sum(a != c for a, c in zip(observed, b)) for b in known}
                        nearest = min(distances.values())
                        winners = [b for b in known if distances[b] == nearest]
                        diagnostics[f'manifest_nearest_distance_{nearest}'] += 1
                        if len(winners) > 1:
                            diagnostics['manifest_nearest_tied'] += 1
                        elif winners[0] == barcode and nearest == 1:
                            diagnostics['expected_unique_nearest_one_mismatch_DIAGNOSTIC_ONLY'] += 1
                        elif nearest == 0:
                            diagnostics['exact_other_manifest_barcode'] += 1
                    continue
                umi = sequence[:umi_length]
                if len(umi) != umi_length or set(umi) - set('ACGT'):
                    outcomes['invalid_umi'] += 1
                    continue
                insert, insert_quality = sequence[end:], quality[end:]
                identifier = name + '_' + umi
                def write(branch, seq, qual):
                    handles[branch].write(f'{identifier}\n{seq}\n+\n{qual}\n')
                if prefix_distance(insert, motif, 2) <= 2:
                    outcomes['technical_prefix'] += 1
                    prefixes[insert[:60]] += 1
                    if not insert.startswith(motif):
                        candidates['edited_prefix_boundary_uncertain'] += 1
                        continue
                    tail, tail_quality = insert[len(motif):], insert_quality[len(motif):]
                    tail_prefixes[tail[:40]] += 1
                    if len(tail) < 40:
                        candidates['exact_prefix_short_tail'] += 1
                    elif 'N' in tail[:5] or min(ord(x) - 33 for x in tail_quality[:5]) < 20:
                        candidates['exact_prefix_low_quality_boundary'] += 1
                    elif motif[:12] in tail[:60]:
                        candidates['exact_prefix_repeated_motif_first60'] += 1
                    else:
                        candidates['exact_single_prefix_pilot'] += 1
                        write('excluded_exact_untrimmed', insert, insert_quality)
                        write('excluded_exact_trim', tail, tail_quality)
                elif len(insert) < 20:
                    outcomes['short_insert'] += 1
                else:
                    outcomes['production_retained'] += 1
                    write('production_retained', insert, insert_quality)
    finally:
        for handle in handles.values():
            handle.close()
    if sum(candidates.values()) != outcomes['technical_prefix']:
        raise ValueError('Technical-prefix classification does not conserve reads')
    return {'sampled_reads': sum(outcomes.values()), 'exclusive_outcomes': dict(outcomes),
            'barcode_diagnostics': dict(diagnostics), 'expected_barcode_offsets': dict(offsets),
            'mean_barcode_Q': {k: quality_sums[k] / n for k, n in barcode_quality.items()},
            'top_observed_barcodes': barcodes.most_common(20), 'technical_prefix_classes': dict(candidates),
            'top_excluded_insert_prefixes': prefixes.most_common(30), 'top_exact_prefix_tails': tail_prefixes.most_common(30)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--barcode-orientation', choices=('forward', 'reverse_complement'), required=True)
    p.add_argument('--max-reads', type=int, default=100000)
    p.add_argument('--seed', type=int, default=1729)
    p.add_argument('--motif', default='CCCTATAGTGAGTCGTAT')
    p.add_argument('--reference-config', type=Path, help='Enable BWA pilot; otherwise audit only')
    p.add_argument('--genome', default='hg38')
    p.add_argument('--map-threads', type=int, default=6)
    p.add_argument('--sort-threads', type=int, default=1)
    args = p.parse_args()
    if args.outdir.exists() or args.max_reads < 1 or args.map_threads < 1 or args.sort_threads < 0:
        p.error('Use a new directory and positive read/CPU limits')
    if len(args.motif) < 12 or set(args.motif) - set('ACGT'):
        p.error('Motif requires at least 12 A/C/G/T bases')
    manifest = read_manifest(args.manifest.resolve())
    for row in manifest:
        if not row['sample'] or row['sample'] in {'.', '..'} or any(x in row['sample'] for x in '/\\\t\r\n'):
            p.error('Unsafe sample identifier')
        if not Path(row['fastq']).is_file():
            p.error(f"Missing FASTQ: {row['fastq']}")
    if len({row['fastq'] for row in manifest}) != len(manifest):
        p.error('This demultiplexed-library audit requires distinct FASTQs per sample')
    if len({len(row['barcode']) for row in manifest}) != 1:
        p.error('Barcode diagnostics require equal barcode lengths')
    oriented = lambda b: b if args.barcode_orientation == 'forward' else reverse_complement(b)
    known = [oriented(row['barcode']) for row in manifest]
    reference = None
    if args.reference_config:
        reference = load_reference_config(args.reference_config, args.genome, 'bwa')
        validate_index(reference['genome_index'], 'bwa')
    args.outdir.mkdir(parents=True)
    report = {'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
              'sampling': 'Uniform original-FASTQ reservoir; one complete sequential scan per library',
              'interpretation': 'No barcode rescue. Manifest barcode catalogue may omit other libraries from the sequencing run. Exact full-prefix removal is experimental; mapping and strict acceptance do not establish physical DSB identity. Repeated/edited prefixes are excluded from the correction pilot.',
              'reference': reference, 'python': sys.version, 'samples': {}}
    metrics = []
    for row in manifest:
        sample = row['sample']; folder = args.outdir / sample; folder.mkdir()
        print(f'{sample}: uniform original FASTQ sampling (full file scan)', flush=True)
        sampling = reservoir(row['fastq'], folder / 'original.sampled.fastq.gz', args.max_reads, args.seed)
        audit = audit_cohort(folder / 'original.sampled.fastq.gz', folder, oriented(row['barcode']), known, args.motif)
        audit['sampling'] = sampling
        if audit['sampled_reads'] != sampling['sampled_reads']:
            raise ValueError('Audit does not conserve the raw cohort')
        audit['mapping'] = {}; endpoints = {}
        sizes = {'production_retained': audit['exclusive_outcomes']['production_retained'],
                 'excluded_exact_untrimmed': audit['technical_prefix_classes'].get('exact_single_prefix_pilot', 0),
                 'excluded_exact_trim': audit['technical_prefix_classes'].get('exact_single_prefix_pilot', 0)}
        for branch, size in sizes.items():
            if not reference or not size:
                metrics.append({'sample': sample, 'branch': branch, 'raw_sampled_reads': audit['sampled_reads'], 'branch_input_reads': size, 'strict_reads': None, 'molecules': None})
                continue
            print(f'{sample}: mapping {branch}; {size:,} reads', flush=True)
            dest = folder / branch
            mapped = map_sample(sample, folder / f'{branch}.fastq.gz', reference['genome_index'], dest,
                                aligner='bwa', min_mapq=20, threads=args.map_threads, sort_threads=args.sort_threads,
                                genome=reference['genome'], blacklist_bed=reference['blacklist_bed'], blacklist_genome=reference['blacklist_genome'])
            qc = json.loads(Path(mapped['MappingQC']).read_text())
            if qc['primary_records'] != size:
                raise ValueError('Mapping primary count differs from input cohort')
            dedup = deduplicate_bam(mapped['BamFilteredPath'], dest / f'{sample}.families.tsv', dest / f'{sample}.counts.bed',
                                    method='directional', umi_length=8, min_mapq=20, five_prime_policy='strict', output_qc=dest / f'{sample}.dedup.json')
            audit['mapping'][branch] = {'mapping': qc, 'deduplication': dedup}
            if branch.startswith('excluded_exact'):
                endpoints[branch] = accepted_reads(mapped['BamFilteredPath'])
            metrics.append({'sample': sample, 'branch': branch, 'raw_sampled_reads': audit['sampled_reads'], 'branch_input_reads': size,
                            'strict_reads': dedup['accepted_reads'], 'molecules': dedup['deduplicated_molecules']})
        if len(endpoints) == 2:
            audit['excluded_endpoint_comparison'] = compare_reads(endpoints['excluded_exact_untrimmed'], endpoints['excluded_exact_trim'])
        (folder / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')
        report['samples'][sample] = audit
        (args.outdir / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    with (args.outdir / 'audit.metrics.tsv').open('x') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics[0]), delimiter='\t'); writer.writeheader(); writer.writerows(metrics)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, EOFError, UnicodeError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
