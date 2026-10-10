#!/usr/bin/env python3
"""Annotate both full-library prefix branches with identical explicit metadata."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.design import validate_design
from helixbusters.genomes import load_reference_config, GenomeFilter


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sensitivity-dir', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--gtf', type=Path, required=True)
    p.add_argument('--reference-config', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--genome', default='hg38')
    p.add_argument('--windows', default='10000,50000,100000')
    p.add_argument('--differential', action='store_true', help='Requires a working DESeq2 environment')
    p.add_argument('--contrast', nargs=2, default=['CHRONIC', 'ACUTE'])
    args = p.parse_args()
    if args.outdir.exists():
        p.error('Use a new output directory')
    source = args.sensitivity_dir.resolve()
    pilot = json.loads((source / 'pilot.json').read_text())
    if not pilot['parameters'].get('all_reads') or pilot['parameters']['genome'] != args.genome:
        p.error('Require complete-library sensitivity results with matching genome')
    with args.metadata.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'sample', 'group', 'replicate', 'donor'}.issubset(reader.fieldnames or []):
            p.error('Metadata requires sample/group/replicate/donor columns')
        rows = list(reader)
    samples = [row['sample'] for row in rows]
    if len(set(samples)) != len(samples) or set(samples) != set(pilot['samples']):
        p.error('Metadata sample set differs from the sensitivity cohort')
    design = validate_design(rows, 'paired')
    reference = load_reference_config(args.reference_config, args.genome, 'bwa')
    masking = GenomeFilter(args.genome, reference['blacklist_bed'], reference['blacklist_genome'])
    for sample in samples:
        for branch in ('baseline', 'candidate_trim'):
            src = source / sample / branch
            qc = json.loads((src / f'{sample}.mapping.json').read_text())
            if qc['genome_filter']['blacklist_sha256'] != masking.sha256:
                p.error(f'Blacklist differs from mapping for {sample}/{branch}')
            for suffix in ('counts.bed', 'families.tsv', 'q20.bam'):
                if not (src / f'{sample}.{suffix}').is_file():
                    p.error(f'Missing input for {sample}/{branch}: {suffix}')
    if not args.gtf.is_file():
        p.error('GTF does not exist')
    if args.differential:
        subprocess.run(['Rscript', '--vanilla', '-e', 'library(DESeq2)'], check=True)
    args.outdir.mkdir(parents=True)
    output = args.outdir.resolve()
    design_file = output / 'design.summary.json'
    design_file.write_text(json.dumps(design, indent=2) + '\n')
    scripts = Path(__file__).resolve().parent
    for branch in ('baseline', 'candidate_trim'):
        folder = output / branch; folder.mkdir()
        counts, headers, families = [], [], []
        for sample in samples:
            src = source / sample / branch
            bed, family, bam = src / f'{sample}.counts.bed', src / f'{sample}.families.tsv', src / f'{sample}.q20.bam'
            if not all(path.is_file() for path in (bed, family, bam)):
                raise ValueError(f'Missing sensitivity inputs for {sample}/{branch}')
            header = folder / f'{sample}.header.json'
            with pysam.AlignmentFile(bam, 'rb') as handle:
                header.write_text(json.dumps(list(zip(handle.references, handle.lengths)), indent=2) + '\n')
            counts.append(str(bed)); families.append(str(family)); headers.append(str(header))
        shared = ['--samples', *samples, '--counts', *counts, '--headers', *headers, '--design-file', str(design_file)]
        print(f'{branch}: windows and GENCODE annotation', flush=True)
        subprocess.run([sys.executable, str(scripts / 'analyze_regions.py'), *shared,
                        '--molecules', *families, '--windows', args.windows], cwd=folder, check=True)
        annotation = [sys.executable, str(scripts / 'annotate_regions.py'), *shared,
                      '--gtf', str(args.gtf.resolve()), '--gtf-genome', args.genome, '--genome', args.genome,
                      '--analysis-dir', str(folder), '--outdir', str(folder / 'Annotation')]
        if reference.get('blacklist_bed'):
            annotation += ['--blacklist', reference['blacklist_bed']]
        subprocess.run(annotation, cwd=folder, check=True)
        if args.differential:
            subprocess.run([sys.executable, str(scripts / 'differential_dsb.py'), '--analysis-dir', str(folder),
                            '--outdir', str(folder / 'Differential'), '--design', 'paired', '--contrast', *args.contrast], check=True)
    (output / 'sensitivity.analysis.json').write_text(json.dumps({
        'parameters': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        'reference': reference, 'design': design,
        'interpretation': 'Matched exploratory sensitivity branches. No MACS discovery in this comparison. Branch-specific abundance filters can yield different tested sets; compare shared tested features and effect sizes, not counts of significant genes alone.'}, indent=2) + '\n')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
