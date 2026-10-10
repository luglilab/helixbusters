#!/usr/bin/env python3
"""Compare relative DSB signal using explicit biological metadata and contrasts."""
import argparse
from pathlib import Path
import subprocess
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.differential import run_comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, help='Explicit donor metadata; sample/group/replicate must match source analysis')
    parser.add_argument('--design', choices=('paired', 'unpaired'), required=True)
    parser.add_argument('--contrast', nargs=2, metavar=('NUMERATOR', 'DENOMINATOR'), required=True)
    parser.add_argument('--minimum-count', type=int, default=5)
    parser.add_argument('--minimum-samples', type=int, default=2)
    parser.add_argument('--fdr', type=float, default=.05)
    parser.add_argument('--iterations', type=int, default=50)
    parser.add_argument('--seed', type=int, default=1729)
    args = parser.parse_args()
    metadata_path = args.metadata or args.analysis_dir / 'analysis.samples.tsv'
    metadata = pd.read_csv(metadata_path, sep='\t', keep_default_na=False, dtype=str)
    if not {'sample', 'group', 'replicate', 'donor'}.issubset(metadata):
        parser.error('Metadata requires sample, group, replicate and donor columns')
    summary = run_comparison(args.analysis_dir, args.outdir, metadata, args.design, *args.contrast,
                             args.minimum_count, args.minimum_samples, args.fdr, args.iterations, args.seed)
    print(f'Completed {len(summary["families"])} feature families: {args.outdir}')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
