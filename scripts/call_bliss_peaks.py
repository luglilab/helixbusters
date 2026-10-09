#!/usr/bin/env python3
"""Run one independently deduplicated BLISS library through MACS3."""
import argparse
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from helixbusters.peaks import call_sample_peaks
from helixbusters.regions import load_sites
from helixbusters.reporting import write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('sample', 'counts', 'header', 'molecules'):
        parser.add_argument(f'--{name}', required=True)
    parser.add_argument('--peak-width', type=int, default=100)
    parser.add_argument('--peak-qvalue', type=float, default=0.01)
    parser.add_argument('--effective-genome-size', type=int, required=True)
    parser.add_argument('--nolambda', action='store_true')
    args = parser.parse_args()
    header, sites = load_sites([args.counts], [args.header])
    peaks, provenance = call_sample_peaks(args.sample, args.molecules, sites[0], header,
        args.peak_width, args.peak_qvalue, args.effective_genome_size, args.nolambda, '.')
    Path('provenance.json').rename(f'{args.sample}.provenance.json')
    write_json(f'{args.sample}.versions.json', {'MACS3_CALLPEAK': {'macs3': provenance['version']}})
    # Empty libraries still emit all declared process outputs.
    if not provenance['molecules']:
        Path(f'{args.sample}_peaks.xls').touch(exist_ok=False)
        Path(f'{args.sample}_summits.bed').touch(exist_ok=False)
        Path('macs3.log').write_text('Empty molecular library; peak calling skipped.\n')
    print(f'{args.sample}: {len(peaks)} peaks')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        sys.exit(1)
