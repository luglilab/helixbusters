"""Shared per-library MACS3 invocation on validated molecular ends."""
from collections import defaultdict
from pathlib import Path
import subprocess

from helixbusters.regions import union_regions
from helixbusters.reporting import validate_label, write_json


def read_peaks(path, header):
    peaks = []
    with Path(path).open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                raise ValueError('Invalid peak row')
            peaks.append((fields[0], int(fields[1]), int(fields[2])))
    union_regions(peaks, header)
    return peaks


def call_sample_peaks(sample, molecule_path, sample_sites, header, peak_width,
                      peak_qvalue, effective_genome_size, nolambda, folder):
    validate_label(sample)
    if peak_width < 2 or peak_width % 2 or not 0 < peak_qvalue < 1 or effective_genome_size < 1:
        raise ValueError('Invalid peak parameters')
    version = subprocess.check_output(['macs3', '--version'], text=True).strip()
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
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    peakfile = folder / f'{sample}_peaks.narrowPeak'
    protected = [peakfile, folder / f'{sample}_peaks.xls', folder / f'{sample}_summits.bed',
                 folder / 'macs3.log', folder / 'provenance.json',
                 folder / f'{sample}.provenance.json', folder / f'{sample}.versions.json']
    if any(path.exists() for path in protected):
        raise FileExistsError(f'Peak outputs already exist in {folder}; use a new output directory')
    command = ['macs3', 'callpeak', '-t', str(molecule_path), '-f', 'BED', '-g', str(effective_genome_size),
               '-n', sample, '--outdir', str(folder), '--nomodel', '--shift', str(-peak_width // 2),
               '--extsize', str(peak_width), '--keep-dup', 'all', '-q', str(peak_qvalue),
               '--min-length', str(peak_width), '--max-gap', str(peak_width)]
    if nolambda:
        command.append('--nolambda')
    if expected:
        with (folder / 'macs3.log').open('x') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    else:
        peakfile.touch(exist_ok=False)
    peaks = read_peaks(peakfile, header)
    provenance = {'sample': sample, 'version': version, 'command': command, 'molecules': expected,
                  'width': peak_width, 'qvalue': peak_qvalue,
                  'effective_genome_size': effective_genome_size, 'nolambda': nolambda}
    write_json(folder / 'provenance.json', provenance)
    return peaks, provenance
