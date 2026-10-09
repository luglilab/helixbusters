"""Sparse molecule counts and replicate-supported hotspot intervals."""
from bisect import bisect_left
from collections import defaultdict
import csv
import json
from pathlib import Path

from helixbusters.reporting import bed_rows


def load_sites(paths, headers):
    """Validate reference compatibility and preserve integer molecule counts."""
    reference = json.loads(Path(headers[0]).read_text())
    result = []
    for path, header in zip(paths, headers, strict=True):
        if json.loads(Path(header).read_text()) != reference:
            raise ValueError('Samples use incompatible reference headers')
        sites = defaultdict(list)
        with Path(path).open() as handle:
            for _, chrom, start, end, count in bed_rows(handle, reference):
                sites[chrom].append((start, count))
        result.append(dict(sites))
    return reference, result


def window_regions(sites, header, width):
    """Return observed bins only; excluded/unobserved genome bins are omitted."""
    if width < 1:
        raise ValueError('Window width must be positive')
    occupied = defaultdict(set)
    for sample in sites:
        for chrom, rows in sample.items():
            occupied[chrom].update(start // width * width for start, _ in rows)
    return [(chrom, start, min(start + width, length)) for chrom, length in header
            for start in sorted(occupied[chrom])]


def consensus_regions(peaks, minimum):
    """Exact overlap support by distinct sample, avoiding transitive unions."""
    if not 1 <= minimum <= len(peaks):
        raise ValueError('Consensus threshold exceeds the number of biological replicates')
    events = defaultdict(lambda: defaultdict(list))
    for sample, intervals in peaks.items():
        for chrom, start, end in intervals:
            if start < 0 or end <= start:
                raise ValueError('Invalid peak interval')
            events[chrom][start].append((sample, 1))
            events[chrom][end].append((sample, -1))
    result = []
    for chrom in sorted(events):
        active = defaultdict(int)
        previous = None
        for position in sorted(events[chrom]):
            support = tuple(sorted(sample for sample, count in active.items() if count > 0))
            if previous is not None and previous < position and len(support) >= minimum:
                if result and result[-1][0] == chrom and result[-1][2] == previous and result[-1][3] == support:
                    result[-1] = (chrom, result[-1][1], position, support)
                else:
                    result.append((chrom, previous, position, support))
            for sample, delta in events[chrom][position]:
                active[sample] += delta
            previous = position
    return result


def union_regions(intervals, header):
    """Disjoint common counting universe across independent conditions."""
    grouped = defaultdict(list)
    lengths = dict(header)
    for chrom, start, end in intervals:
        if chrom not in lengths or not 0 <= start < end <= lengths[chrom]:
            raise ValueError('Peak lies outside the reference')
        grouped[chrom].append((start, end))
    result = []
    for chrom, _ in header:
        for start, end in sorted(grouped[chrom]):
            if result and result[-1][0] == chrom and start <= result[-1][2]:
                result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
            else:
                result.append((chrom, start, end))
    return result


def write_matrix(prefix, regions, samples, sites):
    """Count exact molecular ends in half-open intervals, without CPM rounding."""
    indexed = []
    for sample in sites:
        chroms = {}
        for chrom, rows in sample.items():
            cumulative = [0]
            for _, count in rows:
                cumulative.append(cumulative[-1] + count)
            chroms[chrom] = ([position for position, _ in rows], cumulative)
        indexed.append(chroms)
    with Path(f'{prefix}.counts.tsv').open('x') as matrix, Path(f'{prefix}.regions.bed').open('x') as bed:
        writer = csv.writer(matrix, delimiter='\t')
        writer.writerow(['region', 'chrom', 'start', 'end', *samples])
        for i, (chrom, start, end) in enumerate(regions, 1):
            name = f'region_{i}'
            counts = []
            for sample in indexed:
                positions, cumulative = sample.get(chrom, ([], [0]))
                counts.append(cumulative[bisect_left(positions, end)] - cumulative[bisect_left(positions, start)])
            writer.writerow([name, chrom, start, end, *counts])
            bed.write(f'{chrom}\t{start}\t{end}\t{name}\n')
