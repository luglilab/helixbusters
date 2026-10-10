"""Exact condition membership of replicate-supported consensus intervals."""
from collections import Counter
import csv
from pathlib import Path

from helixbusters.regions import consensus_regions


def write_peak_overlap(condition_intervals, outdir='PeakOverlap'):
    """Report disjoint half-open segments, not transitive overlap clusters."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    folder = Path(outdir)
    folder.mkdir(parents=True, exist_ok=False)
    names = sorted(condition_intervals)
    segments = consensus_regions(condition_intervals, 1)
    counts, bases = Counter(), Counter()
    with (folder / 'condition_membership.bed').open('x') as handle:
        for i, (chrom, start, end, members) in enumerate(segments, 1):
            handle.write(f'{chrom}\t{start}\t{end}\tsegment_{i}\t{",".join(members)}\n')
            counts[members] += 1
            bases[members] += end - start
    # Show all three combinations for a two-condition comparison, even if empty.
    combinations = [(name,) for name in names]
    if len(names) > 1:
        combinations.append(tuple(names))
    combinations = sorted(set(combinations) | set(counts), key=lambda x: (len(x), x))
    with (folder / 'intersections.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['conditions', 'segments', 'covered_bp'])
        writer.writerows((','.join(c), counts[c], bases[c]) for c in combinations)
    with (folder / 'set_sizes.tsv').open('x') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['condition', 'segments', 'covered_bp'])
        writer.writerows((name, sum(n for c,n in counts.items() if name in c),
                         sum(n for c,n in bases.items() if name in c)) for name in names)
    fig, (bars, dots) = plt.subplots(2, 1, figsize=(max(5, len(combinations)*1.2), 4.5),
                                    gridspec_kw={'height_ratios':[3,1]}, sharex=True, layout='constrained')
    bars.bar(range(len(combinations)), [counts[c] for c in combinations], color='#4879ae')
    bars.set_ylabel('Disjoint consensus segments')
    bars.set_title('Condition overlap (exploratory)')
    for x,c in enumerate(combinations):
        bars.text(x, counts[c], str(counts[c]), ha='center', va='bottom')
        dots.scatter([x]*len(names), range(len(names)), color='#dddddd', s=35)
        selected = [i for i,name in enumerate(names) if name in c]
        dots.plot([x]*len(selected), selected, 'o-', color='black')
    dots.set_yticks(range(len(names)), names); dots.set_xticks([])
    dots.set_xlim(-.6, max(.6,len(combinations)-.4)); dots.set_ylim(-.5,max(.5,len(names)-.5))
    fig.savefig(folder / 'conditions_upset.png', dpi=180)
    fig.savefig(folder / 'conditions_upset.pdf'); plt.close(fig)
    (folder / 'README.md').write_text(
        'Membership uses exact genomic overlap of condition consensus intervals. '\
        'Disjoint segments are split when condition membership changes; touching intervals do not overlap. '\
        'Segment counts are not original peak counts or independent statistical tests. '\
        'Covered bases provide a fragmentation-resistant summary. Unique means detected only in that '\
        'condition at these parameters, not significantly differential.\n')
    return {'segments':len(segments), 'covered_bp':sum(bases.values()),
            'unit':'disjoint segments of exact condition membership',
            'interpretation':'exploratory detection overlap; not differential enrichment'}
