"""Anchored, bounded edit distance for technical-prefix exclusion."""


def prefix_distance(sequence, motif, maximum=2):
    """Return edit distance to a prefix, capped at maximum + 1.

    Leading bases cannot be skipped freely. All substitutions, insertions and
    deletions cost one. The band bounds work for reads without the motif.
    """
    if sequence.startswith(motif):
        return 0
    sequence = sequence[:len(motif) + maximum]
    cap = maximum + 1
    previous = {j: j for j in range(min(len(sequence), maximum) + 1)}
    for i, base in enumerate(motif, 1):
        current = {}
        if i <= maximum:
            current[0] = i
        for j in range(max(1, i - maximum), min(len(sequence), i + maximum) + 1):
            current[j] = min(cap, previous.get(j, cap) + 1,
                             current.get(j - 1, cap) + 1,
                             previous.get(j - 1, cap) + (base != sequence[j - 1]))
        if not current or min(current.values()) > maximum:
            return cap
        previous = current
    return min(previous.values(), default=cap)
