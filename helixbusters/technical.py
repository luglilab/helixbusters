"""Anchored, bounded edit distance for technical-prefix exclusion."""

T7_REVERSE = 'CCCTATAGTGAGTCGTATTA'
T7_FORWARD = 'TAATACGACTCACTATAGGG'


def full_t7_rescue_reason(sequence, quality):
    """Return exclusion reason, or None for the experimentally tested exact20 rule."""
    if not sequence.startswith(T7_REVERSE):
        return 'not_exact_full20_prefix'
    tail, tail_quality = sequence[20:], quality[20:]
    if len(tail) < 40:
        return 'short_tail'
    if 'N' in tail[:5] or min(ord(q) - 33 for q in tail_quality[:5]) < 20:
        return 'low_quality_boundary'
    if any(motif[:12] in tail[:60] for motif in (T7_FORWARD, T7_REVERSE)):
        return 'residual_T7_in_either_orientation'
    return None


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
