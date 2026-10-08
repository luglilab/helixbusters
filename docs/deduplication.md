# BLISS coordinates and UMI deduplication

This stage accepts **one biological library per coordinate-sorted BAM**, with
UMIs already extracted. It does not change FASTQ quality filtering or mapping.
Deduplication is independent of group names, treatment, and comparison design.
Do not pool biological replicates before this stage. Multiple `SM` values in
the BAM read-group header are rejected; unlabelled pooled libraries cannot be
detected automatically.

## Coordinate contract

Outputs use zero-based, half-open BED intervals `[position, position + 1)`.
The position is the terminal aligned base at the biological 5-prime end:

| Alignment | Position |
| --- | --- |
| Forward | `reference_start` |
| Reverse | `reference_end - 1` |

For example, an alignment spanning `[100, 120)` yields `[100, 101)` on `+`
and `[119, 120)` on `-`. This represents a base adjacent to a labelled end,
not an interbase cut coordinate or a count of complete double-strand breaks.
Internal insertions/deletions are handled using the reference span from CIGAR.
Chromosome and contig names are preserved exactly, including X, Y, MT and
alternative contigs. Any biological contig exclusion must be explicit upstream.

The default `five_prime_policy="strict"` requires a terminal `M`, `=` or `X`
at the biological 5-prime CIGAR end. Reads with 5-prime soft/hard clipping or
terminal insertions/deletions are excluded and counted in QC. Clipping at the
3-prime end is allowed. `five_prime_policy="aligned"` explicitly permits using
the alignment boundary for an uncertain 5-prime end; it never extrapolates
clipped bases into the reference. Both policies exclude spliced CIGARs (`N`).
These policies do not reconstruct sequence removed by upstream 5-prime trimming.

Unmapped, secondary, supplementary and QC-failed records are excluded. MAPQ 255
(unavailable) is excluded; other reads are filtered by `min_mapq` (default 0,
because the existing pipeline supplies an already filtered BAM). **This is not
a Q30 base-quality filter.** For an unfiltered BAM, choose `min_mapq` explicitly.
Coordinate-only duplicate flags are retained so that distinct UMIs are not lost.

Single-end is the default. A paired-end BAM requires an explicit
`read_selection="read1"` or `"read2"` according to the end bearing the BLISS label.
Only that mate is counted. This is not paired-fragment deduplication and does not
repair the existing paired-end FASTQ preprocessing implementation.

## UMI contract and grouping

The `RX` tag takes precedence when present. Otherwise the UMI is read from the
last underscore-delimited suffix of the read name, allowing underscores in the
original name. UMIs must contain only uppercase A/C/G/T and have `umi_length`
bases (default 8). Missing UMIs, ambiguous bases, and invalid lengths are skipped
with separate QC counts. If all otherwise eligible records have unusable UMIs,
the run fails instead of publishing an apparently successful empty result.

Grouping always stays within **chromosome + position + strand**:

- `method="exact"` (default): each distinct UMI is one inferred molecule.
- `method="directional"`: connect UMI A to UMI B when their Hamming distance is
  at most `max_distance` (default 1) and `count(A) >= 2 * count(B) - 1`. Traverse
  directed paths from roots in decreasing abundance, assigning each UMI once.
  Ties use lexical ordering so output is deterministic. The representative UMI
  is the root and its support is the sum of read counts in the family.

The Python implementation uses the published UMI-tools directional rule. It
does not invoke `umi_tools dedup`, and is not claimed to reproduce its full BAM
behavior: clipping conventions and tie handling differ. Directional paths are
transitive, so the root and a distant descendant may differ by more than the
edge threshold. Neighboring singletons can also merge under this rule. Separate
biological molecules with related UMIs can therefore be merged; validation on
real data is still required. `max_distance=0` reduces to exact grouping; 2 is
available for sensitivity analysis, not the default for 8-base UMIs.

There is **no grouping across nearby genomic positions**, including the 8-nt
neighborhood described in the original BLISS paper. Exact and directional counts
are a baseline for evaluating that additional biological assumption later.

## Usage

After mapping with the existing class:

```python
helixbusters.generate_umi_output_for_samples()  # exact baseline
# Or choose correction explicitly:
helixbusters.generate_umi_output_for_samples(method="directional", umi_length=8)
```

These calls use the same per-sample filenames and replace successful previous
outputs. To retain side-by-side comparisons, use separate output directories or
the independent function with distinct paths:

```python
from pathlib import Path
from helixbusters.deduplication import deduplicate_bam

for method in ("exact", "directional"):
    output = Path("results") / method
    qc = deduplicate_bam(
        "sample.q20.bam",
        output / "families.tsv",
        output / "counts.bed",
        method=method,
        umi_length=8,
        output_molecules=output / "molecules.bed",
        output_sites=output / "sites.tsv",
        output_qc=output / "deduplication.json",
    )
    print(method, qc["deduplicated_molecules"])
```

No BAM index is required. Actual coordinate order is checked, even if the header
says `SO:coordinate`. Counts are buffered until later alignment starts cannot
contribute to an earlier end position; this handles reverse reads of different
lengths without retaining the entire BAM. Memory depends on local end density
and alignment spans. Output is ordered by BAM reference order, position, strand,
then family abundance and lexical UMI. Outputs are staged in temporary files
and replaced only after successful processing; publication of several files is
not a single filesystem transaction.

## Outputs and migration notes

| Output | Columns / meaning |
| --- | --- |
| Existing `*_Chromosome-Location-Strand-UMI-PCR.txt` | `chrom start end strand representative_umi read_support`, no header |
| Existing `*_Chromosome-Location-UMI-Count.bed` | `chrom start end molecules`, no header, summed across strands |
| New `*_sites.tsv` | `chrom start end strand molecules reads`, with header |
| New `*_molecules.bed` | BED6, one row per inferred molecule: `chrom start end representative_umi 0 strand` |
| New `*_deduplication.json` | Accepted/skipped records, exact and corrected molecule counts, inferred duplicate reads, parameters |

The legacy `PCR` filename is retained for API compatibility. Its last column is
**observed read support**, not a measured number of PCR cycles or exclusively
PCR-derived duplicates; BLISS also involves other amplification steps.

Previous negative-strand coordinates were incorrect. X/Y names were changed and
non-numeric contigs dropped; those behaviors are removed. The old count BED could
contain indistinguishable rows from opposite strands: it now has one summed row
per position. Use `sites.tsv` when strand is required.

BED4 count values are not molecule multiplicities understood by MACS3. The new
BED6 has an actual record per inferred molecule; downstream tools must preserve
distinct molecules at identical coordinates (for MACS3, `--keep-dup all`).
This change supplies counts and molecule intervals, not a deduplicated BAM.

QC `duplicate_reads` means `accepted_reads - deduplicated_molecules`;
`umi_groups_merged` means `exact_molecules - deduplicated_molecules`.
Each rejected read receives its first applicable exclusion reason, so
`input_reads = accepted_reads + sum(skipped.values())` for successful runs.

## Validation

The tests create real, unindexed BAMs with known coordinates and UMIs. They cover
both strands, coordinate zero, indels, clipping, contig preservation, filtering,
mate selection, RX precedence, malformed UMIs, directional chains, balanced
UMIs, deterministic ties, sample separation, and failure without replacing
existing outputs. Randomized small graphs are compared with an independent
all-pairs directional-graph oracle.

From the repository root, in a Python 3.12 environment with pandas and pysam:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
```

Matplotlib is needed only when generating alignment quality plots. These tests
establish software behavior, not biological validation on the example FASTQ:
that requires a reference genome and mapped reads.

## References

- [pysam coordinates and alignment API](https://pysam.readthedocs.io/en/latest/api.html)
- [UMI-tools grouping methods](https://umi-tools.readthedocs.io/en/stable/the_methods.html)
- [Original BLISS paper](https://www.nature.com/articles/ncomms15058)
