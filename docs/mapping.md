# Mapping BLISS reads

Mapping consumes reads **after** UMI/barcode extraction. The UMI suffix in each
read name is retained for downstream deduplication. Plain FASTQ and `.fastq.gz`
inputs are passed directly to the aligner. The example raw FASTQ cannot be used
unchanged: its 8-base UMI and 8-base sample barcode must first be removed from
the sequence, retaining the UMI in the identifier.

## Existing API

To select hg19/hg38/mm10/mm39 and remove noncanonical, mitochondrial and
blacklisted reads, configure the [genome and blacklist](genomes.md) on the
Helixbusters instance. The examples below inherit that configuration.

```python
helixbusters.run_bwa_mapping(quality=20, threads=8, sort_threads=1, sort_memory="768M")

# Alternative, using a Bowtie2 index prefix:
helixbusters.run_bowtie2_mapping(
    quality=20, threads=8, sort_threads=1, bowtie2_mode="end-to-end"
)
```

`quality` is the minimum **alignment MAPQ**, not FASTQ Q30. The default remains
20; selecting it is not a conclusion that 20 is biologically optimal for BLISS.
MAPQ scales differ between aligners, so compare retained loci and biological
results as well as percentages when evaluating BWA against Bowtie2.

The class honors `PathReadForwardTrimmed` and `PathReadReverseTrimmed` when set.
Otherwise it uses the existing output filenames (`trimmed.fastq.gz` for SE or
`trimmed_R1.fastq.gz` / `trimmed_R2.fastq.gz` for PE). All sample FASTQ paths are
checked before mapping the first sample. Missing inputs now fail explicitly.
Group names and treatment comparisons do not affect mapping.

Classic BWA-MEM uses its existing alignment settings and may soft-clip reads.
Bowtie2 defaults explicitly to `end-to-end`; `bowtie2_mode="local"` is available
for a controlled comparison. No extra ATAC/Tn5 coordinate shifts are applied.

## Independent per-sample API

```python
from helixbusters.mapping import map_sample

outputs = map_sample(
    sample="Sample1",
    read1="Sample1.trimmed.fastq.gz",
    genome_index="references/human_bwa_index",
    output_dir="results/Sample1/bwa_q20",
    aligner="bwa",
    min_mapq=20,
    threads=8,
    sort_threads=1,
    sort_memory="768M",
)
```

Provide `read2` for paired-end data. This stage aligns and filters both mates;
subsequent deduplication still requires explicit `read_selection="read1"` or
`"read2"` according to the library design. This does not fix the earlier
paired-end FASTQ extraction implementation.

Both APIs accept `aligner_executable` and `samtools_executable` for explicit
executable paths; otherwise programs are found on PATH. Index prefixes must
point to a complete classic BWA index (`.amb/.ann/.bwt/.pac/.sa`) or a complete
Bowtie2 `.bt2` or `.bt2l` set. This does not support BWA-MEM2 index formats.
Index existence is checked; an index is not built or scientifically validated
against the requested genome build automatically.

`threads` controls the aligner. Sorting runs concurrently; `sort_threads` means
samtools' additional threads, and `sort_memory` is its approximate memory limit
per thread. Account for both programs when allocating CPUs and memory. The
pipeline does not interpret `threads` as a strict total process-tree CPU cap.

## Filtering and outputs

| Output | Content |
| --- | --- |
| `sample.all.bam` and `.bai` | Coordinate-sorted complete aligner output, including unmapped and non-primary records |
| `sample.q20.bam` and `.bai` | Filtered alignments; number in the name is the configured MAPQ threshold |
| `sample.mapping.json` | Counts, exclusions, MAPQ histogram, inputs, parameters, command vectors and version information |
| `sample.bwa.log` or `sample.bowtie2.log` | Aligner stderr for the latest attempt |
| `sample.sort.log` | samtools sort stderr for the latest attempt |

The class records these paths in `BamAllPath`, `BamFilteredPath`, `MappingQC`,
`MappingLog` and `SortLog`. BWA and Bowtie2 use the same BAM names for backward
compatibility. Use separate output directories to retain side-by-side aligner
or parameter comparisons; a successful rerun replaces the corresponding files.

The filtered BAM excludes unmapped, secondary, supplementary, QC-failed, unknown
MAPQ (255), and below-threshold alignments. It retains coordinate-only duplicate
flags, both selected input mates, and 5-prime clipping. It does not require
proper-pair flags. When a genome is selected, canonical nuclear and blacklist
filters are also applied as described in [genomes.md](genomes.md). The strict 5-prime policy in the
[deduplication stage](deduplication.md) makes the later decision about ambiguous
ends. MAPQ zero is retained only when the threshold permits it; MAPQ alone is
not proof that an alignment is unique or a genomic end biologically correct.

QC counts are **alignment records**, not molecules, fragments, or DSBs.
`primary_records` includes unmapped primary records; `mapped_primary_records`
excludes unmapped records but includes QC-failed records. Their MAPQ histogram
is calculated before filtering. `retained_ambiguous_five_prime` counts retained
records rejected by the strict end-coordinate policy. Each excluded alignment
has one reason, so `alignment_records = retained_records + sum(excluded_by_reason)`.

## Execution and failure behavior

Commands are argument vectors with no shell interpolation. Both the aligner's
and sorter's exit codes are checked. Failure or interruption stops and reaps
the remaining child processes. Aligner and sorting errors remain in log files.

All BAMs, BAM indexes and the JSON report are prepared in a temporary directory
on the output filesystem. They replace final outputs only after mapping,
filtering, and both index operations succeed. A processing failure preserves
previous successful BAMs and reports. Log files always describe the latest
attempt, including a failed one. Replacement of multiple files is not a single
filesystem transaction, so an OS error during publication can leave mixed files.
The command history includes temporary sorting paths that are removed afterward.

## Tests and scope of validation

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -v
```

Tests cover command arguments, complete index detection, primary/MAPQ filtering,
QC accounting, failures on either side of the pipe, preservation of successful
outputs after failures, and core integration. Executable fixtures test the BWA
command path without claiming to validate the real BWA alignment algorithm.

When `bowtie2`, `bowtie2-build` and `samtools` are on PATH, an additional real
integration test creates a small deterministic genome and gzipped FASTQ, builds
the index, aligns both strands, sorts, filters, indexes and deduplicates. Its
expected result is three molecules at known end coordinates. This establishes
execution and coordinate behavior on synthetic data, not the best biological
mapper, MAPQ threshold, or FASTQ quality policy for an experimental dataset.

## References

- [BWA manual](https://bio-bwa.sourceforge.net/bwa.shtml)
- [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools sorting](https://www.htslib.org/doc/samtools-sort.html)
- [SAM/BAM specification, including MAPQ](https://samtools.github.io/hts-specs/SAMv1.pdf)
