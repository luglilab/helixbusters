# Genome selection and exclusion regions

Helixbusters supports these explicit assembly selections:

| Selection | Accepted assembly alias | Species | Canonical nuclear chromosomes |
| --- | --- | --- | --- |
| `hg19` | `GRCh37` | human | 1–22, X, Y |
| `hg38` | `GRCh38` | human | 1–22, X, Y |
| `mm10` | `GRCm38` | mouse | 1–19, X, Y |
| `mm39` | `GRCm39` | mouse | 1–19, X, Y |

Aliases are case-insensitive. Select one name, e.g. `mm39` or `GRCm39`.
The assembly determines the species, canonical set, and expected chromosome
lengths. Both UCSC names (`chr1`, `chrX`) and Ensembl names (`1`, `X`) are
recognized for nuclear chromosomes. BAM names are preserved in outputs.

## Configure resources once, then choose the build

Copy [references.example.json](references.example.json) to a local reference
catalog and replace the placeholder paths with your **existing** BWA/Bowtie2
index prefixes and build-specific blacklist BEDs. Entries for unused builds
may be omitted. Paths can be absolute or relative to the catalog file.

```python
from helixbusters import Helixbusters

hb = Helixbusters.from_reference_config(
    samplesheet="samplesheet.csv",
    genome="hg38",                 # or hg19, mm10, mm39
    reference_config="references.json",
    aligner="bwa",                 # or bowtie2
)
hb.read_column_from_excel()
hb.create_sample_output_folders("results/hg38")
hb.process_infofile(umi_length=8, threads=8)
hb.run_bwa_mapping(quality=20, threads=8)
hb.generate_umi_output_for_samples(method="directional")
```

Choosing a different build selects that entry's index and blacklist; it does
not convert coordinates between assemblies. Use a separate results directory
for each build. When selecting `aligner="bowtie2"`, run `run_bowtie2_mapping()`;
using the other mapper with a catalog-selected index is rejected.

The example demonstrates connecting the stages, not validation of the existing
FASTQ extraction/trimming stage or an endorsement of a particular Q30 policy.
That earlier stage remains separate from the tested mapping and deduplication
changes. You can instead supply already trimmed FASTQs as described in
[mapping.md](mapping.md).

The catalog selects local resources. Helixbusters does **not** download full
genomes, build indexes, or fetch blacklists implicitly. The example paths are
placeholders, not installed resources. Only the selected entry must be usable.

Alternatively, configure one reference directly:

```python
hb = Helixbusters(
    samplesheet="samplesheet.csv",
    genome="GRCm39",
    genome_index="/references/mm39/bwa/genome",
    blacklist_bed="/references/mm39/mm39.validated.bed.gz",
    blacklist_genome="mm39",
)
```

The explicit `blacklist_genome` declaration is required because BED files do
not generally identify their assembly. A conflicting declaration is an error,
and a missing blacklist never silently disables filtering for a selected build.
The declaration is not proof that an incorrectly labelled BED came from the
claimed assembly: resource provenance must be checked when preparing the catalog.
If `species` is also supplied, it must agree with the build.

Old calls supplying only `species` and `genome_index`, without `genome`, retain
the previous mapping behavior with **no assembly/blacklist filter** for backward
compatibility. Use an explicit build or the reference catalog for this workflow.

## What is removed and when

Mapping runs against the supplied reference, allowing nuclear, mitochondrial,
and alternative loci to compete for alignment. After primary-alignment and MAPQ
filtering, the selected-build filter removes:

1. Mitochondrial alignments named `chrM`, `chrMT`, `MT`, or `M`.
2. Any alignment outside the selected nuclear canonical set: alternate loci,
   unplaced/unlocalized scaffolds, decoys, and other extra contigs.
3. Any read with **at least one aligned reference base** inside a blacklist
   interval, regardless of strand and regardless of whether its 5-prime end
   is inside the interval.

Blacklist overlap uses CIGAR alignment blocks and zero-based half-open BED
coordinates. Soft clips and reference gaps (`D`/`N`) do not by themselves count
as aligned read bases. For example, with blacklist `[100,200)`, a read ending
at 100 or starting at 200 does not overlap. A read spanning 99–101 does.

Filtering is per read, including in paired-end data. It does not discard the
other mate solely because its partner is excluded or rewrite pairing flags.
Deduplication must still explicitly select the BLISS-bearing mate. See
[deduplication.md](deduplication.md).

`BamFilteredPath` is the cleaned BAM consumed by deduplication, so subsequent
molecule and site outputs use only retained alignments. `BamAllPath` remains a
complete diagnostic BAM. The filtered BAM preserves the reference header,
including excluded contig definitions; those contigs have no retained reads.

## Compatibility checks and QC

The BAM must contain exactly one recognized name for every canonical nuclear
chromosome, with its expected assembly length. Missing chromosomes, duplicated
aliases such as both `1` and `chr1`, or incorrect lengths stop processing before
the filtered BAM is written. This deliberately excludes partial/chromosome-only
references and accession-only naming such as `NC_000001.11`; use a complete
reference with UCSC/Ensembl names. The check occurs after mapping, when the BAM
header is available. It detects common assembly mistakes but does not authenticate
reference sequence content, masking choices, or patch versions.

Lengths are from UCSC chromosome-size tables:
[hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes),
[hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes),
[mm10](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes),
[mm39](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes).
Mitochondrial lengths are not checked because mitochondrial alignments are
excluded and hg19/GRCh37 distributions can differ in that reference sequence.

Blacklist BED/BED.gz files are checked for valid integer coordinates,
nonnegative starts, positive widths, and canonical chromosome bounds. Canonical
intervals are merged and indexed for overlap queries. Noncanonical intervals
are counted and ignored because those reads are already excluded. A BED with
no applicable canonical intervals is rejected, including accession-name files
that would otherwise silently miss all reads.

The mapping QC JSON records `mitochondrial`, `noncanonical`, and `blacklist`
exclusion reasons, plus genome, canonical set, blacklist path, declared build,
SHA-256 of the original BED/BED.gz bytes, interval counts, and overlap policy.
Each rejected record receives only its first applicable reason. Earlier
unmapped/flag/MAPQ filters take precedence; these exclusion counters therefore
are not independent totals of every mitochondrial or blacklisted alignment.

## Blacklist sources

The [Boyle Lab blacklist repository](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)
provides assembly-specific lists for hg19, hg38 and mm10. These are reasonable
resources to evaluate for removing problematic genomic regions; their effect
on BLISS signal should still be inspected in the QC and biological analyses.

For GRCm39/mm39, use a dedicated mm39 exclusion set, such as the resource
described by [excluderanges](https://github.com/dozmorovlab/excluderanges) and its
[publication](https://doi.org/10.1093/bioinformatics/btad198). Do not relabel the
mm10 BED as mm39. The original Boyle repository does not supply an mm39 list.
An upstream [issue reports negative starts in an mm39 BED export](https://github.com/dozmorovlab/excluderanges/issues/7).
Such files fail validation here; no silent clipping, coordinate shifting or
automatic liftover is performed. Prepare a verified BED and retain its source
and any conversion steps alongside the reference catalog.

## Validation

Synthetic BAM tests cover all four builds, UCSC/Ensembl names, mitochondrial
aliases, alternative contigs, X/Y retention, half-open overlap boundaries,
CIGAR gaps/clips, malformed BEDs, wrong-build headers, missing blacklists and
catalog mismatches. Filtered output is passed through deduplication to check
that excluded reads never contribute molecules. No full human/mouse reference
mapping or real-data blacklist sensitivity analysis has been performed yet.
