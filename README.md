![Helixbusters](https://github.com/luglilab/helixbusters/blob/master/logo.png)

# Helixbusters

Helixbusters analyzes BLISS-like sequencing data. The Nextflow DSL2 workflow
coordinates UMI/barcode extraction, mapping, alignment filtering, UMI-based
deduplication, bigWig tracks and MultiQC reports for samples and conditions;
the Python package implements mapping, deduplication and reporting.

## Install

From a workstation or HPC login node with Conda/Mamba:

```bash
git clone https://github.com/luglilab/helixbusters.git
cd helixbusters
conda env create --file environment.yml
conda activate helixbusters
python -m pip install --no-deps -e .
module load nextflow/26.04.6  # on the HPC
nextflow -version
```

The environment contains Python 3.12, Excel/CSV readers, `pysam`, BWA,
Bowtie2, samtools, cutadapt, UMI-tools, deepTools, pyBigWig and MultiQC.
It has no machine-specific Conda
prefix. Nextflow is loaded separately through the HPC module system; see
[HPC installation](docs/hpc.md).

## Start a Nextflow run

The workflow currently supports **single-end sequencing**. Sequencing layout
is separate from whether biological samples are paired across conditions.
The Excel/CSV samplesheet requires `Sample`, `Replicate`, `Group`,
`PathReadForward`, and `SampleBarcodeForward`. Add `Donor` when donor identity
is known; it is required for paired designs.

Convert the workbook on a machine where the FASTQs are visible:

```bash
python scripts/samplesheet_to_manifest.py \
  /path/on/cluster/samplesheet_helixbuster.xlsx \
  /path/on/cluster/samples.tsv \
  --check-fastq
```

The TSV manifest uses `sample`, `replicate`, `group`, `fastq`, `barcode`, and
optional `donor`. Absolute HPC paths are preserved when converting on a
workstation; omit `--check-fastq` if those inputs are only accessible on the
cluster. The converter does not verify remote files. Use a new manifest filename
to preserve an existing metadata table.

`group` defines the condition to aggregate. `replicate` identifies a biological
replicate within that condition. One library per condition/replicate is supported;
technical replicate libraries are not currently supported. Replicate labels such
as R1 are never interpreted automatically as donor identities.

Prepare `references.json` with an existing reference index and matching
blacklist. For BWA, `scripts/prepare_igenome.py` can fetch iGenomes and the
Boyle-Lab blacklist and write the catalog for hg19, hg38 or mm10. See the
[Nextflow guide](docs/nextflow.md). Then submit a pilot:

```bash
nextflow run main.nf \
  --manifest samples.tsv \
  --genome hg38 \
  --reference_config references.json \
  --aligner bwa \
  --design unspecified \
  --barcode_orientation forward \
  --umi_length 8 \
  --map_threads 8 \
  --sort_threads 1 \
  --queue YOUR_PARTITION \
  --outdir results/pilot_hg38 \
  -profile slurm
```

Choose the correct build (`hg19`, `hg38`, `mm10`, or `mm39`) from the experiment
reference and replace `YOUR_PARTITION` with a permitted Slurm partition. The
example above uses `hg38` only as syntax. Paired-end input is
rejected. Set the barcode orientation and UMI length from the library design,
not from this example. See the full [Nextflow guide](docs/nextflow.md)
before running all samples.

## Experimental design

| Option | Meaning and validation |
|---|---|
| `--design paired` | Requires a `donor` for every sample, exactly once in every declared condition. Incomplete pairing is rejected. |
| `--design unpaired` | Declares independent biological samples. Donor metadata is optional, but declared donors cannot recur across conditions. |
| `--design unspecified` | Default for existing manifests; mapping and QC can proceed, but pairing must be resolved before statistical inference. |

To validate a paired samplesheet during conversion, pass `--design paired`
to `samplesheet_to_manifest.py`, then pass the same option to Nextflow.
Different treatment labels do not make measurements from the same donor
independent. The workflow records intended designs (`~ donor + condition` or
`~ condition`) and experimental units in `MultiQC/design.summary.json`,
`design.metadata.tsv`, and the MultiQC report. Differential fitting is opt-in.

## Verify the read layout

Inspect original FASTQs before choosing extraction parameters:

```bash
python scripts/inspect_fastq_layout.py \
  --manifest samples.tsv --max-reads 100000 --output fastq_layout.json
```

This read-only diagnostic checks barcode positions in both orientations, read
lengths, initial sequence composition and shared input paths. It samples the
first records of each file, not a random sample or a complete read count. See
[FASTQ layout diagnostics](docs/fastq_layout.md).

## Preparation, mapping and molecule counting

`EXTRACT_UMI` runs `scripts/prepare_bliss_reads.py` automatically. It requires
an exact barcode immediately after the initial UMI, removes both, and retains
the UMI in the read name. `--barcode_orientation forward` uses the manifest
sequence; `reverse_complement` transforms it. Avoid inverting barcodes that
are already listed in the observed orientation. UMIs containing N are excluded,
and `--minimum_insert_length` defaults to 20 bases after extraction.

Optional technical-prefix exclusion uses:

```bash
--technical_prefix CCCTATAGTGAGTCGTAT --technical_prefix_max_errors 2
```

This is an **assay-specific example**, not a universal BLISS adapter filter.
The default motif is empty, and the default error allowance is 0. Values 1 or 2
allow anchored substitutions, insertions and deletions. Matching reads are
excluded entirely; the filter does not trim or rescue them. It does not target
internal motif occurrences. Confirm the motif and compare molecule yield in a
pilot before applying it to a new dataset.

Preparation reports count original FASTQ reads and exclusive exclusion reasons.
The workflow stops if a sample retains no reads. Mapping applies the
selected build's canonical nuclear chromosome set, removes mitochondrial and
blacklist-overlapping reads, and excludes low/unknown MAPQ and nonprimary
alignments. `--mapq` is alignment MAPQ, not per-base Q30.

Directional UMI deduplication is performed per sample and genomic end/strand.
Strict 5-prime acceptance excludes reads with ambiguous biological ends,
including clipping at that end. Molecules are not deduplicated across donors.
See [genome configuration](docs/genomes.md),
[mapping](docs/mapping.md), and [deduplication](docs/deduplication.md).

Reference preparation is a separate step: the Nextflow run expects existing
indexes and blacklists for the same build.

## Outputs and interpretation

| Directory | Contents |
|---|---|
| `prepared/` | FASTQs after UMI/barcode removal and preparation filters |
| `SingleReplicate/<sample>/` | `mapping/`, `deduplication/`, `qc/` and `bigwig/` |
| `MergedReplicate/<group>/` | Condition-pooled filtered BAMs, QC and molecule tracks |
| `MultiQC/` | `multiqc_report.html`, `multiqc_data/`, design metadata and reporting tables |

Molecule-end tracks include raw counts and CPM. Condition tracks also include
the equal-weight mean of replicate CPM when all replicates have nonzero counts.
Coverage CPM tracks use filtered alignments **before UMI deduplication** and
are mapping diagnostics, distinct from molecule-end signal. The report separates
Samtools plots for SingleReplicate and MergedReplicate and hides noncanonical
contigs in chromosome plots while preserving the underlying QC files.

Use per-sample molecule counts for biological inference; pooled tracks are
descriptive. CPM and total recovered molecules do not establish breaks per cell
or absolute damage burden without suitable calibration and experimental metadata.
Residual technical sequence and 5-prime clipping must be assessed before biological
interpretation. Changing a filter should use a new output directory to retain
the previous analysis. Keep Nextflow work/cache directories for `-resume`.

## Focused diagnostic tools

| Script | Purpose |
|---|---|
| `inspect_five_prime.py` | Read-only diagnosis of 5-prime CIGAR operations and clipped sequences in a BAM |
| `inspect_technical_prefix.py` | Count anchored motif variants in prepared FASTQ |
| `pilot_technical_filter.py` | Compare baseline and two-error exclusion on the same prepared-read subset, with mapping and molecule counting |

These tools are separate from the main workflow; use
`python scripts/<name>.py --help` for arguments. Diagnostic reports refuse
existing output filenames.

## Downstream analysis status

Optional genomic windows and MACS3 hotspots are available in `main.nf`:

```bash
# Append to the usual Nextflow command; both options default to false.
--run_windows --window_sizes 1000,5000,10000 \
--run_peak_calling --min_reps_consensus 2 --peak_width 100 --peak_qvalue 0.01
```

Windows use raw integer molecular end counts, independently deduplicated per
sample. Only bins observed in at least one sample are exported; other samples
receive zero for these bins. Bins are anchored at coordinate zero and clipped
at chromosome ends. `Analysis/windows_<width>.counts.tsv` and `.regions.bed`
provide a common region universe and retain every sample as a separate column.
`Analysis/analysis.samples.tsv` records condition, biological replicate and donor.

Peak calling requires **MACS3 in the active worker environment**. The project
environment specification pins MACS3 3.0.5; existing environments require an
explicit update. Each sample's molecular
BED6 contains one record per UMI family; MACS3 uses `--nomodel --keep-dup all`
to preserve independent molecules at identical coordinates. A 100-bp smoothing
width uses shift -50, extension 100 and minimum peak length/maximum gap 100.
This is exploratory hotspot discovery without an experimental control; the
smoothed intervals are not individual break coordinates. Default effective
genome sizes are 2,913,022,398 for hg38 and 2,652,783,500 for mm10. Other
assemblies require `--effective_genome_size`; this value is an approximation
and can be overridden for the reference/read length used.

For a sensitivity comparison without a control, add `--peak_nolambda` to
forward `--nolambda` to MACS3 and use the global background instead of local
lambda. It defaults to false. This can expose regional background biases;
additional calls are not evidence of biological specificity. Keep the same
q-value, smoothing and consensus threshold and use a separate output directory.
The exact command and background choice are recorded in provenance and
`Analysis/analysis.summary.json`. See the
[MACS3 documentation](https://macs3-project.github.io/MACS/docs/callpeak.html).
For sparse datasets, `--window_sizes 1000,5000,10000,50000,100000` also exports
50- and 100-kb windows without changing the MACS3 smoothing width.

Per-sample peaks, logs and provenance appear in `SingleReplicate/<sample>/peaks`.
`MergedReplicate/<condition>/peaks/<condition>.consensus.bed` contains exact
segments supported by at least `--min_reps_consensus` distinct biological
replicates of that condition. Its six columns are chromosome, start, end,
identifier, support count and comma-separated supporting samples (the last
column is not a strand). Overlap chains do not count as support across the
whole union. A condition with fewer replicates than the threshold fails before
mapping; use threshold 1 explicitly for descriptive singleton analysis.
Conditions are discovered independently and are never used as MACS3 controls
for each other. Their consensus intervals form a disjoint common universe in
`Analysis/peaks_consensus.{regions.bed,counts.tsv}` for counting all samples.

`MACS3_CALLPEAK` runs once per sample, with its own cache, logs, provenance and
software-version JSON. Its input is the deduplicated molecular BED6 rather
than a read-coverage BAM; biological pairing does not change this input format.
Each task requests 1 CPU/8 GB/4 hours under the reporting defaults. Tasks may
run concurrently according to scheduler availability. The process declares
`environment.yml` for Nextflow-managed Conda (`-with-conda`); without that flag,
the existing activated environment is used. No container image is configured.
The downstream reporting task consumes staged peaks, builds independent
condition consensus and counts their common regions. Changing window sizes
alone does not change the per-sample peak-calling task command.
Region counts and summaries use one separate 1-CPU/8-GB reporting task.
MultiQC includes a region discovery table. These options can be added on a
resumed run with the same work directory; use a new output directory to
preserve previously published results.

MACS3 q-values describe enrichment under its background model, not differences
between conditions. Low retained molecule counts and residual technical
artifacts must be considered before interpreting hotspots biologically.

## Clipping and molecular-depth diagnostics

For read-only investigation of residual clipping and usable molecular depth,
run `scripts/audit_bliss_recovery.py --single-replicate-dir /path/to/SingleReplicate
--outdir /path/to/new/audit`. It scans each filtered BAM once and uniformly samples
up to 50,000 primary mapped reads (seed 1729). The output compares clipped and
strict-accepted reads for clipping lengths, base quality and exact technical-motif
segments. Family-table counts must conserve accepted reads and deduplicated
molecules. Observed-family thinning curves do not extrapolate library complexity
or re-infer directional UMI clusters.

`scripts/pilot_clipped_prefix.py` is a separate **experimental pilot**, not a
production preprocessing option. It remaps the same filtered-BAM cohort with
and without a narrowly defined correction: only biological 5-prime soft clips
of 12–17 bases exactly matching the prefix of `CCCTATAGTGAGTCGTAT`, clip mean
quality >=30, next five bases >=Q20, and remaining insert >=40 bases. It restores
sequencing orientation and preserves authoritative UMIs. Other reads are
unchanged; both branches retain MAPQ/blacklist filtering and strict endpoints.
Supply `--single-replicate-dir`, `--samples`, `--reference-config` and a new
`--outdir`. The JSON and metrics report molecule yield, unchanged-control
concordance and candidate endpoint shifts. Source aligned boundaries are not
independently validated DSB coordinates; a yield increase alone does not authorize
adopting this correction. The cohort is conditional on original filtered mapping,
so its yield cannot be interpreted as a whole-library gain.

## Differential relative DSB signal and robustness

Enable an explicit contrast using `--run_differential --design paired
--contrast CHRONIC,ACUTE` (or `--design unpaired` for independent biological
samples), together with windows and/or a GTF. The manifest must explicitly
identify donors for paired analysis. At least three biological samples per
contrast condition are required; technical libraries are unsupported.
The numerator is first: positive effects indicate CHRONIC > ACUTE in this
example. Other conditions are not pooled into either contrast group.

The optional module needs R and DESeq2 in the active environment. Check first:

```bash
Rscript --vanilla -e 'stopifnot(requireNamespace("DESeq2", quietly=TRUE)); packageVersion("DESeq2")'
```

If absent, add the optional dependencies to the existing environment using
`conda env update -n helixbusters -f environment.differential.yml` without
`--prune`. The workflow checks DESeq2 before upstream processing when enabled;
ordinary mapping/reporting runs have no new R dependency.

`Analysis/Differential/` contains one subdirectory for each window width and
for promoter/gene-body counts, when available. DESeq2 fits a negative-binomial
Wald model to **integer per-sample molecule counts**, with `~ donor + condition`
for paired samples or `~ condition` for independent samples. The abundance
filter is >=5 molecules in >=2 samples, regardless of condition; customize with
`--differential_min_count` and `--differential_min_samples`. Families with fewer
than 20 eligible features are explicitly skipped. Primary normalization uses
DESeq2 `poscounts` size factors estimated within each family; a separate fit
using total retained library molecules is a normalization sensitivity check.
No gene length normalization is applied to the model: the same feature's length
is constant across samples. Densities remain descriptive annotation columns.

- `results.tsv`: complete feature universe, raw counts, eligibility, unshrunk
  log2 fold change, standard error, p-value, BH FDR, dispersion, Cook's distance,
  convergence and robustness diagnostics.
- `significant.tsv` and `CONDITION.higher_relative_signal.tsv`: primary FDR
  <=0.05 results, with the corresponding effect direction for condition lists.
- `dispersion_fit.pdf`, `diagnostics.{pdf,png}`, size factors, model logs,
  warnings and `sessionInfo.txt`: inspect before biological interpretation.
- Equal-depth support frequencies: 50 seeded molecule subsamples without
  replacement to the smallest library, including molecules outside each
  feature family in a residual category. Support means >=2 molecules in >=2
  biological samples within a condition; frequencies are not p-values.
- For paired samples, donor-specific CPM log ratios and leave-one-donor-out
  median ratios: descriptive influence diagnostics, without model refitting.

Set `--differential_fdr`, `--robustness_iterations` and `--robustness_seed` when
needed. Missing p-values are not evidence of no effect. BH correction is **within
each feature family**, not across all widths and contexts. Declare a primary
family (for example 10-kb windows) before interpreting the others as sensitivity
analyses. MACS consensus regions are excluded from these tests because discovery
uses the same samples; combined gene counts are excluded as an overlapping
alternative to promoter/gene-body summaries. This module does not establish
absolute DSB burden per cell or enrichment over a matched genomic background.
Normalization assumptions and disagreement between the two fits remain material
limitations; three donor pairs provide limited statistical power.

To analyze an existing run without remapping:

```bash
python scripts/differential_dsb.py \
  --analysis-dir /path/to/completed/Analysis \
  --outdir /path/to/new/Differential \
  --metadata /path/to/confirmed_paired_metadata.tsv \
  --design paired --contrast CHRONIC ACUTE
```

Metadata must contain `sample`, `group`, `replicate`, `donor` and match all source
samples, conditions and replicate labels. It can supply previously absent donor
identities, but cannot silently change groups, replicate labels or already
declared donors. Original
inputs are read-only; an existing output directory is refused. The reusable
Slurm wrapper `scripts/run_dsb_comparison.sh` requests 1 CPU/8 GB/4 hours.
MultiQC receives the differential summary during full workflow execution;
standalone comparisons write the corresponding custom-content JSON for a later
MultiQC run.

## GTF annotation and DSB feature allocation

Supply a complete gene/transcript/exon GTF for the same genome assembly:
For the current hg38 experiments, use the pinned GENCODE human release 50
comprehensive annotation on reference chromosomes (`gencode.v50.annotation.gtf.gz`,
GRCh38.p14; [official release page](https://www.gencodegenes.org/human/)).

```bash
# Add these options to the existing pipeline command:
--gtf /path/to/hg38.annotation.gtf.gz --gtf_genome hg38 \
--promoter_upstream 2000 --promoter_downstream 500
```

Annotation can run with or without window analysis and peak calling. It uses
retained deduplicated 1-bp DSB ends, not read-coverage BAMs or window midpoints.
GTF coordinates are converted from 1-based inclusive to 0-based half-open;
standard `chr1`/`1` aliases are matched without rewriting gene IDs or output
chromosomes. The declared GTF genome must match the mapping genome; recognized
build labels in GTF comments, chromosome bounds and occupied-contig coverage
are checked. A build declaration and coordinate checks cannot independently
prove the provenance of a mislabeled GTF. Original inputs are preserved.

Each molecule is assigned once with priority **promoter > exon > intron >
intergenic**. Promoters include all transcript TSS, using the strand-aware
relative interval `[-upstream, downstream)`; a gene-span TSS is used when
transcript records are absent. Exons are the union across all isoforms and
overlapping genes; introns are gene-body bases outside any exon and promoter.
CDS and UTR bases remain within the exon category. Gene-body totals are exon
plus intron, with promoter bases excluded to preserve exclusive categories.
All biotypes in the supplied GTF contribute. Enhancers cannot be identified
from a gene GTF alone and are not inferred from intergenic DSBs.

Outputs under `Analysis/Annotation` include:

- Per-sample DSB site annotations with molecular counts, associated gene IDs,
  gene names and biotypes.
- Window and peak annotations with exact exclusive feature coverage in bp,
  a dominant feature and a mixed-feature flag. A wide interval can overlap
  several genes/features; gene association does not imply functional targeting.
- Molecule counts and percentages per sample; exon-plus-intron gene-body totals.
- Condition means of sample percentages, sample standard deviations and pooled
  counts. Condition plots show equal-weight means and individual biological
  sample points, avoiding dominance by deeply sequenced replicates. Empty
  libraries have undefined percentages and are excluded from percentage means.
- PNG/PDF feature-distribution plots, MultiQC plots for samples and conditions,
  and GTF checksum, parameters and software-version provenance.

Percentages describe the allocation of retained molecules. They are not
enrichment scores or evidence of a treatment effect. A suitable enrichment
background must account for genomic opportunity, blacklist exclusions and
mappability; no such background or significance test is fitted here. Use a new
output directory on resumed runs to preserve previous results.

## Final per-condition gene signal lists

When GTF annotation is enabled, `Analysis/Annotation/GeneSignal` also contains
independent per-condition gene rankings and replicate-supported candidate lists:

- `CONDITION.promoter.candidate_genes.tsv`: reproducible promoter-associated DSB signal.
- `CONDITION.gene_body.candidate_genes.tsv`: reproducible exon/intron DSB signal,
  excluding sites assigned to any promoter.
- `CONDITION.combined.candidate_genes.tsv`: an additional promoter-plus-body summary.
- Corresponding `ranked_genes.tsv` files preserve all genes with observed signal,
  including those below the support threshold. Wide `genes.CONTEXT.counts.tsv`
  and `genes.CONTEXT.CPM.tsv` files preserve all annotated genes, including zeros.

Each ranking includes gene ID, name, biotype, coordinates, uniquely assignable
length, integer molecule counts and CPM per sample, condition mean/median/SD,
replicate support and descriptive CPM per kb. CPM divides by all retained
deduplicated molecules in each sample, not only gene-assigned molecules.
Condition summaries give biological samples equal weight. Zero-depth samples
have undefined CPM and are excluded from normalized condition summaries.

By default, a candidate requires at least `--gene_min_molecules 2` in at least
`--min_reps_consensus` biological samples of that condition (default two).
Use `--gene_min_reps` to set an independent gene-support threshold. A group
with too few replicates has an empty candidate list; thresholds are never
silently relaxed. Rankings place supported candidates first, then sort by
median sample CPM, mean sample CPM and stable gene ID/chromosome tie-breaks.
The support filter is exploratory and does not establish significance.

Gene counting uses only genes overlapping the highest-priority site category:
promoter, then exon, then intron. A site with multiple candidate genes is
excluded from integer gene counts and recorded in `ambiguous_gene_sites.tsv`.
`gene_assignment.samples.tsv` audits unique, ambiguous and intergenic molecules;
their sum conserves the retained sample total. This avoids assigning the same
DSB to several genes. Promoter and body lists are alternatives to the combined
summary; they are not additional independent molecules.

### Feature-length-adjusted gene density

Alongside the original CPM rankings, each condition/context now has
`density_ranked_genes.tsv` and `density_candidate_genes.tsv`. These sort by
median sample **CPM per effective kb**, then mean density, with the same
replicate/molecule support thresholds. Short genes cannot bypass the support
filter.

`effective_assignable_bp` subtracts the union of blacklist intervals from the
uniquely assigned GTF feature partition and restricts lengths to canonical
chromosomes. Promoter and body lengths follow the same feature priority and
ambiguity rules as gene counting. `excluded_assignable_bp` records the difference
from `uniquely_assignable_bp`. Counts and original rankings are preserved;
`mean_CPM_per_kb` keeps its original unmasked denominator. New mean/median
`CPM_per_effective_kb` columns use the blacklist-adjusted denominator, including
sample-level densities.

The annotation task uses the mapping environment's resolved blacklist and
checks its checksum, refusing a modified blacklist or inputs containing
blacklisted/noncanonical DSB ends. A standalone invocation can pass
`--blacklist /path/to/build_matched.bed.gz`; without a mask, lengths remain
unadjusted and provenance records this explicitly. No additional packages are
needed. Mappability and read-span-dependent loss near blacklist boundaries
are not modeled. Effective length is a point-based descriptive opportunity,
not a fully calibrated callable genome or enrichment background.

For the same gene in ACUTE and CHRONIC, length is constant. The differential
model must use preserved integer counts, not length-divided densities.

These are descriptive DSB-associated gene candidates, not differentially
damaged genes or differentially expressed genes. Larger features have more
opportunity to accumulate ends. CPM per effective kb is blacklist-adjusted but
not corrected for mappability and does not establish enrichment.
Condition comparisons require a separate replicate-aware differential model,
confirmed donor structure, suitable normalization, abundance filtering and
multiple-testing correction. No p-values or FDR are fabricated here.
MultiQC includes a summary of candidates and replicate thresholds.

## Genomic-window PCA

With `--run_windows`, the pipeline also generates exploratory per-sample PCA
for every requested `--window_sizes` value. Use `--run_pca false` to disable it.
For example:

```bash
nextflow run main.nf [your existing options] \
  --run_windows --window_sizes 10000,50000,100000 -resume
```

Use a new `--outdir` to preserve previous reports. Results are published under
`Analysis/PCA`: `windows_PCA.png`, `windows_PCA.pdf`, sample coordinates,
selected regions/loadings, parameters and software versions. MultiQC includes
a table of explained variance and correlations with original library depth.
Matplotlib is declared in `environment.yml`; the existing activated environment
must provide it when Nextflow-managed Conda is not used.

PCA uses `log2(1 + CPM)`, with centered features and no variance scaling. A
condition-independent filter requires at least five pooled molecules and
nonzero counts in at least two samples; up to 2,000 windows are selected by
variance. The second row repeats PCA after sampling molecules without
replacement to the smallest library size (seed 1729), using the same features.
This is a depth sensitivity check, not a stability analysis or differential
test. Sample group, replicate and donor metadata are preserved; pairing is
not inferred. Window sizes with fewer than three samples, empty libraries,
insufficient features or no variance are reported as skipped.

The analysis task remains at 1 CPU, with BLAS threads limited to one. Sparse
counts and feature selection can retain depth-related patterns after CPM or
equal-depth sampling. PCA separation does not establish a treatment effect
or absolute DSB burden.

## Controlled aligner pilot

Before changing the production aligner, compare BWA, Bowtie2 end-to-end and
Bowtie2 local on identical prepared single-end reads:

```bash
python scripts/pilot_aligner_comparison.py \
  --prepared-dir /path/to/completed_run/prepared \
  --samples HD1_ACUTE HD3_CHRONIC \
  --reference-config /path/to/references_hg38.json \
  --outdir /path/to/new_aligner_pilot \
  --max-reads 100000 --seed 1729 --map-threads 6 --sort-threads 1
```

Both BWA and Bowtie2 indices must already be configured for the same reference.
The script samples uniformly over each complete prepared FASTQ and reuses that
subset for all three branches. It reads the full FASTQs once for sampling;
alignment branches run sequentially. Allow 8 CPUs and 16 GB RAM with these
thread settings. Existing output directories are rejected.

`aligner_comparison.tsv` reports mapping, strict endpoint acceptance and molecule
counts. `aligner_comparison.json` includes parameters, mapping provenance and
pairwise agreement of read-level chromosome, 5-prime coordinate and strand.
Each sample also has `bowtie2_recovered_reads.tsv` with MAPQ, CIGAR, edit count
and sequence prefixes for reads accepted only by Bowtie2, plus per-branch
`five_prime.json` diagnostics. Exact T7 motif flags are diagnostic, not an
artifact classification. Reference headers are compared; matching headers do
not establish identical reference nucleotide sequences.

MAPQ values are aligner-dependent. More accepted reads alone do not demonstrate
more reliable DSB signal, and sampled molecule counts do not estimate full-library
complexity. Review coordinate concordance and recovered-read sequences before
selecting an aligner or changing the strict endpoint filter.

## Validation

Activate the project environment and put Nextflow on PATH, then run:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests
```

Tests cover preparation, design validation, reference/mapping rules,
deduplication, reporting and lightweight Nextflow wiring. Integration tests
require their external tools and skip when unavailable. Local Nextflow checks
have used 24.10.3; HPC compatibility has also been exercised during runs with
26.04.6. A successful run validates execution, not biological specificity.
