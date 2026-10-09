# Single-end Nextflow workflow

The initial DSL2 workflow handles the single-end layout in the SP036 example:
samplesheet conversion, UMI/sample-barcode extraction, mapping, canonical
chromosome and blacklist filtering, per-sample UMI deduplication, mapping QC,
bigWig generation and condition-level aggregation. It keeps mapping and
deduplication QC reports and publishes outputs below `results/`.
Paired-end support and scheduler-specific resource tuning are not part of this
first workflow yet.

## Install

Create the environment from the repository root and install the package:

```bash
conda env create -f environment.yml
conda activate helixbusters
python -m pip install --no-deps -e .
module load nextflow/26.04.6
nextflow -version
```

The reporting workflow requires Nextflow >=24.10 (explicit input arity keeps
single-replicate conditions as lists). On the HPC, Nextflow is supplied by the module system (available versions
range from 21.04.3 to 26.04.6); use `module load nextflow/26.04.6`. The Python
Conda environment does not install or pin Nextflow. If your site initializes
modules differently, follow its shell setup instructions.

## Prepare a manifest

Declare the experimental design with `--design paired`, `--design unpaired`,
or `--design unspecified` (backward-compatible default). Replicate numbers and
sample names are never used to infer donor identity.

For paired experiments, add a `Donor` column to the Excel/CSV samplesheet. The
converter preserves it as `donor` in the TSV manifest. Use the same stable donor
identifier across conditions, while `Replicate` identifies a biological library
within each condition. Every donor must occur exactly once in each declared
condition. Incomplete pairing is rejected; it requires an explicitly designed
downstream model rather than automatic dropping of samples.

```bash
python scripts/samplesheet_to_manifest.py samples.xlsx samples.paired.tsv --design paired
nextflow run main.nf --manifest samples.paired.tsv --design paired ...
```

For independent samples, use `--design unpaired`. Donor metadata is optional,
but declared donor identities cannot repeat across conditions. For experiments
whose pairing is not yet established, keep `unspecified`; this is valid for
mapping/QC but must be resolved before statistical inference. Different treatment
labels do not make repeated measurements from the same donor independent.

The worker preflight validates design before extraction/mapping, writes
`MultiQC/design.summary.json` and `design.metadata.tsv`, and includes experimental
units in a dedicated MultiQC table. Donor metadata is also propagated in workflow
sample metadata. Formulas `~ donor + condition` and `~ condition` are recorded
as intended designs, not fitted models. This change does not implement differential
analysis, MACS3, or consensus peaks; condition merging remains by `group` for both
designs, and biological replicate support must remain distinct from donor pairing.

The workflow uses a tab-separated manifest rather than reading Excel inside
Nextflow. Convert and validate the workbook on the cluster, where FASTQ files
are visible:

```bash
python scripts/samplesheet_to_manifest.py \
  /path/on/cluster/samplesheet_helixbuster.xlsx \
  /path/on/cluster/samples.tsv \
  --check-fastq
```

This checks required column names, unique sample IDs, A/C/G/T barcode content,
single-end layout and (with `--check-fastq`) that every FASTQ path exists. The
FASTQ check also verifies that `.gz` files are gzip-compressed and that each
input begins with a structurally valid FASTQ record. It reads only the first
record, not the full file. This first workflow intentionally rejects paired-end columns; it must not silently
interpret a paired library as single-end. The FASTQ path in the manifest is
absolute (relative workbook paths are resolved against the workbook folder).

## Configure the selected reference

Copy `docs/references.example.json` on the cluster and replace its paths with
the installed index prefixes and matching blacklist. Paths in the reference
catalog must be absolute because Nextflow stages the catalog into task work
directories. Select the actual build from the experiment metadata/reference
provider; the samplesheet does not identify hg19 versus hg38 (or a mouse build).

For BWA, you can prepare the index and generate a catalog automatically from
the Illumina iGenomes UCSC archive. The helper also downloads the corresponding
Boyle-Lab ENCODE blacklist for hg19, hg38 or mm10. The current Boyle-Lab list
does not include mm39, and the iGenomes table does not list mm39; neither is
silently substituted with an older build. These files are ENCODE blacklists;
cite the Boyle-Lab paper when reporting analyses that use them.

```bash
python scripts/prepare_igenome.py \
  --genome hg38 \
  --cache-dir /project/references/igenomes \
  --output-config /project/references.json \
  --rate-limit 10M
```

It downloads the official iGenomes archive, extracts only its classic BWA
index files, and writes a `references.json` pointing to the cached index and
downloaded blacklist. It records source URLs and SHA-256 values in provenance
JSON files beside the cached resources. The iGenome archive can be large; use a
cluster filesystem with adequate temporary and persistent space and outbound
HTTPS access. This prepares the reference once; the regular Nextflow mapping command then uses
`--reference_config /project/references.json` as before. For Bowtie2 or mm39, use a manually
prepared reference catalog; for mm39, supply a blacklist that is explicitly
matched to mm39.

The downloaded Boyle-Lab files are `hg19-blacklist.v2.bed.gz`,
`hg38-blacklist.v2.bed.gz`, and `mm10-blacklist.v2.bed.gz`. The helper validates
gzip/BED structure; Helixbusters validates chromosome names and coordinates
against the selected build when mapping begins. Source: [Boyle-Lab blacklist
lists](https://github.com/Boyle-Lab/Blacklist/tree/master/lists) and Amemiya,
Kundaje & Boyle, [The ENCODE Blacklist](https://doi.org/10.1038/s41598-019-45839-z).

## Run locally or on Slurm

Start with one sample and a small representative FASTQ subset. For a local
smoke test:

```bash
nextflow run main.nf \
  --manifest samples.tsv \
  --genome hg38 \
  --reference_config references.json \
  --aligner bwa \
  --map_threads 4 \
  --sort_threads 1 \
  --outdir results/pilot_hg38 \
  -profile local
```

Choose `hg19`, `hg38`, `mm10`, or `mm39` only after confirming the assembly.
For an HPC using Slurm, submit from a login node using the site's supported
Nextflow/Java setup and a project/scratch work directory:

```bash
nextflow run main.nf \
  --manifest samples.tsv \
  --genome hg38 \
  --reference_config references.json \
  --aligner bwa \
  --map_threads 8 \
  --sort_threads 1 \
  --queue YOUR_PARTITION \
  --work_dir /path/to/scratch/helixbusters_work \
  --outdir /path/to/project/helixbusters_results \
  -profile slurm
```

Replace `YOUR_PARTITION` with a permitted queue name. Add any account/resource
directives required by the site's Nextflow configuration. On this HPC,
`--queue <partition>` is required because Slurm has no default partition. Find
the permitted partition names with:

```bash
sinfo -o '%P %a %l %D'
```

Use the partition name from the first column, without a trailing `*` (which
marks the default partition when one exists). Then include, for example,
`--queue compute` in the Nextflow command. The defaults in
`nextflow.config` are generic starting values and must be adjusted to local
scheduler policy. Do not launch the entire workflow directly on the login node.
The mapping task reserves `map_threads + sort_threads + 1` CPUs because BWA
and samtools sort run concurrently and `samtools sort -@ N` adds N workers to
its main thread. If Slurm rejects the node request, lower the thread values to
fit a node in the selected partition; `--map_threads 8 --sort_threads 1`
reserves 10 CPUs.
Nextflow creates `work/`, `results/`, `timeline.html`, `report.html`, `trace.txt`
and `dag.html` relative to the launch directory unless paths are specified.

## Stages and outputs

1. `EXTRACT_UMI` streams the original FASTQ through `prepare_bliss_reads.py`.
   It requires an exact barcode immediately after the first `--umi_length`
   bases, removes both UMI and barcode, and appends the UMI to the read name.
   `--barcode_orientation` is `forward` by default; the AA023 EXP1 libraries
   require `reverse_complement`. UMIs containing N are excluded. Inserts must
   have at least `--minimum_insert_length` bases (default 20).
   An optional `--technical_prefix` excludes entire reads whose post-barcode
   sequence starts with that motif. `--technical_prefix_max_errors` defaults
   to 0 (exact matching); 1 or 2 allow anchored Levenshtein edits, including
   substitutions, insertions and deletions. It does not trim or rescue these
   reads, and it does not filter internal occurrences. Parameters and exclusion
   counts are recorded in preparation JSON and MultiQC.
   No fixed technical-sequence length is removed beyond UMI and barcode.
   Each sample writes preparation JSON with counts and parameters and a
   MultiQC preparation table, starting from original input records. Exclusion
   counts are exclusive; reads with multiple issues receive the first reason.
2. `MAP_READS` runs BWA-MEM or Bowtie2 and samtools through the tested Python
   mapping API. Its filtered BAM excludes low/unknown MAPQ, mitochondrial,
   noncanonical and blacklist-overlapping alignments. MAPQ is alignment
   confidence, not a per-base Q30 filter.
3. `DEDUPLICATE` calls the coordinate/strand-aware deduplication implementation
   and writes families, counts, sites, molecules and QC per sample.

Outputs are copied to `outdir/prepared`, `outdir/SingleReplicate/<sample>`,
`outdir/MergedReplicate/<group>` and `outdir/MultiQC`. This replaces the
previous flat `mapping/` and `deduplication/` output layout for new runs. Keep `work/` until the run completes successfully; it
contains task logs and allows Nextflow resume:

```bash
nextflow run main.nf ... -resume
```

MultiQC presents Samtools plots in two separate sections, in order:
`SingleReplicate — Samtools` for `sample__*.txt` (individual all/filtered BAMs),
then `MergedReplicate — Samtools` for `condition__*.txt` (filtered condition pools).
The generated `helixbusters_multiqc_config.json` selects inputs for each section
and excludes noncanonical contigs from chromosome plots. The underlying complete
QC files are retained. Per-sample and per-condition Helixbusters tables remain
separate; pooled BAMs are descriptive outputs, not additional biological replicates.

For AA023 EXP1, the oligo order form specifies eight degenerate bases and the
FASTQ diagnostic finds the reverse-complement barcode at offset 8. Use:

```bash
--barcode_orientation reverse_complement --umi_length 8 \
--minimum_insert_length 20 --technical_prefix CCCTATAGTGAGTCGTAT
```

The HD1_ACUTE pilot on identical 100,000 prepared reads retained 71,741 reads
with the two-edit filter, preserving 11,331 of the baseline 11,335 molecules.
For the full-library comparison, add `--technical_prefix_max_errors 2` and use
a new output directory (`AcuteChronic_technical_filter2`). This is exclusion,
not read rescue. Keep strict 5-prime acceptance, and compare absolute molecule
yield as well as percentages for every biological sample. The pilot does not
establish performance in other samples or resolve all residual clipping.

The motif is observed after the barcode and matches a segment of the STD BLISS
BOTTOM oligo. The prefix filter is conservative and assay-specific: its counts
do not prove that these molecules are adapter dimers. The complete AA023 "new
adapters" structure remains unconfirmed. Inspect preparation retention and a
pilot mapping before interpreting full-library DSB signals. Do not enable this
motif filter for unrelated libraries without evidence. Keep previous results and
use a new output directory for the corrected analysis.

The first pilot must confirm that the barcode is physically positioned directly
after the UMI in R1 and that the barcode sequence is correct. The extraction
pattern describes that layout; a successful Nextflow run cannot establish that
the library design was interpreted correctly. Retain `report.html`,
`timeline.html`, `trace.txt`, and each `*.mapping.json`/`*.dedup.json` with the
analysis.

## Reference downloads and shared iGenomes

The preparer requires `curl`, uses the HTTPS S3 object endpoint, and retains
resumable archives under `<cache-dir>/downloads/*.partial`. Rerun the same
command after an interruption. `--attempts` defaults to 8; diagnostics include
the curl error. HTTP errors, certificate failures and unsupported byte ranges
stop immediately. `--base-url` selects another archive mirror with the same
organism/source/build layout; use a separate cache when changing mirrors.
SHA-256 records provenance; it is not a comparison with a provider checksum.
The archive remains cached after extraction, so allow space for both archive
and index. `100K` limits transfers to roughly 100 KiB/s and can make a large
reference download take days. Output catalogs must have a new filename.

Like nf-core's `igenomes_base`, an existing shared installation can be used:

```bash
python scripts/prepare_igenome.py \
  --genome hg38 \
  --igenomes-base /shared/igenomes \
  --cache-dir ./references/igenomes \
  --blacklist /shared/blacklists/hg38-blacklist.v2.bed.gz \
  --blacklist-genome hg38 \
  --output-config ./references_hg38.json
```

This mode reads `Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa`
(or a version subdirectory) and leaves the shared installation untouched.
Without `--blacklist`, it still downloads the build-matched blacklist.
Keep UCSC hg38 separate from NCBI/Ensembl GRCh38 installations: provider
contents and contig naming can differ even for the same assembly generation.

Nextflow configuration is separated into `conf/base.config`,
`conf/local.config`, and `conf/slurm.config`. Profiles retain the existing
resource defaults. Site-specific overrides can be supplied with
`-c /path/to/site.config`; reference selection still uses the validated JSON
catalog through `--reference_config`.

## Mapping reports and bigWig tracks

The manifest must contain `sample`, `group`, `replicate`, `barcode` and `fastq`.
`group` is the biological condition; `replicate` identifies a biological
replicate within that condition. Each `(group, replicate)` must have exactly
one library. Replicate identifiers may be reused across different conditions.
Sample, group and replicate labels must start with a letter or digit and
contain only letters, digits, `_`, `.` or `-`; labels are never silently renamed.
Technical-library pooling is not supported by this workflow.

The new worker dependencies are `multiqc`, `deeptools` (`bamCoverage`) and
`pybigwig`, listed in `environment.yml`. Update the existing environment after
any running analysis has finished, rather than changing its packages mid-run:

```bash
conda env update -n helixbusters -f environment.yml
conda activate helixbusters
```

`CHECK_ENVIRONMENT` verifies worker executables and the selected index and
blacklist before UMI extraction starts, and records dependency versions in
`MultiQC/environment.json`. Reference paths must be absolute and visible on
all compute nodes. `SAMPLE_QC` requests 2 CPUs, 8 GB RAM and 4 hours;
`CONDITION_QC` requests 2 CPUs, 16 GB and 12 hours. MultiQC requests 1 CPU,
8 GB and 4 hours. Override these generic values in a site config if necessary.
No packages are installed automatically by the workflow.

```text
<outdir>/
  prepared/
  SingleReplicate/<sample>/
    mapping/         # Complete and filtered BAMs, indices, mapping JSON/logs
    deduplication/   # UMI families, molecule/site counts and deduplication JSON
    qc/              # Before/after samtools QC, summary TSV/JSON, provenance
    bigwig/          # Filtered mapping coverage and deduplicated 5-prime ends
  MergedReplicate/<group>/
    mapping/         # Merged filtered BAM/BAI, retaining PCR duplicates
    qc/              # Pooled summary, per-replicate table, merged samtools QC
    bigwig/          # Pooled mapping coverage and raw/CPM/mean end tracks
  MultiQC/
    multiqc_report.html
    multiqc_data/
    helixbusters_samples.tsv
    helixbusters_conditions.tsv
    environment.json
    reporting_versions.txt
```

The report contains two Helixbusters tables (samples and conditions), plus
samtools `flagstat`, `stats` and `idxstats` from complete and filtered sample
BAMs and merged filtered condition BAMs. Samtools entries are named
`sample__<sample>.all`, `sample__<sample>.filtered` and
`condition__<group>.filtered`, avoiding stage/condition name collisions.
See [MultiQC samtools support](https://docs.seqera.io/multiqc/modules/samtools)
and [custom content](https://docs.seqera.io/multiqc/custom_content).

Mapping rate is `mapped_primary_records / primary_records`; retained yield
is `retained_records / primary_records`. Both refer to aligner input **after
UMI/barcode extraction**, not the original FASTQ read count. Filtering removes
unmapped, secondary, supplementary, QC-failing, unknown/low-MAPQ,
mitochondrial, noncanonical and blacklist-overlapping records in that order.
Exclusion categories are mutually exclusive and order-dependent: a low-MAPQ
mitochondrial read is counted as low MAPQ, not mitochondrial. They are not
independent fractions of all reads overlapping the blacklist or mitochondria.
The report includes retention and UMI duplication; the latter is
`duplicate_reads / accepted_reads` in the deduplication stage. Ambiguous
5-prime reads excluded from deduplication are also reported. Samtools
duplicate-flag counts are not UMI duplication estimates; use the Helixbusters
deduplication metrics for molecular duplication.

BigWig files have distinct signal definitions:

| Track suffix | Signal and normalization |
| --- | --- |
| `.coverage.CPM.bw` | Filtered BAM alignment coverage, before UMI deduplication; deepTools CPM with exact scaling and no read extension |
| `.ends.raw.bw` | Number of independently UMI-deduplicated molecules at each 1-bp 5-prime position, strands summed |
| `.ends.CPM.bw` | End counts divided by retained deduplicated molecule total, multiplied by 1e6 |
| `.ends.mean.CPM.bw` | Conditions only: equal-weight mean of individual biological replicate end CPMs |

`--coverage_bin_size` controls mapping coverage bins (default 50 bp). End
tracks retain 1-bp resolution and are written in sparse batches; no dense
whole-genome arrays are allocated. Chromosome names, order and lengths come
from BAM headers and are not rewritten. Raw end counts are retained in
`.ends.counts.bed`. See [bamCoverage](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html)
for its coverage and CPM definitions.

Condition end tracks are produced **after independent library deduplication**.
Raw counts are summed; pooled CPM uses the sum of molecule totals. The mean
CPM track gives each biological replicate equal weight, including zero signal
at positions observed only in other replicates. Libraries with zero molecules
have an empty raw/CPM track and undefined normalization; if any replicate is
empty, the condition mean CPM track is omitted and this is recorded in
provenance. A condition with one replicate has identical pooled and mean CPM.
Condition QC percentages are ratios of pooled counts; replicate means and
sample standard deviations are additionally reported, without significance
tests. Undefined percentages/SDs are missing, not zero.

The condition BAM is a merge of **filtered, non-deduplicated** sample BAMs for
mapping inspection. It must not be deduplicated across biological replicates.
Keep individual samples as experimental units for downstream statistical
analysis. CPM tracks compare relative distributions and do not establish
absolute changes in break burden without an appropriate experimental
normalization strategy.

Allow disk space in both `work/` and the published results for an additional
merged BAM per condition, approximately the sum of filtered replicate BAM
sizes, plus tracks and QC. Use a new output directory when first running this
version; original results are not migrated or deleted. Then resume using the
same work directory as usual. Reporting scripts refuse existing output files.

Focused validation (in the configured environment):

```bash
python -m unittest tests.test_reporting tests.test_reporting_integration tests.test_nextflow_reporting
```

The Nextflow test uses three synthetic samples in two conditions and stub
process outputs to validate metadata, staging, publishing and the numeric
11-CPU request for `--map_threads 5 --sort_threads 5`. It does not map reads.
The reporting tests round-trip real bigWigs and parse the custom sections with
MultiQC when those dependencies are available; missing tools produce explicit
skips. Existing mapping/filtering/deduplication tests remain applicable.

The reporting integration test starts from actual toy-genome BAM files and
executes library QC, independently deduplicated end tracks, a condition BAM
merge, deepTools coverage and a combined MultiQC report. It does not use human
experimental data or claim to validate read alignment against hg38.

Optional downstream analysis: add `--run_windows --window_sizes 1000,5000,10000`
for sparse integer molecule matrices, and/or `--run_peak_calling
--min_reps_consensus 2` for per-sample MACS3 hotspots and condition consensus.
MACS3 must be available on workers only when peak calling is enabled.
See README downstream analysis for smoothing, effective genome size, output
locations and statistical limits. These steps do not fit a differential model.
