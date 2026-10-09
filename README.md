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
`design.metadata.tsv`, and the MultiQC report. It does not fit a differential model.

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

Genomic-window comparisons, MACS3 hotspot discovery, replicate consensus peaks
(`--min_reps_consensus`), and differential analysis are **planned, not implemented
in `main.nf`**. Experimental-design metadata is available for these future analyses.
Peak significance and differential significance between conditions are distinct.

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
