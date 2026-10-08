# First Nextflow workflow

The initial DSL2 workflow handles the single-end layout in the SP036 example:
samplesheet conversion, UMI/sample-barcode extraction, mapping, canonical
chromosome and blacklist filtering, and per-sample UMI deduplication. It keeps
mapping and deduplication QC reports and publishes outputs below `results/`.
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

On the HPC, Nextflow is supplied by the module system (available versions
range from 21.04.3 to 26.04.6); use `module load nextflow/26.04.6`. The Python
Conda environment does not install or pin Nextflow. If your site initializes
modules differently, follow its shell setup instructions.

## Prepare a manifest

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
Nextflow creates `work/`, `results/`, `timeline.html`, `report.html`, `trace.txt`
and `dag.html` relative to the launch directory unless paths are specified.

## Stages and outputs

1. `EXTRACT_UMI` uses the regex `(?P<umi_1>.{N})<sample barcode>` at the start
   of each read, removes the UMI and barcode from the sequence, and appends the
   UMI to the read name using UMI-tools.
2. `MAP_READS` runs BWA-MEM or Bowtie2 and samtools through the tested Python
   mapping API. Its filtered BAM excludes low/unknown MAPQ, mitochondrial,
   noncanonical and blacklist-overlapping alignments. MAPQ is alignment
   confidence, not a per-base Q30 filter.
3. `DEDUPLICATE` calls the coordinate/strand-aware deduplication implementation
   and writes families, counts, sites, molecules and QC per sample.

Outputs are copied to `outdir/prepared`, `outdir/mapping`, and
`outdir/deduplication`. Keep `work/` until the run completes successfully; it
contains task logs and allows Nextflow resume:

```bash
nextflow run main.nf ... -resume
```

The first pilot must confirm that the barcode is physically positioned directly
after the UMI in R1 and that the barcode sequence is correct. The extraction
pattern describes that layout; a successful Nextflow run cannot establish that
the library design was interpreted correctly. Retain `report.html`,
`timeline.html`, `trace.txt`, and each `*.mapping.json`/`*.dedup.json` with the
analysis.
