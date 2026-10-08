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
