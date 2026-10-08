![Helixbusters](https://github.com/luglilab/helixbusters/blob/master/logo.png)

# Helixbusters

Helixbusters analyzes BLISS-like sequencing data. The Nextflow DSL2 workflow
coordinates UMI/barcode extraction, mapping, alignment filtering and UMI-based
deduplication; the Python package implements mapping and deduplication.

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
Bowtie2, samtools, cutadapt and UMI-tools. It has no machine-specific Conda
prefix. Nextflow is loaded separately through the HPC module system; see
[HPC installation](docs/hpc.md).

## Start a Nextflow run

The first workflow supports **single-end** samplesheets, including the layout
in the SP036 example. Convert the workbook on a machine where the FASTQs are
visible:

```bash
python scripts/samplesheet_to_manifest.py \
  /path/on/cluster/samplesheet_helixbuster.xlsx \
  /path/on/cluster/samples.tsv \
  --check-fastq
```

Copy `docs/references.example.json` to a cluster-side `references.json` and
replace the placeholders with the existing reference index and matching
blacklist paths. Then submit a pilot through the local scheduler profile:

```bash
nextflow run main.nf \
  --manifest samples.tsv \
  --genome hg38 \
  --reference_config references.json \
  --aligner bwa \
  --map_threads 8 \
  --sort_threads 1 \
  --outdir results/pilot_hg38 \
  -profile slurm
```

Choose the correct build (`hg19`, `hg38`, `mm10`, or `mm39`) from the experiment
reference. The example above uses `hg38` only as syntax. Paired-end input is
rejected in this first workflow. See the full [Nextflow guide](docs/nextflow.md)
before running all samples.

## Analysis behavior

The workflow extracts the UMI and sample barcode from the start of R1 using
`(?P<umi_1>.{N})<barcode>`, retaining the UMI in the read name. Confirm this
matches the library design before interpreting results. Mapping applies the
selected build's canonical nuclear chromosome set, removes mitochondrial and
blacklist-overlapping reads, and records QC. The mapping `quality` setting is
alignment MAPQ, not per-base Q30. See [genome configuration](docs/genomes.md),
[mapping](docs/mapping.md), and [deduplication](docs/deduplication.md).

The pipeline does not download or build genome indexes or blacklists. All
reference catalog paths should identify existing resources for the same build.
