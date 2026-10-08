# Install and prepare a first run on an HPC

This guide prepares an isolated Conda/Nextflow environment and validates the
supplied samplesheet before scheduling analysis. Run CPU- and memory-intensive work in
an allocated compute job, following the scheduler rules at your institution;
do not run a full mapping job on a login node.

## 1. Get the project onto the cluster

From a cluster login node, clone the repository into a suitable project or
software directory. If outbound GitHub access is disabled, transfer a clean
source archive from your workstation using the approved HPC transfer route.

```bash
git clone https://github.com/luglilab/helixbusters.git
cd helixbusters
```

If you already have a checkout, update it with `git pull` only when you are
ready to use the latest repository version.

## 2. Create the environment

Load the Conda/Mamba module used by your HPC, if required, then run:

```bash
conda env create --file environment.yml
conda activate helixbusters
python -m pip install --no-deps -e .
module load nextflow/26.04.6
nextflow -version
```

The Conda environment includes Python 3.12, pandas, NumPy, Excel readers,
pysam, BWA, Bowtie2, samtools, cutadapt and UMI-tools. Nextflow is provided by
the HPC module system. If the Conda environment already exists, update it with
`conda env update --file environment.yml --prune` and repeat the editable
install after pulling code changes.

The cluster has Nextflow modules from `nextflow/21.04.3` through
`nextflow/26.04.6`. Use the newest available module for this workflow:

```bash
module load nextflow/26.04.6
nextflow -version
```

If module names are not visible in a fresh shell, run `module avail nextflow`
and use the cluster's module initialization instructions.

Confirm that the commands are available:

```bash
python --version
python -c 'import pandas, openpyxl, pysam, helixbusters; print("Python imports OK")'
bwa 2>&1 | head -n 2
bowtie2 --version | head -n 1
samtools --version | head -n 1
cutadapt --version
umi_tools --version
nextflow -version
```

## 3. Stage the samplesheet and FASTQs

Copy `samplesheet_helixbuster.xlsx` to a cluster-visible project directory.
The paths inside `PathReadForward` (and `PathReadReverse` for paired-end data)
must be paths on the HPC filesystem, not paths from your Mac. Do not edit the
sample/group labels just to encode a time point: `Group` may be `ACUTE` and
`CHRONIC`, while a separate experimental design records which group is baseline
or stimulated where relevant.

The current example workbook contains six **single-end** samples: three per
group, and one eight-base sample barcode per row. It does not by itself tell us
which reference build is correct. Confirm the organism/build and the actual
FASTQ locations with the experiment/reference owner before mapping.

Check that the workbook is readable and every input path exists:

```bash
python - samplesheet_helixbuster.xlsx <<'PY'
import os
import sys
import pandas as pd

sheet = pd.read_excel(sys.argv[1])
required = {"Sample", "Replicate", "Group", "PathReadForward", "SampleBarcodeForward"}
missing = required - set(sheet.columns)
if missing:
    raise SystemExit(f"Missing required columns: {sorted(missing)}")
if sheet["Sample"].isna().any() or sheet["Sample"].duplicated().any():
    raise SystemExit("Sample names must be present and unique")
paired = {"PathReadReverse", "SampleBarcodeReverse"} <= set(sheet.columns)
columns = ["PathReadForward"] + (["PathReadReverse"] if paired else [])
bad = [(row.Sample, column, getattr(row, column))
       for row in sheet.itertuples(index=False)
       for column in columns
       if not isinstance(getattr(row, column), str)
       or not os.path.isfile(getattr(row, column))]
if bad:
    for sample, column, path in bad:
        print(f"Missing input: sample={sample}, column={column}, path={path}")
    raise SystemExit(1)
print(f"Samplesheet OK: {len(sheet)} samples; {'paired-end' if paired else 'single-end'}")
print("Groups:", sheet["Group"].value_counts(dropna=False).to_dict())
PY
```

Run this check from a compute allocation if metadata or the FASTQ filesystem
is not accessible on the login node. Keep raw input files read-only and write
results to a project/scratch location with enough capacity for decompressed
FASTQs and BAMs.

## 4. Configure the local reference catalog

Copy `docs/references.example.json` to a private cluster-side file such as
`references.json`. Replace index prefixes and blacklist paths with resources
already installed on the HPC. Keep each blacklist's `genome` label equal to its
selected build. Use `hg19`, `hg38`, `mm10`, or `mm39`; do not choose based on
the human/mouse species alone. See [genome configuration](genomes.md) for the
canonical-chromosome policy and compatibility checks.

The project does not download or build genome indexes, and it does not fetch
blacklists. A complete reference index and the matching blacklist must exist
before mapping can start.

## 5. Submit a small pilot

Use your HPC's scheduler to request the CPUs, memory, walltime and scratch space
appropriate for the sample size. Activate the environment inside the submitted
job. Start with one sample and a small FASTQ subset to validate barcode/UMI
orientation, output paths and resource use before scaling to all six samples.

Follow the [Nextflow guide](nextflow.md) to convert the workbook to a validated
manifest and run the DSL2 workflow. Start with one sample and a reduced FASTQ;
then scale to the full dataset once UMI/barcode orientation, reference paths,
resource use and QC are correct. The workflow currently supports single-end
samplesheets only and rejects paired-end layouts.

For this workbook the sequencing is single-end, so downstream BLISS end
coordinates are derived from that read. The filtered BAM applies the selected
build's canonical nuclear chromosome and blacklist rules and excludes
mitochondrial alignments. The `quality`/`min_mapq` option is alignment MAPQ;
it is independent of per-base Q30 filtering. Keep the BAM mapping JSON report
with the run for read-retention counts and parameter provenance.

## Troubleshooting

- **`Missing optional dependency 'openpyxl'`:** recreate/update the Conda
  environment from the repository's `environment.yml`.
- **FASTQ not found:** the workbook still contains workstation paths or the
  input filesystem is not mounted in the job; replace them with cluster paths.
- **Incomplete index or missing blacklist:** request the correct build-specific
  resource paths from the reference administrator and update `references.json`.
- **Conda cannot solve/download packages:** use the site's supported Conda
  mirror or environment-module procedure; do not silently drop mapping tools.
- **Build mismatch:** check that the index and blacklist both match the chosen
  assembly. Header chromosome lengths are checked by the pipeline, but cannot
  authenticate the actual reference sequence.
