```markdown
# Helixbusters

**Helixbusters** is a Python-based pipeline for processing next-generation sequencing (NGS) data with Unique Molecular Identifier (UMI) extraction, adapter trimming, and BWA-based read alignment. It supports both single-end and paired-end sequencing and produces output files for UMI counts and PCR duplicates.

## Features

- **UMI extraction**: Extract UMIs from single-end or paired-end FASTQ files.
- **Adapter trimming**: Trim adapters from reads using `cutadapt` for both single-end and paired-end data.
- **BWA mapping**: Align the trimmed reads to the genome using BWA and filter alignments based on mapping quality.
- **UMI counting**: Generate output files for UMI counts per location and PCR duplicates.

## Dependencies

The following Python packages and external tools are required to run the pipeline:

- **Python packages**:
  - `pandas`
  - `pysam`
  - `subprocess`

- **External tools**:
  - `bwa`
  - `samtools`
  - `cutadapt`
  - `umi_tools`

You can install the required Python packages with:

```bash
pip install pandas pysam
```

Ensure that `bwa`, `samtools`, `cutadapt`, and `umi_tools` are installed and available in your system's PATH.

## Installation

Clone the repository and install the dependencies:

```bash
git clone https://github.com/yourusername/helixbusters.git
cd helixbusters
pip install -r requirements.txt
```

## Usage

1. **Prepare the samplesheet**: Create an Excel or CSV file that contains the information about your samples, including file paths, sample barcodes, and modality (single-end or paired-end sequencing). The required columns are:
   - `Sample`, `Replicate`, `Group`, `PathReadForward`, `SampleBarcodeForward`, `PathReadReverse`, `SampleBarcodeReverse`.

2. **Run the pipeline**:

```python
from helixbusters.core import Helixbusters

# Initialize Helixbusters with your samplesheet, species, mismatch, and genome index path.
helixbusters = Helixbusters(samplesheet="path/to/samplesheet.xlsx", species="human", mismatch=1, genome_index="path/to/genome_index")

# Read the samplesheet
helixbusters.read_column_from_excel()

# Create output folders for each sample
helixbusters.create_sample_output_folders(output_folder="output/directory")

# Process UMI extraction and trimming
helixbusters.process_infofile(umi_length=8, threads=4)

# Run BWA mapping and filter BAM files by quality
helixbusters.run_bwa_mapping(quality=20, threads=10)

# Generate UMI output files for all samples
helixbusters.generate_umi_output_for_samples()
```

## Pipeline Steps

1. **UMI extraction**: The pipeline extracts UMIs from single-end or paired-end FASTQ files, using parallel execution for large datasets.
2. **Adapter trimming**: After UMI extraction, adapters are trimmed from the reads using `cutadapt`, which can handle both single-end and paired-end sequencing.
3. **BWA mapping**: The trimmed reads are aligned to a reference genome using `bwa mem`, followed by sorting and filtering based on mapping quality using `samtools`.
4. **UMI counting**: For each sample, UMI counts are generated per chromosome and location, and PCR duplicates are identified.

## Example Workflow

Here’s an example workflow for processing human samples using paired-end sequencing data:

```python
helixbusters = Helixbusters(
    samplesheet="data/samplesheet.xlsx", 
    species="human", 
    mismatch=1, 
    genome_index="/path/to/bwa/index"
)
helixbusters.read_column_from_excel()
helixbusters.create_sample_output_folders(output_folder="results")
helixbusters.process_infofile(umi_length=8, threads=8)
helixbusters.run_bwa_mapping(quality=30, threads=8)
helixbusters.generate_umi_output_for_samples()
```

## Output

The pipeline generates several outputs for each sample:
- **Trimmed reads**: FASTQ files after UMI extraction and adapter trimming.
- **BAM files**: Aligned and filtered BAM files.
- **UMI output files**:
  - `Chromosome-Location-Strand-UMI-PCR.txt`: Contains UMI counts per chromosome, location, and strand.
  - `Chromosome-Location-UMI-Count.bed`: Contains unique UMI counts per location.

## License

This project is licensed under the MIT License.

