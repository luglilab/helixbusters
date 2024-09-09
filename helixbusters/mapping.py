import os
import gzip
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def download_and_index_genome_with_gff_biopython(genome_name, base_dir):
    """
    Downloads the specified genome reference, GFF file, and BWA index files (if available), and processes the genome with Biopython.
    The files will be organized in subdirectories named after the genome (e.g., 'hg38', 'hg19', 'mm10', 'GRCm38').

    Args:
        genome_name (str): The genome to download ('hg38', 'hg19', 'mm10', 'GRCm38').
        base_dir (str): Base directory where subdirectories for each genome will be created.

    Supported genomes:
        - 'hg38' for human genome version hg38
        - 'hg19' for human genome version hg19
        - 'mm10' for mouse genome version mm10
        - 'GRCm38' for mouse genome version GRCm38
    """

    # UCSC and Ensembl download URLs for genome and GFF files
    genome_urls = {
        'hg38': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
        'hg19': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz',
        'mm10': 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz',
        'GRCm38': 'ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'
    }

    gff_urls = {
        'hg38': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.gtf.gz',
        'hg19': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.gtf.gz',
        'mm10': 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.gtf.gz',
        'GRCm38': 'ftp://ftp.ensembl.org/pub/release-100/gff3/mus_musculus/Mus_musculus.GRCm38.100.gff3.gz'
    }

    # Pre-built BWA index URLs for some genomes
    bwa_index_urls = {
        'hg38': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.bwa.tar.gz',
        'hg19': 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.bwa.tar.gz',
        'mm10': 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.bwa.tar.gz',
        'GRCm38': None  # BWA index may not be pre-built for some genomes
    }

    # Ensure the genome name is valid
    if genome_name not in genome_urls:
        raise ValueError(f"Unsupported genome: {genome_name}. Supported genomes are 'hg38', 'hg19', 'mm10', 'GRCm38'.")

    # Set up subdirectory for this genome
    genome_dir = os.path.join(base_dir, genome_name)
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)

    # Set up paths for downloading genome, GFF, and BWA index files
    genome_url = genome_urls[genome_name]
    gff_url = gff_urls[genome_name]
    bwa_index_url = bwa_index_urls.get(genome_name)

    genome_fa_gz = os.path.join(genome_dir, f"{genome_name}.fa.gz")
    genome_fa = os.path.join(genome_dir, f"{genome_name}.fa")
    gff_gz = os.path.join(genome_dir, f"{genome_name}.gff.gz")
    gff_file = os.path.join(genome_dir, f"{genome_name}.gff")
    bwa_index_dir = os.path.join(genome_dir, "bwa_index")

    # Step 1: Download the genome if it hasn't been downloaded already
    if not os.path.exists(genome_fa_gz):
        print(f"Downloading {genome_name} genome...")
        download_file(genome_url, genome_fa_gz)

    # Step 2: Decompress the genome if it hasn't been decompressed already
    if not os.path.exists(genome_fa):
        print(f"Decompressing {genome_name} genome...")
        decompress_gzip(genome_fa_gz, genome_fa)

    # Step 3: Download the GFF file if it hasn't been downloaded already
    if not os.path.exists(gff_gz):
        print(f"Downloading {genome_name} GFF file...")
        download_file(gff_url, gff_gz)

    # Step 4: Decompress the GFF file if it hasn't been decompressed already
    if not os.path.exists(gff_file):
        print(f"Decompressing {genome_name} GFF file...")
        decompress_gzip(gff_gz, gff_file)

    # Step 5: Download the BWA index files (if available)
    if bwa_index_url:
        if not os.path.exists(bwa_index_dir):
            print(f"Downloading pre-built BWA index for {genome_name}...")
            bwa_index_tar_gz = os.path.join(genome_dir, f"{genome_name}.bwa.tar.gz")
            download_file(bwa_index_url, bwa_index_tar_gz)
            unpack_archive(bwa_index_tar_gz, bwa_index_dir)
            print(f"Pre-built BWA index for {genome_name} downloaded and extracted.")
    else:
        print(f"No pre-built BWA index available for {genome_name}.")

    # Step 6: Process the FASTA file with Biopython (optional step)
    print(f"Processing {genome_name} genome with Biopython...")
    with open(genome_fa, "r") as handle:
        genome_records = list(SeqIO.parse(handle, "fasta"))
        print(f"Number of sequences in {genome_name}: {len(genome_records)}")

    print(f"{genome_name} genome, GFF, and BWA index are ready in {genome_dir}.")


def download_file(url, output_path):
    """
    Download a file from a given URL and save it to a specified path.

    Args:
        url (str): The URL of the file to download.
        output_path (str): The path where the downloaded file will be saved.
    """
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {url} to {output_path}")
    else:
        raise Exception(f"Failed to download {url}")


def decompress_gzip(input_gz_path, output_path):
    """
    Decompress a .gz file and save the output to a specified file.

    Args:
        input_gz_path (str): Path to the .gz file.
        output_path (str): Path where the decompressed file will be saved.
    """
    with gzip.open(input_gz_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            f_out.write(f_in.read())
    print(f"Decompressed {input_gz_path} to {output_path}")


def ensure_bwa_index(genome_fa):
    """
    Ensure that the reference genome is indexed with BWA.
    """
    try:
        index_files = [f"{genome_fa}.{ext}" for ext in ['bwt', 'pac', 'ann', 'amb', 'sa']]
        if all(os.path.exists(f) for f in index_files):
            print("BWA index files found. Skipping indexing.")
        else:
            print(f"Indexing genome {genome_fa} with BWA...")
            bwa_index_cmd = ['bwa', 'index', genome_fa]
            subprocess.run(bwa_index_cmd, check = True)
            print(f"Indexing completed for {genome_fa}.")
    except subprocess.CalledProcessError as e:
        print(f"Error during BWA indexing: {e}")


def map_reads_with_bwa(trimmed_reads, genome_fa, output_dir, paired_end=False):
    """
    Maps trimmed reads to the reference genome using BWA.
    """
    try:
        ensure_bwa_index(genome_fa)

        os.makedirs(output_dir, exist_ok = True)

        sample_name = os.path.basename(trimmed_reads[0] if paired_end else trimmed_reads).split('.')[0]
        output_sam = os.path.join(output_dir, f"{sample_name}_aligned.sam")

        if paired_end:
            read1, read2 = trimmed_reads
            bwa_cmd = ['bwa', 'mem', genome_fa, read1, read2]
        else:
            bwa_cmd = ['bwa', 'mem', genome_fa, trimmed_reads]

        print(f"Mapping reads with BWA for sample: {sample_name}...")
        with open(output_sam, 'w') as sam_file:
            subprocess.run(bwa_cmd, stdout = sam_file, check = True)

        print(f"Reads aligned to {genome_fa}. Output SAM file: {output_sam}")
        return output_sam
    except subprocess.CalledProcessError as e:
        print(f"Error running BWA: {e}")


def read_fastq(filepath):
    """
    Reads the trimmed FASTQ file using Biopython.
    """
    try:
        with open(filepath, "r") as handle:
            fastq_records = list(FastqGeneralIterator(handle))
            print(f"Loaded {len(fastq_records)} records from {filepath}")
            return fastq_records
    except Exception as e:
        print(f"An error occurred while reading the FASTQ file: {e}")
