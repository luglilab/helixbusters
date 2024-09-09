import os
import subprocess

# Genome URLs for fasta and gff files
genomes = {
    'GRCh38': {
        'fasta': 'ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
        'gff': 'ftp://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz'
    },
    'hg38': {
        'fasta': 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
        'gff': 'ftp://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz'
    },
    'GRCm38': {
        'fasta': 'ftp://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz',
        'gff': 'ftp://ftp.ensembl.org/pub/release-104/gff3/mus_musculus/Mus_musculus.GRCm38.104.gff3.gz'
    },
    'mm10': {
        'fasta': 'ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz',
        'gff': 'ftp://ftp.ensembl.org/pub/release-104/gff3/mus_musculus/Mus_musculus.GRCm38.104.gff3.gz'
    }
}


def download_file(url, output_dir):
    """Download a file using wget."""
    file_name = os.path.join(output_dir, os.path.basename(url))
    if not os.path.exists(file_name):
        subprocess.run(['wget', '-P', output_dir, url])
    else:
        print(f"{file_name} already exists. Skipping download.")
    return file_name


def extract_gz(file_path):
    """Extract .gz file."""
    if file_path.endswith(".gz"):
        subprocess.run(['gunzip', file_path])


def bwa_index(fasta_file):
    """Generate BWA index for the given FASTA file."""
    subprocess.run(['bwa', 'index', fasta_file])


def get_genome_data(genome_build, output_dir="./genome_data"):
    """Download and process genome data by generating the BWA index."""
    if genome_build not in genomes:
        raise ValueError(f"{genome_build} is not a valid genome build.")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Get URLs for the selected genome build
    fasta_url = genomes[genome_build]['fasta']
    gff_url = genomes[genome_build]['gff']

    # Download FASTA and GFF files
    fasta_file = download_file(fasta_url, output_dir)
    gff_file = download_file(gff_url, output_dir)

    # Extract .gz files
    extract_gz(fasta_file)
    extract_gz(gff_file)

    # Generate BWA index
    fasta_file_unzipped = fasta_file[:-3]  # remove .gz extension
    bwa_index(fasta_file_unzipped)

    print(f"BWA index created for {fasta_file_unzipped}")
