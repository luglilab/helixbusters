from helixbusters.utils import (
    read_excel_column,
    run_cutadapt_single_end,
    run_cutadapt_paired_end,
    download_and_index_genome_with_gff_biopython,
    map_reads_with_bwa,
    read_fastq
)
import os

class Helixbusters:
    def __init__(self, samplesheet, species, mismatch):
        self.samplesheet = samplesheet
        self.species = species
        self.mismatch = mismatch
        self.modality = None
        self.infofile = None

        # Validate species attribute
        if species.lower() not in ['mouse', 'human']:
            raise ValueError("Species must be 'mouse' or 'human'.")

        # Validate mismatch attribute
        if not 0 <= mismatch <= 3:
            raise ValueError("Mismatch must be between 0 and 3.")

    def read_column_from_excel(self):
        """
        Reads and returns the content of the samplesheet using the read_excel_column function from utils.
        """
        self.infofile, self.modality = read_excel_column(self.samplesheet)

    def create_sample_output_folders(self, output_folder):
        """
        Creates subfolders for each sample in the 'Sample' column within the output_folder.
        Adds a new 'OutputPath' column to self.infofile with the path to each sample's output folder.

        Args:
            output_folder (str): The base folder where subfolders will be created.
        """
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        output_paths = []
        for index, row in self.infofile.iterrows():
            sample_name = row['Sample']
            sample_output_path = os.path.join(output_folder, sample_name)
            if not os.path.exists(sample_output_path):
                os.makedirs(sample_output_path)
            output_paths.append(sample_output_path)

        self.infofile['OutputPath'] = output_paths

    def trim_reads(self, min_len=20, threads=10):
        """
        Trim reads based on the modality determined from the samplesheet.
        Calls the appropriate function for single-end or paired-end reads.
        Uses the 'OutputPath' column from self.infofile for the output directory.
        Fetches adapter sequences from 'SampleBarcodeForward' and 'SampleBarcodeReverse' (for paired-end).
        """
        if self.infofile is None or self.modality is None:
            raise ValueError("No sample information or modality found. Please ensure to run read_column_from_excel first.")

        for _, row in self.infofile.iterrows():
            output_path = row['OutputPath']

            if self.modality == 'single-end':
                read_path = row['PathReadForward']
                adapter_seq = row['SampleBarcodeForward']
                print(f"Trimming single-end reads for sample: {row['Sample']} to {output_path}")
                run_cutadapt_single_end(read_path, output_path, adapter_seq, min_len, threads)

            elif self.modality == 'paired-end':
                read1_path = row['PathReadForward']
                read2_path = row['PathReadReverse']
                adapter_seq1 = row['SampleBarcodeForward']
                adapter_seq2 = row['SampleBarcodeReverse']
                print(f"Trimming paired-end reads for sample: {row['Sample']} to {output_path}")
                run_cutadapt_paired_end(read1_path, read2_path, output_path, adapter_seq1, adapter_seq2, min_len, threads)

            else:
                raise ValueError("Unsupported modality: must be 'single-end' or 'paired-end'.")

    def download_and_index_genome(self, base_dir):
        """
        Downloads the human or mouse genome reference based on the species provided.
        Indexes the genome using BWA.
        """
        if self.species == 'human':
            genome = 'hg38'
        elif self.species == 'mouse':
            genome = 'mm10'
        else:
            raise ValueError(f"Unsupported species: {self.species}")

        print(f"Downloading and indexing genome for {self.species}: {genome}")
        download_and_index_genome_with_gff_biopython(genome, base_dir)

    def map_reads(self, genome_fa, output_dir):
        """
        Maps the trimmed reads of all samples in the samplesheet to the reference genome using BWA.

        Args:
            genome_fa (str): Path to the reference genome FASTA file.
            output_dir (str): Directory where the SAM/BAM output will be saved.
        """
        if self.infofile is None or self.modality is None:
            raise ValueError("No sample information or modality found. Please ensure to run read_column_from_excel first.")

        for _, row in self.infofile.iterrows():
            sample_name = row['Sample']
            output_path = row['OutputPath']

            # Get the trimmed reads for mapping
            if self.modality == 'single-end':
                trimmed_reads = os.path.join(output_path, 'trimmed.fastq.gz')
                print(f"Mapping single-end reads for sample: {sample_name}")
                map_reads_with_bwa(trimmed_reads, genome_fa, output_dir, paired_end=False)

            elif self.modality == 'paired-end':
                trimmed_reads_r1 = os.path.join(output_path, 'trimmed_R1.fastq.gz')
                trimmed_reads_r2 = os.path.join(output_path, 'trimmed_R2.fastq.gz')
                print(f"Mapping paired-end reads for sample: {sample_name}")
                map_reads_with_bwa((trimmed_reads_r1, trimmed_reads_r2), genome_fa, output_dir, paired_end=True)

            else:
                raise ValueError("Unsupported modality: must be 'single-end' or 'paired-end'.")

        print(f"All samples have been mapped to the reference genome.")
