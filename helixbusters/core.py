from helixbusters.utils import (
    read_excel_column,
    extract_umi_parallel,
    run_cutadapt_single_end,
    run_cutadapt_paired_end,
    process_bam_and_generate_umi_outputs,
    plot_alignment_quality
)
import os
import subprocess
import pysam
import pandas as pd


class Helixbusters:
    def __init__(self, samplesheet, species, mismatch, genome_index):
        self.samplesheet = samplesheet
        self.species = species.lower()
        self.mismatch = mismatch
        self.genome_index = genome_index  # Path to the genome index
        self.modality = None
        self.infofile = None

        # Validate species attribute
        if self.species not in ['mouse', 'human']:
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

    def process_infofile(self, umi_length=8, threads=1):
        """
        Iterates over the rows of self.infofile, performing UMI extraction and trimming for each sample.

        Args:
        - umi_length (int): Length of the UMI sequence. Defaults to 8.
        - threads (int): Number of threads to use for UMI extraction and trimming. Defaults to 1.
        """
        for index, row in self.infofile.iterrows():
            fastq1_path = row['PathReadForward']
            fastq2_path = row.get('PathReadReverse', None)  # Use PathReadReverse if available (for paired-end)
            barcode_forward = row['SampleBarcodeForward']  # SampleBarcodeForward used for trimming
            barcode_reverse = row.get('SampleBarcodeReverse', None)  # For paired-end data
            output_path = row['OutputPath']

            # Ensure output directory exists
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            print(f"Processing sample {row['Sample']} (row {index + 1}/{len(self.infofile)})")

            # Step 1: UMI extraction
            extract_umi_parallel(
                fastq1_path=fastq1_path,
                fastq2_path=fastq2_path,
                adapter1=barcode_forward,
                adapter2=barcode_reverse,
                output_path=output_path,
                umi_length=umi_length,
                threads=threads,
                sample_barcode=barcode_forward  # Trimming will use SampleBarcodeForward
            )

            # Step 2: Post-processing trimming using cutadapt based on modality
            if self.modality == 'single-end':
                trimmed_r1_path = os.path.join(output_path, "trimmed.fastq.gz")
                run_cutadapt_single_end(
                    read_path=os.path.join(output_path, "final_output.fastq.gz"),
                    output_path=output_path,
                    adapter_seq=barcode_forward,
                    threads=threads
                )
                self.infofile.at[index, 'PathReadForwardTrimmed'] = trimmed_r1_path
            elif self.modality == 'paired-end':
                trimmed_r1_path = os.path.join(output_path, "trimmed_R1.fastq.gz")
                trimmed_r2_path = os.path.join(output_path, "trimmed_R2.fastq.gz")
                run_cutadapt_paired_end(
                    read1_path=os.path.join(output_path, "final_output.fastq.gz"),
                    read2_path=os.path.join(output_path, "final_output_R2.fastq.gz"),
                    output_path=output_path,
                    adapter_seq1=barcode_forward,
                    adapter_seq2=barcode_reverse,
                    threads=threads
                )
                self.infofile.at[index, 'PathReadForwardTrimmed'] = trimmed_r1_path
                self.infofile.at[index, 'PathReadReverseTrimmed'] = trimmed_r2_path

    def run_bwa_mapping(self, quality=20, threads=10):
        """
        Run BWA mapping and Samtools sorting for each sample using the trimmed reads.
        Adds the paths of BAM files ('BamAllPath' and 'BamFilteredPath') to self.infofile.

        Args:
            quality (int): Minimum mapping quality.
            threads (int): Number of threads to use.
        """
        if self.infofile is None or self.modality is None:
            raise ValueError(
                "No sample information or modality found. Please ensure to run read_column_from_excel first.")

        for index, row in self.infofile.iterrows():
            sample = row['Sample']
            output_path = row['OutputPath']
            aux_path = output_path  # Assuming aux files are in the same folder as output

            # Paths for BAM files
            bam_all = os.path.join(output_path, f"{sample}.all.bam")
            bam_filtered = os.path.join(output_path, f"{sample}.q{quality}.bam")

            if self.modality == 'single-end':
                # Single-end mode: using trimmed.fastq.gz
                trimmed_r1_path = os.path.join(output_path, "trimmed.fastq.gz")
                if not os.path.exists(trimmed_r1_path):
                    raise FileNotFoundError(f"Trimmed FASTQ file not found: {trimmed_r1_path}")

                # Single-end command, using the trimmed.fastq.gz file
                bwa_cmd = f"bwa mem -v 1 -t {threads} {self.genome_index} {trimmed_r1_path} | samtools sort --threads {threads} -T {aux_path}/{sample} -o {bam_all}"

            elif self.modality == 'paired-end':
                # Paired-end mode: using trimmed_R1.fastq.gz and trimmed_R2.fastq.gz
                trimmed_r1_path = os.path.join(output_path, "trimmed_R1.fastq.gz")
                trimmed_r2_path = os.path.join(output_path, "trimmed_R2.fastq.gz")

                if not os.path.exists(trimmed_r1_path) or not os.path.exists(trimmed_r2_path):
                    raise FileNotFoundError(f"Trimmed FASTQ files not found: {trimmed_r1_path}, {trimmed_r2_path}")

                # Paired-end command, using the trimmed_R1.fastq.gz and trimmed_R2.fastq.gz files
                bwa_cmd = f"bwa mem -v 1 -t {threads} {self.genome_index} {trimmed_r1_path} {trimmed_r2_path} | samtools sort --threads {threads} -T {aux_path}/{sample} -o {bam_all}"

            else:
                raise ValueError(f"Unsupported modality: {self.modality}")

            # Run BWA and Samtools sorting
            print(f"Running BWA and sorting for sample {sample}...")
            subprocess.run(bwa_cmd, shell=True, check=True)

            # Filter BAM by quality and index BAM files
            print(f"Filtering BAM file for sample {sample} with minimum quality {quality}...")
            view_cmd = f"samtools view --threads {threads} -b -q {quality} {bam_all} > {bam_filtered}"
            subprocess.run(view_cmd, shell=True, check=True)

            # Index both BAM files (all and filtered)
            print(f"Indexing BAM files for sample {sample}...")
            pysam.index(bam_all)
            pysam.index(bam_filtered)

            print(f"BWA mapping and BAM processing complete for sample {sample}.")

            # Store the paths of the BAM files in self.infofile
            self.infofile.at[index, 'BamAllPath'] = bam_all
            self.infofile.at[index, 'BamFilteredPath'] = bam_filtered

    def run_bowtie2_mapping(self, quality=20, threads=10):
        """
        Run Bowtie2 mapping and Samtools sorting for each sample using the trimmed reads.
        Adds the paths of BAM files ('BamAllPath' and 'BamFilteredPath') to self.infofile.

        Args:
            quality (int): Minimum mapping quality.
            threads (int): Number of threads to use.
        """
        if self.infofile is None or self.modality is None:
            raise ValueError(
                "No sample information or modality found. Please ensure to run read_column_from_excel first.")

        for index, row in self.infofile.iterrows():
            sample = row['Sample']
            output_path = row['OutputPath']
            aux_path = output_path  # Assuming aux files are in the same folder as output

            # Paths for BAM files
            bam_all = os.path.join(output_path, f"{sample}.all.bam")
            bam_filtered = os.path.join(output_path, f"{sample}.q{quality}.bam")

            if self.modality == 'single-end':
                # Single-end mode: using trimmed.fastq.gz
                trimmed_r1_path = os.path.join(output_path, "trimmed.fastq.gz")
                if not os.path.exists(trimmed_r1_path):
                    raise FileNotFoundError(f"Trimmed FASTQ file not found: {trimmed_r1_path}")

                # Single-end command, using the trimmed.fastq.gz file
                bowtie2_cmd = f"bowtie2 -x {self.genome_index} -U {trimmed_r1_path} -p {threads} | samtools sort --threads {threads} -T {aux_path}/{sample} -o {bam_all}"

            elif self.modality == 'paired-end':
                # Paired-end mode: using trimmed_R1.fastq.gz and trimmed_R2.fastq.gz
                trimmed_r1_path = os.path.join(output_path, "trimmed_R1.fastq.gz")
                trimmed_r2_path = os.path.join(output_path, "trimmed_R2.fastq.gz")

                if not os.path.exists(trimmed_r1_path) or not os.path.exists(trimmed_r2_path):
                    raise FileNotFoundError(f"Trimmed FASTQ files not found: {trimmed_r1_path}, {trimmed_r2_path}")

                # Paired-end command, using the trimmed_R1.fastq.gz and trimmed_R2.fastq.gz files
                bowtie2_cmd = f"bowtie2 -x {self.genome_index} -1 {trimmed_r1_path} -2 {trimmed_r2_path} -p {threads} | samtools sort --threads {threads} -T {aux_path}/{sample} -o {bam_all}"

            else:
                raise ValueError(f"Unsupported modality: {self.modality}")

            # Run Bowtie2 and Samtools sorting
            print(f"Running Bowtie2 and sorting for sample {sample}...")
            subprocess.run(bowtie2_cmd, shell=True, check=True)

            # Filter BAM by quality and index BAM files
            print(f"Filtering BAM file for sample {sample} with minimum quality {quality}...")
            view_cmd = f"samtools view --threads {threads} -b -q {quality} {bam_all} > {bam_filtered}"
            subprocess.run(view_cmd, shell=True, check=True)

            # Index both BAM files (all and filtered)
            print(f"Indexing BAM files for sample {sample}...")
            pysam.index(bam_all)
            pysam.index(bam_filtered)

            print(f"Bowtie2 mapping and BAM processing complete for sample {sample}.")

            # Store the paths of the BAM files in self.infofile
            self.infofile.at[index, 'BamAllPath'] = bam_all
            self.infofile.at[index, 'BamFilteredPath'] = bam_filtered

    def generate_umi_output_for_samples(self):
        """
        Generates UMI output files for all samples based on their .q20.bam files.
        This method processes each BAM file in the infofile and generates:
        1. Chromosome-Location-Strand-UMI-PCR.txt
        2. Chromosome-Location-UMI-Count.bed
        """
        for index, row in self.infofile.iterrows():
            sample = row['Sample']
            bam_filtered_path = row['BamFilteredPath']  # Path to the .q20.bam file

            if not bam_filtered_path or not os.path.exists(bam_filtered_path):
                print(f"Warning: BAM file not found for sample {sample}. Skipping...")
                continue

            # Define the output file paths
            output_file_umi_pcr = os.path.join(row['OutputPath'], f"{sample}_Chromosome-Location-Strand-UMI-PCR.txt")
            output_file_umi_count = os.path.join(row['OutputPath'], f"{sample}_Chromosome-Location-UMI-Count.bed")

            # Call the utility function to generate UMI outputs
            print(f"Generating UMI output files for sample {sample}...")
            process_bam_and_generate_umi_outputs(bam_filtered_path, output_file_umi_pcr, output_file_umi_count)

            # Optionally, store the paths to the generated files in the infofile
            self.infofile.at[index, 'UMI_PCR_Output'] = output_file_umi_pcr
            self.infofile.at[index, 'UMI_Count_Output'] = output_file_umi_count

        print("UMI output generation completed for all samples.")

    def generate_alignment_quality_plots(self):
        """
        Generate alignment quality distribution plots for all samples in self.infofile.
        The plots are saved in the respective output directories of each sample.
        """
        for index, row in self.infofile.iterrows():
            sample = row['Sample']
            bam_file = row['BamFilteredPath']  # Assuming we are using the filtered BAM file

            if not bam_file or not os.path.exists(bam_file):
                print(f"Warning: BAM file not found for sample {sample}. Skipping...")
                continue

            output_plot_path = os.path.join(row['OutputPath'], f"{sample}_alignment_quality_distribution.png")

            # Generate the alignment quality plot for this sample
            print(f"Generating alignment quality plot for sample {sample}...")
            plot_alignment_quality(bam_file, output_plot_path)

        print("Alignment quality plots generated for all samples.")
