from helixbusters.utils import (
    read_excel_column,
    run_cutadapt_single_end,
    run_cutadapt_paired_end
)
import os
import subprocess
import pysam

class Helixbusters:
    def __init__(self, samplesheet, species, mismatch, genome_index):
        self.samplesheet = samplesheet
        self.species = species.lower()
        self.mismatch = mismatch
        self.genome_index = genome_index  # Now this is a path to the genome index
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

    def trim_reads(self, min_len=20, threads=10):
        """
        Trim reads based on the modality determined from the samplesheet.
        Calls the appropriate function for single-end or paired-end reads.
        Uses the 'OutputPath' column from self.infofile for the output directory.
        Fetches adapter sequences from 'SampleBarcodeForward' and 'SampleBarcodeReverse' (for paired-end).
        Adds the trimmed read paths to new columns 'PathReadForwardTrimmed' (and 'PathReadReverseTrimmed' for paired-end).
        """
        if self.infofile is None or self.modality is None:
            raise ValueError(
                "No sample information or modality found. Please ensure to run read_column_from_excel first.")

        for index, row in self.infofile.iterrows():
            output_path = row['OutputPath']

            if self.modality == 'single-end':
                read_path = row['PathReadForward']
                adapter_seq = row['SampleBarcodeForward']
                trimmed_read_path = os.path.join(output_path, "trimmed.fastq.gz")  # Fixed trimmed filename

                print(f"Trimming single-end reads for sample: {row['Sample']} to {output_path}")
                run_cutadapt_single_end(read_path, output_path, adapter_seq, min_len, threads)

                # Store the path of the trimmed reads in self.infofile
                self.infofile.at[index, 'PathReadForwardTrimmed'] = trimmed_read_path

            elif self.modality == 'paired-end':
                read1_path = row['PathReadForward']
                read2_path = row['PathReadReverse']
                adapter_seq1 = row['SampleBarcodeForward']
                adapter_seq2 = row['SampleBarcodeReverse']

                trimmed_read1_path = os.path.join(output_path, "trimmed_R1.fastq.gz")  # Fixed trimmed filename for R1
                trimmed_read2_path = os.path.join(output_path, "trimmed_R2.fastq.gz")  # Fixed trimmed filename for R2

                print(f"Trimming paired-end reads for sample: {row['Sample']} to {output_path}")
                run_cutadapt_paired_end(read1_path, read2_path, output_path, adapter_seq1, adapter_seq2, min_len,
                                        threads)

                # Store the paths of the trimmed reads in self.infofile
                self.infofile.at[index, 'PathReadForwardTrimmed'] = trimmed_read1_path
                self.infofile.at[index, 'PathReadReverseTrimmed'] = trimmed_read2_path

            else:
                raise ValueError("Unsupported modality: must be 'single-end' or 'paired-end'.")

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

            # Paths for single-end or paired-end reads
            read1 = row['PathReadForwardTrimmed']
            read2 = row.get('PathReadReverseTrimmed', None)

            # Paths for BAM files
            bam_all = os.path.join(output_path, f"{sample}.all.bam")
            bam_filtered = os.path.join(output_path, f"{sample}.q{quality}.bam")

            if self.modality == 'single-end':
                # Single-end command, using PathReadForwardTrimmed
                bwa_cmd = f"bwa mem -v 1 -t {threads} {self.genome_index} {read1} | samtools sort --threads {threads} -T {aux_path}/{sample} -o {bam_all}"
            elif self.modality == 'paired-end':
                # Paired-end command, using PathReadForwardTrimmed and PathReadReverseTrimmed
                bwa_cmd = f"bwa mem -v 1 -t {threads} {self.genome_index} {read1} {read2} | samtools sort --threads {threads} -T {aux_path}/{sample} -o {bam_all}"
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


