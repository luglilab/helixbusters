from helixbusters.core import Helixbusters
helix = Helixbusters(samplesheet="/home/lugli/spuccio/Projects/SP036_Lise/Dev/samplesheet_tp2.csv", species="human",
mismatch=1,genome_index="/mnt/references/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa",
)
helix.read_column_from_excel()
helix.infofile
helix.modality
helix.create_sample_output_folders("/home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/")
helix.process_infofile(umi_length=8, threads=10)

helix.create_sample_output_folders("/home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/")
helix.trim_reads("/home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/")
helix.lowercase_reads("/home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/")
helix.run_bwa_mapping(quality=20, threads=10)


helix.infofile['PathReadForward'][0]


 #/home/lugli/spuccio/Projects/SP036_Lise/RITM0026254_BLISS_VL_EG/test_sample_pipeline/subSPE_013_LIB_S1_L001_R1_001.fastq.gz

fastp -i /home/lugli/spuccio/Projects/SP036_Lise/RITM0026254_BLISS_VL_EG/test_sample_pipeline/subSPE_013_LIB_S1_L001_R1_001.fastq.gz -o /home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/Sample1/fastp.fastq.gz -a CGTGTGAG -U --umi_loc=read1 --umi_len=8  




from helixbusters.utils import (
    read_excel_column,
    run_cutadapt_single_end,
    run_cutadapt_paired_end,
    run_cutadapt_single_end_lowercase,
    run_cutadapt_paired_end_lowercase,
    run_fastp_single_end,run_fastp_paired_end
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

    def fastp_trim_reads(self, min_len=20, threads=10):
        """
        Trim reads based on the modality determined from the samplesheet using fastp.
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

                print(f"Trimming single-end reads for sample: {row['Sample']} to {output_path} using fastp")
                run_fastp_single_end(read_path, output_path, adapter_seq, min_len, threads)

                # Store the path of the trimmed reads in self.infofile
                self.infofile.at[index, 'PathReadForwardTrimmed'] = trimmed_read_path

            elif self.modality == 'paired-end':
                read1_path = row['PathReadForward']
                read2_path = row['PathReadReverse']
                adapter_seq1 = row['SampleBarcodeForward']
                adapter_seq2 = row['SampleBarcodeReverse']

                trimmed_read1_path = os.path.join(output_path, "trimmed_R1.fastq.gz")  # Fixed trimmed filename for R1
                trimmed_read2_path = os.path.join(output_path, "trimmed_R2.fastq.gz")  # Fixed trimmed filename for R2

                print(f"Trimming paired-end reads for sample: {row['Sample']} to {output_path} using fastp")
                run_fastp_paired_end(read1_path, read2_path, output_path, adapter_seq1, adapter_seq2, min_len, threads)

                # Store the paths of the trimmed reads in self.infofile
                self.infofile.at[index, 'PathReadForwardTrimmed'] = trimmed_read1_path
                self.infofile.at[index, 'PathReadReverseTrimmed'] = trimmed_read2_path

            else:
                raise ValueError("Unsupported modality: must be 'single-end' or 'paired-end'.")

    def lowercase_reads(self, min_len=20, threads=10):
        """
        Apply lowercase action to reads based on the modality determined from the samplesheet.
        Calls the appropriate function for single-end or paired-end reads.
        Uses the 'OutputPath' column from self.infofile for the output directory.
        Fetches adapter sequences from 'SampleBarcodeForward' and 'SampleBarcodeReverse' (for paired-end).
        Adds the lowercase read paths to new columns 'PathReadForwardLowercase' (and 'PathReadReverseLowercase' for paired-end).
        """
        if self.infofile is None or self.modality is None:
            raise ValueError(
                "No sample information or modality found. Please ensure to run read_column_from_excel first.")

        for index, row in self.infofile.iterrows():
            output_path = row['OutputPath']

            if self.modality == 'single-end':
                read_path = row['PathReadForward']
                adapter_seq = row['SampleBarcodeForward']
                lowercase_read_path = os.path.join(output_path, "lowercase.fastq.gz")  # Fixed lowercase filename

                print(f"Applying lowercase action to single-end reads for sample: {row['Sample']} to {output_path}")
                run_cutadapt_single_end_lowercase(read_path, output_path, adapter_seq, min_len, threads)

                # Store the path of the lowercase reads in self.infofile
                self.infofile.at[index, 'PathReadForwardLowercase'] = lowercase_read_path

            elif self.modality == 'paired-end':
                read1_path = row['PathReadForward']
                read2_path = row['PathReadReverse']
                adapter_seq1 = row['SampleBarcodeForward']
                adapter_seq2 = row['SampleBarcodeReverse']

                lowercase_read1_path = os.path.join(output_path,
                                                    "lowercase_R1.fastq.gz")  # Fixed lowercase filename for R1
                lowercase_read2_path = os.path.join(output_path,
                                                    "lowercase_R2.fastq.gz")  # Fixed lowercase filename for R2

                print(f"Applying lowercase action to paired-end reads for sample: {row['Sample']} to {output_path}")
                run_cutadapt_paired_end_lowercase(read1_path, read2_path, output_path, adapter_seq1, adapter_seq2,
                                                  min_len, threads)

                # Store the paths of the lowercase reads in self.infofile
                self.infofile.at[index, 'PathReadForwardLowercase'] = lowercase_read1_path
                self.infofile.at[index, 'PathReadReverseLowercase'] = lowercase_read2_path

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







import os
import subprocess

def extract_umi(
    fastq1_path, 
    fastq2_path=None, 
    adapter1=None, 
    adapter2=None, 
    output_path=None, 
    umi_length=8, 
    threads=1
):
    """
    Extract UMIs from single-end or paired-end FASTQ files using umi_tools.
    
    Args:
    - fastq1_path (str): Path to the first FASTQ file (for single-end or paired-end reads).
    - fastq2_path (str, optional): Path to the second FASTQ file (only for paired-end reads).
    - adapter1 (str): Adapter sequence for --bc-pattern for single-end reads or the first read in paired-end reads.
    - adapter2 (str, optional): Adapter sequence for the second read in paired-end reads (if applicable).
    - output_path (str): Path to write the output FASTQ file(s).
    - umi_length (int): Length of the UMI sequence. Defaults to 8.
    - threads (int): Number of threads to use for the extraction. Defaults to 1.
    
    Returns:
    - None
    """

    # Check if all required arguments are provided
    if not fastq1_path or not adapter1 or not output_path:
        raise ValueError("You must provide at least the input fastq path, adapter, and output path.")

    # Construct UMI pattern for single-end or paired-end
    bc_pattern1 = f"(?P<umi_1>.{{{umi_length}}}){adapter1}"
    
    # Construct the umi_tools extract command for single-end or paired-end reads
    if fastq2_path and adapter2:
        # Paired-end reads
        bc_pattern2 = f"(?P<umi_2>.{{{umi_length}}}){adapter2}"
        cmd = [
            "umi_tools", "extract",
            "--bc-pattern", bc_pattern1,
            "--bc-pattern2", bc_pattern2,
            "-I", fastq1_path,
            "--read2-in", fastq2_path,
            "-S", os.path.join(output_path, "output_R1.fastq.gz"),
            "--read2-out", os.path.join(output_path, "output_R2.fastq.gz"),
            "--extract-method=regex",
            "--threads", str(threads)
        ]
    else:
        # Single-end reads
        cmd = [
            "umi_tools", "extract",
            "--bc-pattern", bc_pattern1,
            "-I", fastq1_path,
            "-S", os.path.join(output_path, "output.fastq.gz"),
            "--extract-method=regex",
            "--threads", str(threads)
        ]
    
    # Run the command
    try:
        subprocess.run(cmd, check=True)
        print("UMI extraction completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during UMI extraction: {e}")
        raise

# Example usage for single-end
# extract_umi(fastq1_path="path/to/single_end.fastq.gz", adapter1="CGTGTGAG", output_path="path/to/output", umi_length=8, threads=4)

# Example usage for paired-end
# extract_umi(fastq1_path="path/to/read1.fastq.gz", fastq2_path="path/to/read2.fastq.gz", adapter1="CGTGTGAG", adapter2="AGTGTGCG", output_path="path/to/output", umi_length=8, threads=4)
