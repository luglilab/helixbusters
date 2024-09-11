import os
import subprocess
import pandas as pd

def read_excel_column(file_path):
    """
    Reads a specified column from an Excel or CSV file and performs checks.
    Determines if the modality is single-end or paired-end sequencing.

    Args:
        file_path (str): Path to the file (Excel or CSV).

    Returns:
        tuple: (pandas.DataFrame, str) containing the DataFrame and the modality ("single end" or "paired end")
    """
    required_columns = ['Sample', 'Replicate', 'Group', 'PathReadForward', 'SampleBarcodeForward', 'PathReadReverse', 'SampleBarcodeReverse']

    try:
        # Determine file type based on extension
        _, file_extension = os.path.splitext(file_path)

        if file_extension.lower() == '.xlsx' or file_extension.lower() == '.xls':
            # Read Excel file
            df = pd.read_excel(file_path, header=0)
        elif file_extension.lower() == '.csv':
            # Read CSV file
            df = pd.read_csv(file_path, header=0)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}. Only .xlsx, .xls, and .csv are supported.")

        # Check if all required columns are present
        missing_columns = [col for col in ['Sample', 'Replicate', 'Group', 'PathReadForward', 'SampleBarcodeForward'] if col not in df.columns]
        if missing_columns:
            raise ValueError(f"The following required columns are missing: {', '.join(missing_columns)}")

        # Determine if it's single-end or paired-end based on column presence
        if 'PathReadReverse' in df.columns and 'SampleBarcodeReverse' in df.columns:
            modality = 'paired-end'
            check_columns = ['Sample', 'PathReadForward', 'SampleBarcodeForward', 'PathReadReverse', 'SampleBarcodeReverse']
        else:
            modality = 'single-end'
            check_columns = ['Sample', 'PathReadForward', 'SampleBarcodeForward']

        # Check if required columns are unique
        for col in check_columns:
            if not df[col].is_unique:
                raise ValueError(f"The values in column '{col}' are not unique.")

        return df, modality

    except Exception as e:
        print(f"An error occurred: {e}")
        return pd.DataFrame(), None


def extract_umi(fastq1_path, fastq2_path=None, adapter1=None, adapter2=None,
                output_path=None, umi_length=8, threads=1):
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
            "--verbose=0"
        ]
    else:
        # Single-end reads
        cmd = [
            "umi_tools", "extract",
            "--bc-pattern", bc_pattern1,
            "-I", fastq1_path,
            "-S", os.path.join(output_path, "output.fastq.gz"),
            "--extract-method=regex",
            "--verbose=0"
        ]

    # Run the command
    try:
        subprocess.run(cmd, check=True)
        print("UMI extraction completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during UMI extraction: {e}")
        raise

#
# def run_cutadapt_single_end(read_path, output_path, adapter_seq=None, min_len=20, threads=10):
#     """
#     Run the cutadapt command to trim single-end sequencing reads.
#
#     Args:
#         read_path (str): Path to the input fastq file.
#         output_path (str): Path to the directory where output files will be saved.
#         adapter_seq (str): Adapter sequence to be trimmed (can be None).
#         min_len (int): Minimum length of reads after trimming.
#         threads (int): Number of threads to use for processing.
#     """
#
#     # Ensure that min_len is an integer
#     #min_len = int(min_len)
#
#     # Base command for cutadapt
#     cutadapt_cmd = ["cutadapt", "-m", '10', "-j", '10','--report','minimal','-O', '8']
#
#     # Add adapter sequence if provided
#     if adapter_seq:
#         cutadapt_cmd.extend(["-b", adapter_seq])
#
#     # Add output paths and input file
#     cutadapt_cmd.extend(
#         [f"--untrimmed-output={output_path}/untrimmed.fastq.gz", "-o", f"{output_path}/trimmed.fastq.gz", read_path])
#
#     # Run the cutadapt command
#     try:
#         subprocess.run(cutadapt_cmd, check = True)
#         print(f"Single-end trimming completed. Trimmed reads saved to {output_path}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running cutadapt: {e}")
# #todo
# def run_cutadapt_paired_end(read1_path, read2_path, output_path, adapter_seq1=None, adapter_seq2=None, min_len=20, threads=10):
#     """
#     Run the cutadapt command to trim paired-end sequencing reads.
#
#     Args:
#         read1_path (str): Path to the forward read (R1) fastq file.
#         read2_path (str): Path to the reverse read (R2) fastq file.
#         output_path (str): Path to the directory where output files will be saved.
#         adapter_seq1 (str): Adapter sequence for forward reads (R1) (can be None).
#         adapter_seq2 (str): Adapter sequence for reverse reads (R2) (can be None).
#         min_len (int): Minimum length of reads after trimming.
#         threads (int): Number of threads to use for processing.
#     """
#
#     # Ensure that min_len is an integer
#     min_len = int(min_len)
#
#     # Base command for cutadapt
#     cutadapt_cmd = ["cutadapt", "-m", '10', "-j", '10','--report','minimal','-O', '8']
#
#     # Add adapter sequences if provided
#     if adapter_seq1 and adapter_seq2:
#         cutadapt_cmd.extend(["-b", adapter_seq1, "-B", adapter_seq2])
#     elif adapter_seq1:
#         cutadapt_cmd.extend(["-b", adapter_seq1])
#     elif adapter_seq2:
#         cutadapt_cmd.extend(["-B", adapter_seq2])
#
#     # Add output paths and input files
#     cutadapt_cmd.extend([
#         f"--untrimmed-output={output_path}/untrimmed.fastq.gz",
#         "-o", f"{output_path}/trimmed_R1.fastq.gz",
#         "-p", f"{output_path}/trimmed_R2.fastq.gz",
#         read1_path, read2_path
#     ])
#
#     # Run the cutadapt command
#     try:
#         subprocess.run(cutadapt_cmd, check=True)
#         print(f"Paired-end trimming completed. Trimmed reads saved to {output_path}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running cutadapt: {e}")
#
# def run_fastp_paired_end(read1_path, read2_path, output_path, adapter_seq1=None, adapter_seq2=None, min_len=20, threads=10):
#     """
#     Run the fastp command to trim paired-end sequencing reads.
#
#     Args:
#         read1_path (str): Path to the forward read (R1) fastq file.
#         read2_path (str): Path to the reverse read (R2) fastq file.
#         output_path (str): Path to the directory where output files will be saved.
#         adapter_seq1 (str): Adapter sequence for forward reads (R1) (can be None).
#         adapter_seq2 (str): Adapter sequence for reverse reads (R2) (can be None).
#         min_len (int): Minimum length of reads after trimming.
#         threads (int): Number of threads to use for processing.
#     """
#
#     # Base command for fastp
#     fastp_cmd = [
#         "fastp",
#         "--in1", read1_path,  # Forward input read (R1)
#         "--in2", read2_path,  # Reverse input read (R2)
#         "--out1", f"{output_path}/trimmed_R1.fastq.gz",  # Output for trimmed forward read (R1)
#         "--out2", f"{output_path}/trimmed_R2.fastq.gz",  # Output for trimmed reverse read (R2)
#         "--unpaired1", f"{output_path}/untrimmed_R1.fastq.gz",  # Output for untrimmed forward read (R1)
#         "--unpaired2", f"{output_path}/untrimmed_R2.fastq.gz",  # Output for untrimmed reverse read (R2)
#         "--length_required", str(min_len),  # Minimum length of reads after trimming
#         "--thread", str(threads)  # Number of threads to use
#     ]
#
#     # Add adapter sequences if provided
#     if adapter_seq1:
#         fastp_cmd.extend(["--adapter_sequence", adapter_seq1])
#     if adapter_seq2:
#         fastp_cmd.extend(["--adapter_sequence_r2", adapter_seq2])
#
#     # Run the fastp command
#     try:
#         subprocess.run(fastp_cmd, check=True)
#         print(f"Paired-end trimming completed using fastp. Trimmed reads saved to {output_path}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running fastp: {e}")
#
#
# def run_cutadapt_single_end_lowercase(read_path, output_path, adapter_seq=None, min_len=20, threads=10):
#     """
#     Run the cutadapt command to trim single-end sequencing reads with the --action=lowercase option.
#
#     Args:
#         read_path (str): Path to the input fastq file.
#         output_path (str): Path to the directory where output files will be saved.
#         adapter_seq (str): Adapter sequence to be trimmed (can be None).
#         min_len (int): Minimum length of reads after trimming.
#         threads (int): Number of threads to use for processing.
#     """
#
#     # Base command for cutadapt with --action=lowercase
#     cutadapt_cmd = ["cutadapt", "--action=lowercase", "-m", '10', "-j", '10', '--report', 'minimal','-O', '8']
#
#     # Add adapter sequence if provided
#     if adapter_seq:
#         cutadapt_cmd.extend(["-b", adapter_seq])
#
#     # Add output paths and input file, changing names to include lowercase/unlowercase
#     cutadapt_cmd.extend(
#         [f"--untrimmed-output={output_path}/unlowercase.fastq.gz", "-o", f"{output_path}/lowercase.fastq.gz", read_path])
#
#     # Run the cutadapt command
#     try:
#         subprocess.run(cutadapt_cmd, check=True)
#         print(f"Single-end trimming with lowercase action completed. Output saved to {output_path}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running cutadapt: {e}")
#
# def run_cutadapt_paired_end_lowercase(read1_path, read2_path, output_path, adapter_seq1=None, adapter_seq2=None, min_len=20, threads=10):
#     """
#     Run the cutadapt command to trim paired-end sequencing reads with the --action=lowercase option.
#
#     Args:
#         read1_path (str): Path to the forward read (R1) fastq file.
#         read2_path (str): Path to the reverse read (R2) fastq file.
#         output_path (str): Path to the directory where output files will be saved.
#         adapter_seq1 (str): Adapter sequence for forward reads (R1) (can be None).
#         adapter_seq2 (str): Adapter sequence for reverse reads (R2) (can be None).
#         min_len (int): Minimum length of reads after trimming.
#         threads (int): Number of threads to use for processing.
#     """
#
#     # Ensure that min_len is an integer
#     min_len = int(min_len)
#
#     # Base command for cutadapt with --action=lowercase
#     cutadapt_cmd = ["cutadapt","--action=lowercase", "-m", '10', "-j", '10', '--report', 'minimal', '-O', '8']
#
#     # Add adapter sequences if provided
#     if adapter_seq1 and adapter_seq2:
#         cutadapt_cmd.extend(["-b", adapter_seq1, "-B", adapter_seq2])
#     elif adapter_seq1:
#         cutadapt_cmd.extend(["-b", adapter_seq1])
#     elif adapter_seq2:
#         cutadapt_cmd.extend(["-B", adapter_seq2])
#
#     # Add output paths and input files, changing names to include lowercase/unlowercase
#     cutadapt_cmd.extend([
#         f"--untrimmed-output={output_path}/unlowercase.fastq.gz",
#         "-o", f"{output_path}/lowercase_R1.fastq.gz",
#         "-p", f"{output_path}/lowercase_R2.fastq.gz",
#         read1_path, read2_path
#     ])
#
#     # Run the cutadapt command
#     try:
#         subprocess.run(cutadapt_cmd, check=True)
#         print(f"Paired-end trimming with lowercase action completed. Output saved to {output_path}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running cutadapt: {e}")
#
# def run_fastp_single_end(read_path, output_path, adapter_seq=None, min_len=20, threads=10):
#     """
#     Run the fastp command to trim single-end sequencing reads.
#
#     Args:
#         read_path (str): Path to the input fastq file.
#         output_path (str): Path to the directory where output files will be saved.
#         adapter_seq (str): Adapter sequence to be trimmed (can be None).
#         min_len (int): Minimum length of reads after trimming.
#         threads (int): Number of threads to use for processing.
#     """
#
#     # Base command for fastp
#     fastp_cmd = [
#         "fastp",
#         "--in1", read_path,
#         "--out1", f"{output_path}/trimmed.fastq.gz",
#         "--unpaired1", f"{output_path}/untrimmed.fastq.gz",  # Untrimmed output equivalent
#         "--length_required", str(min_len),
#         "--thread", str(threads)
#     ]
#
#     # Add adapter sequence if provided
#     if adapter_seq:
#         fastp_cmd.extend(["--adapter_sequence", adapter_seq])
#
#     # Run the fastp command
#     try:
#         subprocess.run(fastp_cmd, check=True)
#         print(f"Single-end trimming completed using fastp. Trimmed reads saved to {output_path}.")
#     except subprocess.CalledProcessError as e:
#         print(f"Error running fastp: {e}")



