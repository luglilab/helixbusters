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

def run_cutadapt_single_end(read_path, output_path, adapter_seq=None, min_len=20, threads=10):
    """
    Run the cutadapt command to trim single-end sequencing reads.

    Args:
        read_path (str): Path to the input fastq file.
        output_path (str): Path to the directory where output files will be saved.
        adapter_seq (str): Adapter sequence to be trimmed (can be None).
        min_len (int): Minimum length of reads after trimming.
        threads (int): Number of threads to use for processing.
    """

    # Ensure that min_len is an integer
    #min_len = int(min_len)

    # Base command for cutadapt
    cutadapt_cmd = ["cutadapt", "-m", '10', "-j", '10']

    # Add adapter sequence if provided
    if adapter_seq:
        cutadapt_cmd.extend(["-b", adapter_seq])

    # Add output paths and input file
    cutadapt_cmd.extend(
        [f"--untrimmed-output={output_path}/untrimmed.fastq.gz", "-o", f"{output_path}/trimmed.fastq.gz", read_path])

    # Run the cutadapt command
    try:
        subprocess.run(cutadapt_cmd, check = True)
        print(f"Single-end trimming completed. Trimmed reads saved to {output_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running cutadapt: {e}")
#todo
def run_cutadapt_paired_end(read1_path, read2_path, output_path, adapter_seq1=None, adapter_seq2=None, min_len=20, threads=10):
    """
    Run the cutadapt command to trim paired-end sequencing reads.

    Args:
        read1_path (str): Path to the forward read (R1) fastq file.
        read2_path (str): Path to the reverse read (R2) fastq file.
        output_path (str): Path to the directory where output files will be saved.
        adapter_seq1 (str): Adapter sequence for forward reads (R1) (can be None).
        adapter_seq2 (str): Adapter sequence for reverse reads (R2) (can be None).
        min_len (int): Minimum length of reads after trimming.
        threads (int): Number of threads to use for processing.
    """

    # Ensure that min_len is an integer
    min_len = int(min_len)

    # Base command for cutadapt
    cutadapt_cmd = ["cutadapt", "-m", str(min_len), "-j", str(threads), "--no-indels"]

    # Add adapter sequences if provided
    if adapter_seq1 and adapter_seq2:
        cutadapt_cmd.extend(["-b", adapter_seq1, "-B", adapter_seq2])
    elif adapter_seq1:
        cutadapt_cmd.extend(["-b", adapter_seq1])
    elif adapter_seq2:
        cutadapt_cmd.extend(["-B", adapter_seq2])

    # Add output paths and input files
    cutadapt_cmd.extend([
        f"--untrimmed-output={output_path}/untrimmed.fastq.gz",
        "-o", f"{output_path}/trimmed_R1.fastq.gz",
        "-p", f"{output_path}/trimmed_R2.fastq.gz",
        read1_path, read2_path
    ])

    # Run the cutadapt command
    try:
        subprocess.run(cutadapt_cmd, check=True)
        print(f"Paired-end trimming completed. Trimmed reads saved to {output_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running cutadapt: {e}")



