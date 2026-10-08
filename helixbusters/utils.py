import os
import subprocess
import pandas as pd
import shutil
import gzip
import pysam

def read_excel_column(file_path):
    """
    Reads a specified column from an Excel or CSV file and performs checks.
    Determines if the modality is single-end or paired-end sequencing.

    Args:
        file_path (str): Path to the file (Excel or CSV).

    Returns:
        tuple: (pandas.DataFrame, str) containing the DataFrame and the modality ("single end" or "paired end")
    """
    required_columns = ['Sample', 'Replicate', 'Group', 'PathReadForward', 'SampleBarcodeForward', 'PathReadReverse',
                        'SampleBarcodeReverse']

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
        missing_columns = [col for col in ['Sample', 'Replicate', 'Group', 'PathReadForward', 'SampleBarcodeForward'] if
                           col not in df.columns]
        if missing_columns:
            raise ValueError(f"The following required columns are missing: {', '.join(missing_columns)}")

        # Determine if it's single-end or paired-end based on column presence
        if 'PathReadReverse' in df.columns and 'SampleBarcodeReverse' in df.columns:
            modality = 'paired-end'
            check_columns = ['Sample', 'PathReadForward', 'SampleBarcodeForward', 'PathReadReverse',
                             'SampleBarcodeReverse']
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


def split_fastq(fastq_path, output_dir, lines_per_chunk=4000000):
    """
    Split a large FASTQ file into smaller chunks for parallel processing.

    Args:
    - fastq_path (str): Path to the input FASTQ file.
    - output_dir (str): Directory where the chunked FASTQ files will be stored.
    - lines_per_chunk (int): Number of lines per chunk (FASTQ files typically have 4 lines per read).

    Returns:
    - List of chunk file paths.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Use the split command to chunk the FASTQ file
    chunk_prefix = os.path.join(output_dir, "chunk_")
    split_cmd = ["split", "-l", str(lines_per_chunk), fastq_path, chunk_prefix]

    subprocess.run(split_cmd, check=True)

    # Return the list of chunked FASTQ files
    return [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith("chunk_")]


def decompress_gz(gz_path, output_path):
    """
    Decompress a gzipped file.

    Args:
    - gz_path (str): Path to the input gzipped file.
    - output_path (str): Path to store the decompressed file.

    Returns:
    - Path to the decompressed file.
    """
    with gzip.open(gz_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return output_path


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
    cutadapt_cmd = ["cutadapt", "-m", str(min_len), "-j", str(threads), "--report", "minimal", "-O", '8']

    if adapter_seq:
        cutadapt_cmd.extend(["-g", adapter_seq])

    cutadapt_cmd.extend(
        [f"--untrimmed-output={output_path}/untrimmed.fastq.gz", "-o", f"{output_path}/trimmed.fastq.gz", read_path])

    try:
        subprocess.run(cutadapt_cmd, check=True)
        print(f"Single-end trimming completed. Trimmed reads saved to {output_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running cutadapt: {e}")


def run_cutadapt_paired_end(read1_path, read2_path, output_path, adapter_seq1=None, adapter_seq2=None, min_len=20,
                            threads=10):
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
    cutadapt_cmd = ["cutadapt", "-m", str(min_len), "-j", str(threads), "--report", "minimal", "-O", '8']

    if adapter_seq1:
        cutadapt_cmd.extend(["-g", adapter_seq1])
    if adapter_seq2:
        cutadapt_cmd.extend(["-G", adapter_seq2])

    cutadapt_cmd.extend([
        f"--untrimmed-output={output_path}/untrimmed_R1.fastq.gz",
        f"--untrimmed-paired-output={output_path}/untrimmed_R2.fastq.gz",
        "-o", f"{output_path}/trimmed_R1.fastq.gz",
        "-p", f"{output_path}/trimmed_R2.fastq.gz",
        read1_path, read2_path
    ])

    try:
        subprocess.run(cutadapt_cmd, check=True)
        print(f"Paired-end trimming completed. Trimmed reads saved to {output_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running cutadapt: {e}")


def extract_umi_parallel(fastq1_path, fastq2_path=None, adapter1=None, adapter2=None,
                         output_path=None, umi_length=8, threads=1, chunk_lines=4000000, sample_barcode=None):
    """
    Extract UMIs from single-end or paired-end FASTQ files using umi_tools with parallel execution.
    Automatically decompress gzipped FASTQ files before processing.

    Args:
    - fastq1_path (str): Path to the first FASTQ file (for single-end or paired-end reads).
    - fastq2_path (str, optional): Path to the second FASTQ file (only for paired-end reads).
    - adapter1 (str): Adapter sequence for --bc-pattern for single-end reads or the first read in paired-end reads.
    - adapter2 (str, optional): Adapter sequence for the second read in paired-end reads (if applicable).
    - output_path (str): Path to write the output FASTQ file(s).
    - umi_length (int): Length of the UMI sequence. Defaults to 8.
    - threads (int): Number of threads to use for parallel UMI extraction. Defaults to 1.
    - chunk_lines (int): Number of lines per chunk when splitting FASTQ files. Defaults to 4 million.
    - sample_barcode (str): Barcode sequence to trim from the final output using cutadapt.

    Returns:
    - None
    """
    # Check if input files are gzipped and decompress them if necessary
    if fastq1_path.endswith('.gz'):
        decompressed_fastq1 = decompress_gz(fastq1_path, fastq1_path[:-3])  # Remove .gz extension
    else:
        decompressed_fastq1 = fastq1_path

    if fastq2_path and fastq2_path.endswith('.gz'):
        decompressed_fastq2 = decompress_gz(fastq2_path, fastq2_path[:-3])
    else:
        decompressed_fastq2 = fastq2_path

    # Step 1: Split the decompressed FASTQ file(s) into chunks
    chunked_files1 = split_fastq(decompressed_fastq1, output_dir=os.path.join(output_path, "chunks1"),
                                 lines_per_chunk=chunk_lines)

    if fastq2_path:
        chunked_files2 = split_fastq(decompressed_fastq2, output_dir=os.path.join(output_path, "chunks2"),
                                     lines_per_chunk=chunk_lines)
    else:
        chunked_files2 = [None] * len(chunked_files1)  # For single-end reads

    # Step 2: Create umi_tools extract commands for each chunk
    commands = []
    processed_files = []  # Store the paths of processed files to concatenate later
    for chunk1, chunk2 in zip(chunked_files1, chunked_files2):
        output_chunk = os.path.join(output_path, f"{os.path.basename(chunk1)}_processed.fastq.gz")
        processed_files.append(output_chunk)
        if chunk2:
            output_chunk_r2 = os.path.join(output_path, f"{os.path.basename(chunk1)}_processed_R2.fastq.gz")
            cmd = [
                "umi_tools", "extract",
                "--bc-pattern", f"(?P<umi_1>.{{{umi_length}}}){adapter1}",
                "--bc-pattern2", f"(?P<umi_2>.{{{umi_length}}}){adapter2}",
                "-I", chunk1,
                "--read2-in", chunk2,
                "-S", output_chunk,
                "--read2-out", output_chunk_r2,
                "--extract-method=regex",
                "--verbose=0"
            ]
        else:
            cmd = [
                "umi_tools", "extract",
                "--bc-pattern", f"(?P<umi_1>.{{{umi_length}}}){adapter1}",
                "-I", chunk1,
                "-S", output_chunk,
                "--extract-method=regex",
                "--verbose=0"
            ]
        commands.append(cmd)

    # Step 3: Run the commands in parallel
    try:
        processes = []
        for cmd in commands:
            # Run each umi_tools extract command in parallel
            p = subprocess.Popen(cmd)
            processes.append(p)
            if len(processes) >= threads:  # Limit number of parallel threads
                # Wait for all threads to finish before launching more
                for p in processes:
                    p.wait()
                processes = []

        # Wait for any remaining processes to finish
        for p in processes:
            p.wait()

        print("UMI extraction completed successfully.")

        # Step 4: Concatenate all processed chunks into a single file
        final_output_path = os.path.join(output_path, "final_output.fastq.gz")
        if fastq2_path:
            final_output_r2 = os.path.join(output_path, "final_output_R2.fastq.gz")

        with open(final_output_path, 'wb') as final_output:
            for processed_file in processed_files:
                with open(processed_file, 'rb') as chunk_file:
                    shutil.copyfileobj(chunk_file, final_output)

        print(f"Chunks concatenated into: {final_output_path}")

        # Step 5: Run cutadapt to trim SampleBarcodeForward from the final output
        if sample_barcode:
            if fastq2_path:
                # Paired-end trimming
                run_cutadapt_paired_end(final_output_path, final_output_r2, output_path, adapter_seq1=sample_barcode)
            else:
                # Single-end trimming
                run_cutadapt_single_end(final_output_path, output_path, adapter_seq=sample_barcode, threads=threads)

            print(f"Trimming with cutadapt completed using SampleBarcodeForward {sample_barcode}.")

    except subprocess.CalledProcessError as e:
        print(f"Error during UMI extraction: {e}")
        raise

    finally:
        # Clean up chunk files after processing
        shutil.rmtree(os.path.join(output_path, "chunks1"))
        if fastq2_path:
            shutil.rmtree(os.path.join(output_path, "chunks2"))

        # Optionally remove processed files
        for processed_file in processed_files:
            os.remove(processed_file)

        print("Temporary chunk files and processed files removed.")


def process_bam_and_generate_umi_outputs(input_bam, output_file_umi_pcr,
                                         output_file_umi_count, **options):
    """Deduplicate BLISS ends; see helixbusters.deduplication.deduplicate_bam.

    The existing three positional arguments remain supported. The count BED
    now contains one row per position, summing molecules across strands.
    """
    from helixbusters.deduplication import deduplicate_bam

    return deduplicate_bam(input_bam, output_file_umi_pcr, output_file_umi_count, **options)


def plot_alignment_quality(bam_file, output_plot_path):
    """
    Generate a plot for the distribution of alignment quality from a BAM file.

    Args:
        bam_file (str): Path to the input BAM file.
        output_plot_path (str): Path where the plot will be saved.

    Returns:
        None
    """
    import matplotlib.pyplot as plt

    # Open the BAM file
    alignment_qualities = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                alignment_qualities.append(read.mapping_quality)

    # Create the plot for alignment quality distribution
    plt.figure(figsize=(10, 6))
    plt.hist(alignment_qualities, bins=50, color='blue', alpha=0.7)
    plt.title('Distribution of Alignment Quality')
    plt.xlabel('Mapping Quality')
    plt.ylabel('Frequency')
    plt.grid(True)

    # Save the plot to a file
    plt.savefig(output_plot_path)
    plt.close()

    print(f"Plot saved to {output_plot_path}")
