from helixbusters.utils import read_excel_column
import os
import subprocess
import sys
import gzip
from Bio import SeqIO


class Helixbusters:
    def __init__(self, samplesheet, species, mismatch):
        self.samplesheet = samplesheet
        self.species = species
        self.mismatch = mismatch

        # Validate species attribute
        if species.lower() not in ['mouse', 'human']:
            raise ValueError("Species must be 'mouse' or 'human'.")

        # Validate mismatch attribute
        if not 0 <= mismatch <= 3:
            raise ValueError("Mismatch must be between 0 and 3.")

        self.mismatch = mismatch

    def read_column_from_excel(self):
        """
        Reads and returns the content of the samplesheet using the read_excel_column function from utils.
        """
        column_data = read_excel_column(self.samplesheet)
        return column_data

    def process_samples(self):
        """
        Processes each row from the samplesheet and generates an output string.
        The output string is saved in a new column in the DataFrame.

        Returns:
            pd.DataFrame: Updated DataFrame with the new column containing the output string.
        """
        df = self.read_column_from_excel()

        # Check if DataFrame is empty (which would indicate an error occurred in reading)
        if df.empty:
            raise ValueError("Samplesheet could not be read properly or is empty.")

        # Initialize a new column in the DataFrame to hold the output strings
        df['OutputString'] = ''

        for index, row in df.iterrows():
            # Extract fields from the row
            f1 = row['Sample']  # First field (Sample)
            f2 = row['Sample'].replace('_', '').replace('/', '')  # Remove underscores and slashes
            f3 = row['Group']  # Third field (Group)
            f6 = row['SampleBarcode1']  # Fifth field (SampleBarcode1)

            # Create the output string
            output_string = f'^ 8...8 {f3}[{self.mismatch},0,0] 1...1000 $'

            # Add the output string to the new column
            df.at[index, 'OutputString'] = output_string

        # Return the modified DataFrame with the new column
        return df

    def decompress_and_process_fastq(self, updated_df, output_dir):
        """
        Decompress and process the FASTQ files using information from the DataFrame (updated_df),
        utilizing Python libraries for file processing.

        Args:
            updated_df (pd.DataFrame): DataFrame containing sample information and file paths.
            output_dir (str): Path to the output directory.
        """
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Loop through each row of the DataFrame
        for index, row in updated_df.iterrows():
            # Extract required information for each sample from the DataFrame
            sample = row['Sample']
            path_read_f = row['PathReadF']  # Path to the r1 fastq file
            path_read_r2 = row.get('PathReadR2')  # Assuming PathReadR2 is in the DataFrame for paired-end data
            output_string = row['OutputString']  # This is the output string generated earlier

            # Paths for intermediate files
            r1_fa = os.path.join(output_dir, f'{sample}_r1.fa')
            r1_oneline_fq = os.path.join(output_dir, f'{sample}_r1oneline.fq')
            r2_oneline_fq = os.path.join(output_dir, f'{sample}_r2oneline.fq') if path_read_r2 else None

            print(f"Processing sample {sample}...")

            # Process r1 (forward read)
            self.process_fastq_file(path_read_f, r1_fa, r1_oneline_fq)

            # If paired-end (r2 available), process r2 as well
            if path_read_r2:
                self.process_fastq_file(path_read_r2, None, r2_oneline_fq)

            print(f"Done processing sample {sample} FASTQ files.")

    def process_fastq_file(self, input_path, output_fa_path=None, output_oneline_fq_path=None):
        """
        Process a single FASTQ file. Decompress it and perform transformations.
        If output paths are provided, save the transformed files in FASTA and oneline formats.

        Args:
            input_path (str): Path to the input FASTQ file.
            output_fa_path (str): Path to the output FASTA file (optional).
            output_oneline_fq_path (str): Path to the sorted one-line FASTQ file (optional).
        """
        # Open the gzipped FASTQ file and read sequences
        with gzip.open(input_path, 'rt') as handle:
            records = list(SeqIO.parse(handle, 'fastq'))

        # If output_fa_path is provided, write to FASTA format
        if output_fa_path:
            with open(output_fa_path, 'w') as fa_out:
                for record in records:
                    # Writing to FASTA with only the first two lines of each record
                    fa_out.write(f'>{record.id}\n{str(record.seq)}\n')

        # If output_oneline_fq_path is provided, write the FASTQ file in oneline format and sort by sequence ID
        if output_oneline_fq_path:
            sorted_records = sorted(records, key=lambda rec: rec.id)  # Sort by sequence ID

            with open(output_oneline_fq_path, 'w') as fq_out:
                for record in sorted_records:
                    fq_out.write(record.format('fastq'))

    def run_command(self, command):
        """
        Utility method to run shell commands.

        Args:
            command (str): Command to run as a string.

        Returns:
            None
        """
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running command: {command}\n{e}")
            raise
