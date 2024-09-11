from helixbusters.utils import (
    read_excel_column,
    extract_umi
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
        Iterates over the rows of self.infofile, performing UMI extraction for each sample.

        Args:
        - umi_length (int): Length of the UMI sequence. Defaults to 8.
        - threads (int): Number of threads to use for UMI extraction. Defaults to 1.
        """
        for index, row in self.infofile.iterrows():
            fastq1_path = row['PathReadForward']
            fastq2_path = row.get('PathReadReverse', None)  # Use PathReadReverse if available (for paired-end)
            barcode_forward = row['SampleBarcodeForward']  # Now using SampleBarcodeForward instead of AdapterForward
            barcode_reverse = row.get('SampleBarcodeReverse',
                                      None)  # Now using SampleBarcodeReverse instead of AdapterReverse
            output_path = row['OutputPath']

            # Ensure output directory exists
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            print(f"Processing sample {row['Sample']} (row {index + 1}/{len(self.infofile)})")

            # Call the extract_umi function from helixbusters.utils with the appropriate arguments
            extract_umi(
                fastq1_path=fastq1_path,
                fastq2_path=fastq2_path,
                adapter1=barcode_forward,
                adapter2=barcode_reverse,
                output_path=output_path,
                umi_length=umi_length,
                threads=threads
            )

