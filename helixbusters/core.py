from helixbusters.utils import read_excel_column
import os
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
            f1 = row['Sample']  # First field (Sample)
            f2 = row['Sample'].replace('_', '').replace('/', '')  # Remove underscores and slashes
            f3 = row['Group']  # Third field (Group)
            f6 = row['SampleBarcode1']  # Fifth field (SampleBarcode1)

            output_string = f'^ 8...8 {f3}[{self.mismatch},0,0] 1...1000 $'
            df.at[index, 'OutputString'] = output_string

        return df

    def reverse_complement(self, seq):
        """
        Returns the reverse complement of the DNA sequence.
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in seq[::-1])

    def fastq_to_fasta_file(self, input_fastq_gz, adapter_sequence, output_fasta_file, min_rest_length=20):
        """
        Process a FASTQ file, filter sequences by adapter and length, and output in FASTA format.
        """
        with gzip.open(input_fastq_gz, 'rt') as fastq_file, open(output_fasta_file, 'w') as fasta_file:
            while True:
                header = fastq_file.readline().strip()
                if not header:
                    break  # End of file
                sequence = fastq_file.readline().strip()
                fastq_file.readline()  # Skip '+'
                fastq_file.readline()  # Skip quality scores

                header_parts = header.split(' ')
                read_id = header_parts[0][1:]  # Remove '@'

                if adapter_sequence not in sequence:
                    continue  # Skip if adapter not found

                first_chunk = sequence[:8]
                second_chunk = self.reverse_complement(sequence[8:16])
                rest_of_sequence = sequence[16:]

                if len(rest_of_sequence) < min_rest_length:
                    continue  # Skip if too short

                fasta_header = f">{read_id}:[1,{len(sequence)}]"
                formatted_sequence = f"{first_chunk} {second_chunk} {rest_of_sequence}"
                fasta_file.write(f"{fasta_header}\n{formatted_sequence}\n")

    def decompress_and_process_fastq(self, updated_df, output_dir, min_rest_length=20):
        """
        Decompress and process the FASTQ files using the DataFrame (updated_df),
        and update the DataFrame with paths to the created FASTA files.
        """
        os.makedirs(output_dir, exist_ok=True)
        updated_df['FASTA_Path'] = ''

        for index, row in updated_df.iterrows():
            sample = row['Sample']
            path_read_f = row['PathReadF']  # Path to the r1 fastq file
            sample_barcode1 = row['SampleBarcode1']  # Adapter sequence
            output_string = row['OutputString']

            r1_fa = os.path.join(output_dir, f'{sample}_r1.fa')
            print(f"Processing sample {sample}...")

            self.fastq_to_fasta_file(path_read_f, sample_barcode1, r1_fa, min_rest_length)
            updated_df.at[index, 'FASTA_Path'] = r1_fa

            print(f"Done processing and filtering sample {sample} FASTA files.")

        return updated_df

    def prepare_for_mapping(self, updated_df, output_dir):
        """
        This function processes the FASTA files and joins them with corresponding FASTQ files.
        """
        print('Parse the FASTA files, filtering and trimming ...')
        os.makedirs(output_dir, exist_ok=True)

        for index, row in updated_df.iterrows():
            sample = row['Sample']
            fasta_path = row['FASTA_Path']

            r1_output_path = os.path.join(output_dir, f'{sample}_r1.2b.aln.fq')
            r2_output_path = os.path.join(output_dir, f'{sample}_r2.2b.aln.fq') if 'r2oneline.fq' in updated_df.columns else None

            with open(fasta_path, 'r') as fasta_file:
                fasta_lines = fasta_file.readlines()

            id_genomic = {}
            for i in range(0, len(fasta_lines), 2):
                header = fasta_lines[i].strip().split(':')[:7]
                sequence = fasta_lines[i + 1].strip()
                id_genomic['@'.join(header)] = sequence

            id_genomic_sorted = dict(sorted(id_genomic.items()))
            r1oneline_path = os.path.join(output_dir, 'r1oneline.fq')
            r1_sequences = self.load_fastq_as_dict(r1oneline_path)

            r1_joined = self.join_sequences(id_genomic_sorted, r1_sequences, is_paired_end=False)
            with open(r1_output_path, 'w') as r1_output_file:
                for line in r1_joined:
                    r1_output_file.write(line + '\n')

            if r2_output_path:
                r2oneline_path = os.path.join(output_dir, 'r2oneline.fq')
                r2_sequences = self.load_fastq_as_dict(r2oneline_path)

                r2_joined = self.join_sequences(id_genomic_sorted, r2_sequences, is_paired_end=True)
                with open(r2_output_path, 'w') as r2_output_file:
                    for line in r2_joined:
                        r2_output_file.write(line + '\n')

        print('Done! Ready to be aligned to the reference genome!')

    def load_fastq_as_dict(self, fastq_file_path):
        """
        Load a FASTQ file into a dictionary with sequence IDs as keys and sequences as values.
        """
        sequences = {}
        with open(fastq_file_path, 'r') as fastq_file:
            while True:
                header = fastq_file.readline().strip()
                if not header:
                    break  # End of file
                sequence = fastq_file.readline().strip()
                fastq_file.readline()  # Skip '+'
                quality = fastq_file.readline().strip()

                sequence_id = header.split()[0]
                sequences[sequence_id] = (sequence, quality)

        return sequences

    def join_sequences(self, id_genomic_sorted, fastq_sequences, is_paired_end=False):
        """
        Join the sorted genomic IDs with sequences from the fastq_sequences dictionary.
        """
        joined_sequences = []

        for seq_id, genomic_sequence in id_genomic_sorted.items():
            if seq_id in fastq_sequences:
                original_sequence, quality = fastq_sequences[seq_id]

                if is_paired_end:
                    joined_sequences.append(f"{seq_id}\n{genomic_sequence}")
                else:
                    trimmed_sequence = genomic_sequence[-len(original_sequence):]
                    joined_sequences.append(f"{seq_id}\n{original_sequence}\n+\n{trimmed_sequence}")

        return joined_sequences
