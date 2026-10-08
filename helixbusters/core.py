from helixbusters.utils import (
    read_excel_column,
    extract_umi_parallel,
    run_cutadapt_single_end,
    run_cutadapt_paired_end,
    process_bam_and_generate_umi_outputs,
    plot_alignment_quality
)
import os
import pandas as pd
from helixbusters.mapping import map_sample, validate_sample_name
from helixbusters.genomes import (normalize_genome, genome_species, make_genome_filter,
                                  load_reference_config)


class Helixbusters:
    def __init__(self, samplesheet, species=None, mismatch=1, genome_index=None, *,
                 genome=None, blacklist_bed=None, blacklist_genome=None):
        self.samplesheet = samplesheet
        self.genome = normalize_genome(genome) if genome is not None else None
        if species is None and self.genome is not None:
            species = genome_species(self.genome)
        self.species = species.lower() if isinstance(species, str) else species
        if self.genome is not None and self.species != genome_species(self.genome):
            raise ValueError("Species does not match the selected genome")
        self.blacklist_bed = blacklist_bed
        self.blacklist_genome = blacklist_genome
        make_genome_filter(self.genome, blacklist_bed, blacklist_genome)
        self.reference_aligner = None
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

    @classmethod
    def from_reference_config(cls, samplesheet, genome, reference_config, *,
                              aligner="bwa", mismatch=1):
        """Choose a build and its index/blacklist from a local reference catalog."""
        reference = load_reference_config(reference_config, genome, aligner)
        instance = cls(samplesheet, mismatch=mismatch, **reference)
        instance.reference_aligner = aligner
        return instance

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

    def run_bwa_mapping(self, quality=20, threads=10, **options):
        """Map trimmed reads with classic BWA-MEM; quality is minimum MAPQ.

        Additional options: sort_threads, sort_memory, aligner_executable,
        samtools_executable. See docs/mapping.md for filtering and QC outputs.
        """
        self._run_mapping("bwa", quality, threads, **options)

    def run_bowtie2_mapping(self, quality=20, threads=10, **options):
        """Map trimmed reads with Bowtie2 (end-to-end by default).

        bowtie2_mode="local" is an explicit alternative. quality is MAPQ,
        not the per-base FASTQ quality. Other options match run_bwa_mapping.
        """
        self._run_mapping("bowtie2", quality, threads, **options)

    def _run_mapping(self, aligner, quality, threads, **options):
        if self.reference_aligner is not None and self.reference_aligner != aligner:
            raise ValueError(f"Reference config selected a {self.reference_aligner} index, not {aligner}")
        if self.genome_index is None:
            raise ValueError("A genome_index prefix is required for mapping")
        if {"genome", "blacklist_bed", "blacklist_genome"} & options.keys():
            raise ValueError("Configure genome and blacklist on Helixbusters, not on the mapping method")
        if self.infofile is None or self.modality not in {"single-end", "paired-end"}:
            raise ValueError("Load sample information and set a valid sequencing modality first")
        if not {"Sample", "OutputPath"}.issubset(self.infofile.columns):
            raise ValueError("Sample information must contain Sample and OutputPath")
        if self.infofile["Sample"].duplicated().any():
            raise ValueError("Sample names must be unique to prevent output overwrites")
        jobs = []
        for index, row in self.infofile.iterrows():
            sample = row["Sample"]
            validate_sample_name(sample)
            output = row["OutputPath"]
            if pd.isna(output) or not output:
                raise ValueError(f"Missing OutputPath for sample {sample}")
            read1 = row.get("PathReadForwardTrimmed")
            if read1 is None or pd.isna(read1) or not read1:
                name = "trimmed.fastq.gz" if self.modality == "single-end" else "trimmed_R1.fastq.gz"
                read1 = os.path.join(output, name)
            read2 = None
            if self.modality == "paired-end":
                read2 = row.get("PathReadReverseTrimmed")
                if read2 is None or pd.isna(read2) or not read2:
                    read2 = os.path.join(output, "trimmed_R2.fastq.gz")
            for path in (read1, read2):
                if path is not None and not os.path.isfile(path):
                    raise FileNotFoundError(f"Trimmed FASTQ file not found for sample {sample}: {path}")
            jobs.append((index, sample, output, read1, read2))

        # All sample FASTQ paths have been checked before the first mapping job.
        for index, sample, output, read1, read2 in jobs:
            print(f"Running {aligner} mapping for sample {sample}...")
            outputs = map_sample(sample, read1, self.genome_index, output,
                                 read2=read2, aligner=aligner, min_mapq=quality,
                                 threads=threads, genome=self.genome,
                                 blacklist_bed=self.blacklist_bed,
                                 blacklist_genome=self.blacklist_genome, **options)
            for column, path in outputs.items():
                self.infofile.at[index, column] = path

    def generate_umi_output_for_samples(self, method="exact", umi_length=8,
                                        max_distance=1, min_mapq=0,
                                        five_prime_policy="strict",
                                        read_selection="single-end"):
        """
        Count BLISS molecules from each coordinate-sorted BamFilteredPath.

        Exact grouping is the baseline; directional corrects related UMIs at
        the same position and strand. Paired-end BAMs require explicit read1
        or read2 selection according to the library design. Also writes a
        molecule BED6, stranded site table and deduplication QC JSON.
        """
        if self.infofile is None:
            raise ValueError("Load sample information before generating UMI outputs")
        required = {"Sample", "OutputPath", "BamFilteredPath"}
        if not required.issubset(self.infofile.columns):
            raise ValueError("Sample information must contain Sample, OutputPath and BamFilteredPath")
        if self.infofile["Sample"].duplicated().any():
            raise ValueError("Sample names must be unique to prevent output overwrites")
        # Check all BAM paths before producing any sample outputs.
        for _, row in self.infofile.iterrows():
            path = row["BamFilteredPath"]
            if pd.isna(path) or not path or not os.path.isfile(path):
                raise FileNotFoundError(f"BAM file not found for sample {row['Sample']}: {path}")
        for index, row in self.infofile.iterrows():
            sample = row['Sample']
            bam_filtered_path = row['BamFilteredPath']

            # Define the output file paths
            output_file_umi_pcr = os.path.join(row['OutputPath'], f"{sample}_Chromosome-Location-Strand-UMI-PCR.txt")
            output_file_umi_count = os.path.join(row['OutputPath'], f"{sample}_Chromosome-Location-UMI-Count.bed")
            output_molecules = os.path.join(row['OutputPath'], f"{sample}_molecules.bed")
            output_sites = os.path.join(row['OutputPath'], f"{sample}_sites.tsv")
            output_qc = os.path.join(row['OutputPath'], f"{sample}_deduplication.json")

            # Call the utility function to generate UMI outputs
            print(f"Generating UMI output files for sample {sample}...")
            process_bam_and_generate_umi_outputs(
                bam_filtered_path, output_file_umi_pcr, output_file_umi_count,
                method=method, umi_length=umi_length, max_distance=max_distance,
                min_mapq=min_mapq, five_prime_policy=five_prime_policy,
                read_selection=read_selection, output_molecules=output_molecules,
                output_sites=output_sites, output_qc=output_qc,
            )

            # Optionally, store the paths to the generated files in the infofile
            self.infofile.at[index, 'UMI_PCR_Output'] = output_file_umi_pcr
            self.infofile.at[index, 'UMI_Count_Output'] = output_file_umi_count
            self.infofile.at[index, 'Molecules_Output'] = output_molecules
            self.infofile.at[index, 'Sites_Output'] = output_sites
            self.infofile.at[index, 'Deduplication_QC'] = output_qc

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
