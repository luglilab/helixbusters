from helixbusters.core import Helixbusters
helix = Helixbusters(samplesheet="/home/lugli/spuccio/Projects/SP036_Lise/Dev/samplesheet_tp2.csv", species="human",
mismatch=1,genome_index="/mnt/references/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa",
)
helix.read_column_from_excel()
helix.infofile
helix.modality
helix.create_sample_output_folders("/home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/")
helix.process_infofile(umi_length=8, threads=10)
helix.run_bwa_mapping(quality=20, threads=10)
helix.generate_umi_output_for_samples()


helix.infofile['PathReadForward'][0]
helix.infofile['OutputPath'][0]

 #/home/lugli/spuccio/Projects/SP036_Lise/RITM0026254_BLISS_VL_EG/test_sample_pipeline/subSPE_013_LIB_S1_L001_R1_001.fastq.gz