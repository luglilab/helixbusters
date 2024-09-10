from helixbusters.core import Helixbusters
helix = Helixbusters(samplesheet="/home/lugli/spuccio/Projects/SP036_Lise/Dev/samplesheet_tp.csv", species="human",
mismatch=1,genome_index="/mnt/references/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa",
)
helix.read_column_from_excel()
helix.infofile
helix.modality
helix.create_sample_output_folders("/home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/")
helix.trim_reads("/home/lugli/spuccio/Projects/SP036_Lise/Dev/TestPipe/")
helix.run_bwa_mapping(quality=20, threads=10)


helix.infofile['BamFilteredPath'][0]