echo "Install conda before"
echo "Download test data"
wget -P ./test_data ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip ./test_data/human_g1k_v37.fasta.gz

echo "Configure and run environment"
conda env create -f environment.yml -yes
conda activate consensus_seq_creator

echo "Help run"
python ./src/ConsensusSeqCreator.py -h