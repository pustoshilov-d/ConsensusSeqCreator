# ConsensusSeqCreator

Tool for creating fasta-consensuses (parts of the reference genome with mutations from multi-sample VCF).

*Python 3.10 and pandas are required.*


#### Configure
After installing Conda run `sh configure.sh` 
or do the steps below 

#### 1. Download test data
- `wget -P ./test_data ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz`
- `gunzip ./test_data/human_g1k_v37.fasta.gz`

#### 2. Create and run environment
- `conda env create -f environment.yml -yes`
- `conda activate consensus_seq_creator`

#### 3. Read help
- `python ./src/ConsensusSeqCreator.py -h`

#### 4. Run test
- `python ./src/ConsensusSeqCreator.py -vcf ./test_data/mini_test.vcf -fasta ./test_data/human_g1k_v37.fasta -len 20 -num 10`

Result will be in `./output` or specified in -out directory. 

#### Parameters
##### Required:
- **-vcf --vcf_path** vcs file
- **-fasta --fasta_path** reference genome fasta file
- **-len --consensus_len** length of consensuses (default `70`)
- **-num --consensus_num** number of consensuses (default `10`)

##### Optional:
- **-out -output_dir** output directory (default `output`)
- **-a --allele** allele index (default `1`, only 1 is tested)
- **-s --sample** sample name to use (default fist sample)
- **--mutation_place** place of mutation in consensus (default `center`, {start, center, end})
- **--filter_af_lower** lower value of AF (default `None`)
- **--filter_af_upper** upper value of AF (default `None`)
- **--filter_af_none_include** include or not snp with AF=Non (default `False`)
- **--filter_alt1_only** if few ALT, use only first one (default `True`, only True is tested)
- **--filter_homo** use only homologous alleles (default `True`, only True tested)
- **--filter_pass** only FILTER=PASS (default `True`)
