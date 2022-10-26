DESC = '''\
Tool for creating fasta-consensuses (parts of 
the reference genome with mutations from multi-sample VCF).
Author: pustoshilov-d
        '''
__author__ = "pustoshilov-d"

"""
TODO:
    refactoring: вынести все тексты и цирфы как константы
    feature: добавить поддержку gz
    refactoring: в enrich_af и fasta_parse использовать regex
    debug: проверить поддержку гаплоидных 
    feature: дорабоать поддержку фазирования
    feature: выбор методов решателя, какую использовать аллель
    feature: добавить параметр перемешивать/нет replace в .sample() 
    refactoring: в индексировании фасты добавить Exception при выходе за границы
    refactoring: ввести логирование
    refactoring: ввести нормальный набор Exceptions
"""

from utils import (
    enrich_gt,
    enrich_af,
    filter_homo,
    filter_alt1_only,
    solve_alt,
    mutate_seq,
    save_consensus,
    prepare_dir,
)
from consts import (
    FIRST_SAMPLE_COLUMN
)
import pandas as pd


class ConsensusSeqCreator():
    vcf_path: str
    fasta_path: str
    output_dir: str

    consensus_len: int
    consensus_num: int
    allele: int
    mutation_place: str
    sample: str

    filter_pass: bool
    filter_homo: bool
    filter_alt1_only: bool
    filter_af_lower: float
    filter_af_upper: float
    filter_af_none_include: bool

    snp_df: pd.DataFrame
    snp_full_df: pd.DataFrame

    def __init__(
            self, vcf_path: str, fasta_path: str, output_dir: str, sample: str, mutation_place: str,
            consensus_len: int, consensus_num: int, allele: int, filter_pass: bool, filter_alt1_only: bool,
            filter_homo: bool, filter_af_lower: float, filter_af_upper: float, filter_af_none_include: bool,
            other_args
    ) -> None:
        self.vcf_path = vcf_path
        self.fasta_path = fasta_path
        self.output_dir = output_dir
        self.sample = sample

        self.consensus_len = consensus_len
        self.consensus_num = consensus_num
        self.allele = allele
        self.mutation_place = mutation_place

        self.filter_pass = filter_pass
        self.filter_alt1_only = filter_alt1_only
        self.filter_homo = filter_homo
        self.filter_af_lower = filter_af_lower
        self.filter_af_upper = filter_af_upper
        self.filter_af_none_include = filter_af_none_include

        self.parse_vcf()
        if not self.sample or self.sample not in self.snp_df.columns:
            self.sample = self.snp_df.columns[FIRST_SAMPLE_COLUMN]

        self.enrich_snp_data()
        self.filter_snp()
        prepare_dir(self.output_dir)
        self.split_fasta()

    def parse_vcf(self) -> None:
        snp_df = pd.DataFrame()

        with open(self.vcf_path, "r") as snp_file:
            for line in snp_file:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        snp_df = pd.DataFrame(columns=line[1:].split("\t"))
                    continue

                snp_df.loc[len(snp_df)] = line.split("\t")

        self.snp_df = snp_df

    def enrich_snp_data(self) -> None:
        self.snp_df["AF"] = self.snp_df.apply(enrich_af, axis=1)
        self.snp_df["GT"] = self.snp_df.apply(
            enrich_gt, axis=1, sample=self.sample)

        self.snp_full_df = self.snp_df

    def filter_snp(self) -> None:
        if self.filter_pass:
            self.snp_df = self.snp_df[self.snp_df["FILTER"] == "PASS"]

        if self.filter_homo:
            self.snp_df = self.snp_df[self.snp_df["GT"].apply(filter_homo)]

        if self.filter_alt1_only:
            self.snp_df = self.snp_df[self.snp_df["GT"].apply(
                filter_alt1_only)]

        if not self.filter_af_none_include:
            self.snp_df = self.snp_df[self.snp_df["AF"].notna()]

        if self.filter_af_lower:
            self.snp_df = self.snp_df[self.snp_df["AF"] >=
                                      self.filter_af_lower]

        if self.filter_af_upper:
            self.snp_df = self.snp_df[self.snp_df["AF"] <=
                                      self.filter_af_upper]

        # alt:sample:allele solver
        self.snp_df["ALT"] = self.snp_df.apply(
            solve_alt, axis=1, allele=self.allele)

        # по количеству
        try:
            self.snp_df = self.snp_df.sample(self.consensus_num)
        except:
            pass

        self.snp_df = self.snp_df[["CHROM", "POS",
                                   "REF", "ALT", "GT"]]
        self.snp_df["POS"] = self.snp_df["POS"].astype(int)

    def split_fasta(self) -> None:
        print("\nSNPs to use")
        print(self.snp_df)
        print('\n')
        progress = 0
        snp_df_len = len(self.snp_df)

        cur_chr = ""
        chr_seq = ""
        start_flag = True
        required_chr = True

        cur_snp_df = pd.DataFrame()

        with open(self.fasta_path, "r") as fasta:
            while True:
                line = fasta.readline()

                if line and not line.startswith(">"):
                    if required_chr:
                        chr_seq = chr_seq + line.strip()
                    continue

                if not start_flag:
                    print(str(cur_chr)+" chr in work")

                    for _, snp in cur_snp_df.iterrows():
                        mutated_seq = mutate_seq(snp=snp, chr_seq=chr_seq,
                                                 consensus_len=self.consensus_len,
                                                 mutation_place=self.mutation_place)
                        save_consensus(seq=mutated_seq,
                                       snp=snp, output_dir=self.output_dir)

                        progress += 1
                        print(str(progress)+"/"+str(snp_df_len) +
                              " SNPs are processed")

                if progress == snp_df_len or not line:
                    print("\nDone. Results in dir: " + str(self.output_dir))
                    break

                cur_chr = line[1:].split(" ")[0]
                cur_snp_df = self.snp_df[self.snp_df["CHROM"] == cur_chr]
                required_chr = len(cur_snp_df) != 0
                start_flag = False
                chr_seq = ""


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=DESC)

    parser.add_argument("-vcf", "--vcf_path", required=True, help="vcs file")
    parser.add_argument("-fasta", "--fasta_path", required=True,
                        help="reference genome fasta file")

    parser.add_argument("-len", "--consensus_len", default=70, type=int,
                        help="length of consensuses", required=True)
    parser.add_argument("-num", "--consensus_num", default=10, type=int,
                        help="number of consensuses", required=True)

    parser.add_argument("-out", "--output_dir",
                        default="output", help="output directory")

    parser.add_argument("-a", "--allele", default=1, type=int,
                        help="allele index (only 1 is tested)")
    parser.add_argument("-s", "--sample", default=None,
                        help="sample name to use")
    parser.add_argument("--mutation_place", default="center", choices=["start", "center", "end"],
                        help="place of mutation")

    # filtering
    parser.add_argument("--filter_af_lower", default=None, type=float,
                        help="lower value of AF")
    parser.add_argument("--filter_af_upper", default=None, type=float,
                        help="upper value of AF")
    parser.add_argument("--filter_af_none_include", default=False, type=bool,
                        help="include or not snp with AF=None")
    parser.add_argument("--filter_alt1_only", default=True, type=bool,
                        help="if few ALT, use only first one (only True is tested)")
    parser.add_argument("--filter_homo", default=True, type=bool,
                        help="use only homologous alleles (only True tested)")
    parser.add_argument("--filter_pass", default=True, type=bool,
                        help="only FILTER=PASS")

    parser.add_argument("--other_args", nargs="*")

    args = vars(parser.parse_args())

    creator = ConsensusSeqCreator(**args)
