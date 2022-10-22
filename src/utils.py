import math
import os
import shutil
import pandas as pd

from consts import (
    FASTA_LINE_LEN
)


def enrich_gt(row: pd.Series, sample: str):
    try:
        GT_index = row["FORMAT"].split(":").index("GT")
        return row[sample].split(":")[GT_index]
    except:
        pass


def enrich_af(row: pd.Series):
    try:
        AF_field = [x for x in row["INFO"].split(";") if "AF" in x][0]
        return float(AF_field.split("=")[1])
    except:
        pass


def filter_homo(value: str):
    if "/" in value:
        return len(set(value.split("/"))) == 1
    if "|" in value:
        return len(set(value.split("|"))) == 1
    return False


def filter_alt1_only(value: str):
    try:
        if "/" in value:
            return int(value.split("/")[0]) == 1
        if "|" in value:
            return int(value.split("|")[0]) == 1
        return False
    except:
        return False


def solve_alt(row: pd.Series, allele: int):
    try:
        if "/" in row["GT"]:
            alt_index = int(row["GT"].split("/")[allele-1])
        elif "|" in row["GT"]:
            alt_index = int(row["GT"].split("|")[allele-1])
        else:
            return row["ALT"]

        if alt_index == 0:
            return row["REF"]

        return row["ALT"].split(",")[alt_index-1]

    except Exception as ex:
        print("error", ex)
        return row["ALT"]


def mutate_seq(snp: pd.Series, consensus_len: int, chr_seq: str, mutation_place: str):
    chr_seq_len = len(chr_seq)

    alt_len = len(snp["ALT"])
    ref_len = len(snp["REF"])
    diff_len = ref_len - alt_len
    ref_start = snp["POS"] - 1
    ref_end = ref_start + ref_len
    default_intervals_len = (consensus_len - alt_len)/2

    consensus_start = 0
    consensus_end = chr_seq_len - 1

    # граничные условия
    match mutation_place:
        case "center":
            default_start = ref_start - math.floor(default_intervals_len)
            default_end = ref_end + math.ceil(default_intervals_len)
            shifted_start = consensus_end - consensus_len - diff_len
            shifted_end = consensus_start + consensus_len + diff_len

            if default_start < 0 and not shifted_end >= chr_seq_len:
                consensus_end = shifted_end
            if default_end >= chr_seq_len and not shifted_start < 0:
                consensus_start = shifted_start
            if not default_start < 0 and not default_end >= chr_seq_len:
                consensus_start = default_start
                consensus_end = default_end
        case "start":
            default_end = ref_start + consensus_len + diff_len
            consensus_start = ref_start
            if not default_end >= chr_seq_len:
                consensus_end = default_end
        case "end":
            default_start = ref_end - consensus_len - diff_len
            consensus_end = ref_end
            if not default_start < 0:
                consensus_start = default_start

    first_part = chr_seq[consensus_start:ref_start]
    second_part = chr_seq[ref_end:consensus_end]
    mutated_seq = first_part + snp["ALT"] + second_part
    original_seq = chr_seq[consensus_start:consensus_end]
    # print(consensus_len, len(mutated_seq))
    # print(original_seq)
    # print(mutated_seq)

    return mutated_seq


def prepare_dir(output_dir: str):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)


def save_consensus(seq: str, snp: pd.Series, output_dir: str):
    name = " ".join(snp[["CHROM", "POS", "REF", "ALT"]].astype(
        str).values.tolist())

    fasta_path = os.path.join(output_dir, name+".fasta")

    wrapped_seq = [seq[y-FASTA_LINE_LEN:y]
                   for y in range(FASTA_LINE_LEN, len(seq)+FASTA_LINE_LEN, FASTA_LINE_LEN)]

    with open(fasta_path, "w", encoding="utf-8") as fasta_file:
        fasta_file.write(">"+name+"\n")
        for line in wrapped_seq:
            fasta_file.write(line+"\n")


# df = pd.Series({'REF': 'BB', 'POS': 3, 'ALT': 'A'})
# print(df)

# seq = 'AABFDFSADADVFVRASSC'
# print(len(seq))
# print(mutate_seq(snp=df, mutation_place='center',
#       chr_seq=seq, consensus_len=5))
