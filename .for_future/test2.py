import csv
from sys import argv

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

script = ""
ref = ""
subs = ""
outfile = "new.txt"
out = open(outfile, 'wb')
seq_string = []
# user_descrip=raw_input("Give the file a description: ")
user_descrip = "Give the file a description: "

with open(ref, "r") as infile:
    for k, line in enumerate(infile):
        if ">" in line:
            pass
        else:
            seq_string.append(line.strip("\n"))

seq_str = "".join(seq_string)
print(len(seq_str))
new_indexed_seq = list(seq_str)
with open(subs) as subs_file:
    subs_reader = csv.reader(subs_file, delimiter='\t')
    for line in subs_reader:
        pos = (int(line[0]) - 1)
        base = line[1]
        original_base = seq_str[pos]

        new_indexed_seq[pos] = base

out_string = "".join(new_indexed_seq)

with open(outfile, 'wb') as out:
    out_seq = SeqRecord(
        Seq(out_string), id='anscestral_sequence', description=user_descrip)
    SeqIO.write(out_seq, out, 'fasta')
