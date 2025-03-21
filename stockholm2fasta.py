from Bio import SeqIO
import sys


records = SeqIO.parse(sys.argv[1], "stockholm")

input = sys.argv[1]
output = input.split(".")[0]
out = output + ".aligned.fasta"

count = SeqIO.write(records, out, "fasta")
print("Converted %i records" % count)
