import sys
from Bio import SeqIO

in_fasta = SeqIO.parse(open(sys.argv[1]), 'fasta')
out_gtf = sys.argv[1].rsplit('.', 1)[0] + ".gtf"

with open(out_gtf, "w") as gtf_out:
    for record in in_fasta:
        seq_id = record.id
        length = len(record.seq)
        gene_id = seq_id.split("_i")[0]
        line = (
            f"{seq_id}\ttranscriptome\tCDS\t1\t{length}\t.\t+\t.\t"
            f'gene_id "{gene_id}"; transcript_id "{seq_id}"; exon_number 1\n'
            f"{seq_id}\ttranscriptome\texon\t1\t{length}\t.\t+\t.\t"
            f'gene_id "{gene_id}"; transcript_id "{seq_id}"; exon_number 1\n'
        )
        gtf_out.write(line)