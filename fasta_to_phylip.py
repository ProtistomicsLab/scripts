# A script to convert fasta to phylip.
# Tested and working on Python 3.6 

from Bio import AlignIO
from pathlib import Path
from Bio.AlignIO.PhylipIO import SequentialPhylipWriter
import re
import sys

file = open(sys.argv[1])
file_name=sys.argv[1]
out_name = file_name.split(".")[0]
out_name = out_name + ".phy"


print(out_name)


#fileoutnex = re.sub('fasta', 'nex', file)
#pathway = '/mypath/'
alignment = AlignIO.read(file, 'fasta')
#print(alignment)
# write out in nexus
#AlignIO.write(alignment, Path(pathway, fileoutnex), 'nexus')

# write out in relaxed phylip
#fileoutphy = re.sub('fasta', 'phy', file)
with open(out_name, 'w') as output_handle:        
    SequentialPhylipWriter(output_handle).write_alignment(alignment, id_width=30)
