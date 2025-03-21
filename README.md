# scripts
This is where we store all software tools that do not require installation, compilation or doing anything else before you run them.

## **Important info for beginners:**

* **.py** files are Python scripts, i.e. they require Python to run.  

* **.pl** files are Perl scripts, i.e. they require Perl to run.  

* **.sh** files are bash scripts, i.e. they require nothing to run (since bash is everywhere by default).  

* **.R** files are R scripts, i.e. they require R to run.


## **Currently in the repository:**

* **[branchlen.py](https://github.com/ProtistomicsLab/scripts/blob/main/branchlen.py)** (by Kacper - for extracting branch length data from Newick trees)


* **[fasta_to_phylip.py](https://github.com/ProtistomicsLab/scripts/blob/main/fasta_to_phylip.py)** (by Jadzia - for converting FASTA sequences to Phylip format, used e.g. in PAML)

  usage: fasta_to_phylip.py <input.fa> <output.phy>

* **[getting_completness_of_the_genome_based_on_transcriptome.py](https://github.com/ProtistomicsLab/scripts/blob/main/getting_completness_of_the_genome_based_on_transcriptome.py)** (by Jadzia - for estimating genome completness based on the available transcriptome from the same organism)

  usage: getting_completness_based_on_transcriptome.py [-h] -p <predicted_peptides.pep.fa> -db_f <fasta_to_be_used_as_blast_database.fasta> -t <NUMBER_OF_THREADS> -l <annotations_list.txt> -cds <predicted_CDS.fasta> -prots <predicted_proteins.fasta>

* **[newick_to_nexus.py](https://github.com/ProtistomicsLab/scripts/blob/main/newick_to_nexus.py)** (by Jadzia)
A python script that converts newick to nexus.

usage: newick_to_nexus.py [input] [output]

[input] name of the input file in fasta newick [output] name of the output file in nexus format

* **[stockholm2fasta.py](https://github.com/ProtistomicsLab/scripts/blob/main/stockholm2fasta.py)** (by Jadzia)
A Python script to convert stockholm format to fasta format.

usage: stockholm2fasta.py [input]

Note! The output will be automaticly saved as [input].aligned.fasta
