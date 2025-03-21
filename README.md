# scripts
This is where we store all software tools that do not require installation, compilation or doing anything else before you run them.

## **Important info for beginners:**

* **.py** files are Python scripts, i.e. they require Python to run.  

* **.pl** files are Perl scripts, i.e. they require Perl to run.  

* **.sh** files are bash scripts, i.e. they require nothing to run (since bash is everywhere by default).  

* **.R** files are R scripts, i.e. they require R to run.


## **Currently in the repository:**

* **[branchlen.py](https://github.com/ProtistomicsLab/scripts/blob/main/branchlen.py)** (by Kacper - for extracting branch length data from Newick trees)


* **[fasta_to_phylip.py](https://github.com/ProtistomicsLab/scripts/blob/main/fasta_to_phylip.py)** (by Jadzia)
Python script that converts input in fasta to format into output in phylip format. Usesful while running PAML or CODEML.

usage: fasta_to_phylip.py [input] [output]

[input] name of the input file in fasta format [output] name of the output file in phylip format

* **[getting_completness_of_the_genome_based_on_transcriptome.py](https://github.com/ProtistomicsLab/scripts/blob/main/getting_completness_of_the_genome_based_on_transcriptome.py)** (by Jadzia)
Python script to estimate genome completness based on transcriptome.

usage: getting_completness_based_on_transcriptome.py [-h] -p PREDICTED_PEPTIDES -db_f FASTA_FILE_FOR_DATABASE_FOR_BLASTING -t NUMBER_OF_THREADS -l LIST_OF_ANNOTATED -cds PREDICTED_CODING_SEQ -prots PREDICTED_PROTS

* **[newick_to_nexus.py](https://github.com/ProtistomicsLab/scripts/blob/main/newick_to_nexus.py)** by Jadzia
* A python script that converts newick to nexus.

* usage: newick_to_nexus.py [input] [output]

* [input] name of the input file in fasta newick [output] name of the output file in nexus format
* 
