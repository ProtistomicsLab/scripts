#!/usr/bin/python3.6


import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
import csv
import collections
import subprocess

def bash_command(cmd):

	command_run = subprocess.call(cmd, shell=True, executable='/bin/bash')

	if command_run == 0:
		return True
	else:
		return False


def getting_complete_proteins(prots):


	records = SeqIO.parse(prots, 'fasta')
	q_dict = SeqIO.index(prots, 'fasta')

	full_prots = []

	for record in records:

		if (str(record.seq)).startswith("M") and (str(record.seq)).endswith("*"):

			full_prots.append(record.id)
	to_keep = [q_dict[name] for name in full_prots]

	return SeqIO.write(to_keep, "full_proteins.fasta", "fasta")

def blasting_against_db(db_f, thre):

	bash_command(f'makeblastdb -in {db_f} -dbtype prot -out {db_f}_db')

	return bash_command(f'blastp -db {db_f}_db -query full_proteins.fasta -out full_proteins.blastp -num_threads {thre} -outfmt 6')


def getting_top_hit_based_on_evalue(l):

	with open("full_proteins.blastp", "r") as blast, open(l, 'r') as annotated:

		data = pd.read_csv(blast, sep = "\t", names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
		df = pd.DataFrame(data)


		anno = annotated.readlines()
		anno_into_list = []

		for i in anno:

			i = i.strip()
			anno_into_list.append(i)

		df_50 =	df[(df['evalue'] <= 1e-50)]
		df_100 = df[(df['evalue'] <= 1e-100)]


		result_df_50 = df_50[df_50['sseqid'].isin(anno_into_list)]
		result_df_100 = df_100[df_100['sseqid'].isin(anno_into_list)]

		df_50_l = set(result_df_50["qseqid"].values.tolist())
		df_100_l = set(result_df_100["qseqid"].values.tolist())


		print(f'Threshold 1e-50 gave a DataFrame with {df_50.shape[0]} rows out of which {len(df_50_l)} have annotations.')
		print(f'Threshold 1e-100 gave a DataFrame with {df_100.shape[0]} rows out of which {len(df_100_l)} have annotations')

		q_dict = SeqIO.index("full_proteins.fasta", 'fasta')
		to_keep = [q_dict[name] for name in df_100_l]

		return SeqIO.write(to_keep, "sleceted_prots_with_annotations.aa", "fasta"), len(df_100_l)


def blast_out_parser(blastout):

	data = pd.read_csv(blastout, sep = "\t", names = ['qseqid', 'sseqid', 'pident', 'qlen', 'length', 'evalue', 'bitscore'])
	df = pd.DataFrame(data)

	df['pident_treshold'] = df['pident'] >= 90
	df = df[df.pident_treshold]

	df['len_perc'] = df['length'] / df['qlen'] * 100
	df['len_perc'] = df['pident'] >= 50
	df = df[df.len_perc]

	df_l = set(df["qseqid"].values.tolist())

	return df_l


def calculating_completness(cds, pep, numb):


	bash_command(f'makeblastdb -in {cds} -dbtype nucl -out ./sleceted_nucl_with_annotations_db')
	bash_command(f'makeblastdb -in {pep} -dbtype prot -out ./sleceted_pep_with_annotations_db') 

	bash_command(f'blastp -db sleceted_pep_with_annotations_db -query sleceted_prots_with_annotations.aa -out prots_against_slected.blastp -num_threads {thre} -outfmt "6 qseqid sseqid pident qlen length evalue bitscore"')
	bash_command(f'tblastn -db sleceted_nucl_with_annotations_db -query sleceted_prots_with_annotations.aa -out cds_against_selected.tblastn -num_threads {thre} -outfmt "6 qseqid sseqid pident qlen length evalue bitscore"')



	with open("prots_against_slected.blastp", 'r') as prots, open("cds_against_selected.tblastn", 'r') as cds:

		prots_l = len(blast_out_parser(prots))
		cds_l = len(blast_out_parser(cds))

		numb = numb[1]

		bla = (prots_l * 100) / numb
		bla2 = (cds_l * 100) / numb 

		print(f'Completness based on predicted proteins ---> {bla} %')
		print(f'Completness based on predicted cds ---> {bla2}%')


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Calculates genome completness based on the proteins predicted from transcripts')

	required = parser.add_argument_group('required arguments')

	required.add_argument('-p', '--predicted_peptides', type = str, help = 'Proteins predicted from the transcritpome assembly', required = True)
	required.add_argument('-db_f', '--fasta_file_for_database_for_blasting', type = str, help = 'Full path to the database for blasting', required = True)
	required.add_argument('-t', '--number_of_threads', type = str, help = 'Number of threads for blasting', required = True)
	required.add_argument('-l', '--list_of_annotated', type = str, help = 'Number of threads for blasting', required = True)
	required.add_argument('-cds', '--predicted_coding_seq', type = str, help = 'Coding sequences predicted with AUGUSTUS', required = True)
	required.add_argument('-prots', '--predicted_prots', type = str, help = 'Proteins predicted with AUGUSTUS', required = True)

	args=parser.parse_args()

	prots = args.predicted_peptides
	thre = args.number_of_threads
	l = args.list_of_annotated
	db_f = args.fasta_file_for_database_for_blasting
	cds = args.predicted_coding_seq
	pep = args.predicted_prots

	getting_complete_proteins(prots)
	blasting_against_db(db_f, thre)
	numb = getting_top_hit_based_on_evalue(l)

	calculating_completness(cds, pep, numb)
