'''
This code was written in Python 3.7

Purpose of this code: 
-Given a .txt file containing FASTA-formatted DNA sequences, calculate the codon frequency at each codon position and output it in .csv format for easy handling by Excel/Graphpad

EXAMPLE INPUT FORMAT:
>Example_sequence_1
ATGTTTTAT
>Example_sequence_2
ATGTATTTT


EXAMPLE OUTPUT FORMAT:
Codon #, Example_sequence_1, Example_sequence_2
1, 1.00, 1.00
2, 0.26, 0.73
3, 0.73, 0.26

'''

import Bio
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO
import os.path
import hive_cut_database_manager
import re
import numpy as np
import pandas as pd

#Uses class methods from hive_cut_database_manager.py
#Asks what organism to use for the codon optimization table (either ecoli or human)
while True:
	_Optimization_Organism = input("Default organism codon table to use if organism origin not specified (ecoli/human) ? : ")
	if _Optimization_Organism.lower() == 'ecoli':
		_Optimization_Organism = 'Escherichia_coli'
		_Optimization_Organism_Table = hive_cut_database_manager.species('Escherichia_coli').get_HIVE_CUT_frequency_python_dictionary()
		break
	if _Optimization_Organism.lower() == 'human':
		_Optimization_Organism = 'Homo_sapiens'
		_Optimization_Organism_Table = hive_cut_database_manager.species('Homo_sapiens').get_HIVE_CUT_frequency_python_dictionary()
		break
	else:
		continue

class DNA_sequence:
	'''Instances represent DNA sequences (only ATCG bases) that can optionally be named by the seq_name argument, and organism argument specifies LATIN name of origin organism'''
	def __init__(self, DNA_seq, seq_name=None, organism=_Optimization_Organism):
		try:
			self.DNA_seq = Seq(DNA_seq.upper())
		except TypeError:
			self.DNA_seq = DNA_seq.upper()

		self.protein_seq = self.DNA_seq.translate()

		self.seq_name = seq_name
		self.organism = organism
		#Sets codon frequency table to user chosen default organism Kazusa table
		self.codon_frequency = _Optimization_Organism_Table

	def generate_codon_frequency_dataframe(self):
		'''Function to generate a list of a DNA sequence's codon frequencies Ex. ATGTCA = [1.00, 0.59] '''
		codon_frequency_list = []
		codon_position = 0
		codon_position_list = []
		for position in range(0,len(self.DNA_seq),3):
			codon = self.DNA_seq[position:position+3]
			codon_position = codon_position + 1
			codon_position_list.append(codon_position)
			codon_frequency_list.append(self.codon_frequency[codon])
		codon_frequency_df = pd.DataFrame(data=codon_frequency_list, index=codon_position_list,columns=[self.seq_name])
		return codon_frequency_df


if __name__ == "__main__":
	while True:
	# #Assigns variable from user input of the FASTA-formated .txt file containing DNA sequences to be displayed
		filename = str(input("Fasta-formatted filename containing DNA sequences? (exit to exit): "))
		if filename == 'exit':
			break
		if os.path.isfile(filename):
			break
		else:
			print("Invalid file. Please enter valid file name.")

	if filename != "exit":
		#Parses file containing FASTA formatted DNA sequences
		file_data = list(Bio.SeqIO.parse(filename,"fasta"))

		#Creates list containing the FASTA names of each DNA sequence
		seq_list = [record.id for record in file_data]
		codon_freq_df_records = pd.DataFrame()
		print('Using %s codon frequency table.' % _Optimization_Organism)

		#Iterates through each FASTA record in the file and appends its codon frequency at each codon position to a master dataframe that is exported as a .csv for plotting in other applications
		for record in file_data:
			record = DNA_sequence(str(record.seq), seq_name=record.id, organism=_Optimization_Organism)
			record_codon_freq_df = record.generate_codon_frequency_dataframe()
			codon_freq_df_records= codon_freq_df_records.join(record_codon_freq_df, how='outer')
		output_filename = input('Output filename .csv : ? ')
		codon_freq_df_records.to_csv(output_filename)
		print('Wrote .csv to %s' % output_filename)
	else:
		print("Exiting program")



