'''
This code was written in Python 3.7

Purpose of this code: 
-Given a user-entered protein sequence, generate a user-given number of sequences using random reverse translation and write them into a .txt file in FASTA format

EXAMPLE INPUT FORMAT:
MSAWHIY

EXAMPLE OUTPUT FORMAT:
>1
ATGAGTGCCTGGCACATCTAC
>2
ATGTCGGCATGGCATATCTAC
>3
ATGTCCGCGTGGCATATTTAC
>4
ATGAGCGCGTGGCACATATAC
'''

#Dictionary of amino acid:synonymous codons
synonymous_codon_dict = {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'], 
          'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
          'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'], 
          'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'], 
          'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'], 
          'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'], 
          'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}

import random

def random_reverse_translation(protein_seq):
	"""Given a protein sequence, generates a random reverse sequence"""
	protein_seq = protein_seq.upper()
	random_reverse_seq = ''
	for residue in protein_seq:
		synonymous_codons = synonymous_codon_dict[residue]
		random_number = random.randint(0,len(synonymous_codons)-1)
		random_reverse_seq += synonymous_codons[random_number]
	return random_reverse_seq

if __name__ == '__main__':
	input_seq = input('Protein sequence to generate random reverse translated DNA sequence: ')
	number_of_random_seqs = int(input('Number of random seqs to generate? : '))
	output_filename = input('Output filename? : ')
	
	#Generates random protein sequences and writes to file in .fasta format
	for i in range(0,number_of_random_seqs):
		name = '>' + str(i+1)
		with open(output_filename,'a') as f:
			f.write(name)
			f.write('\n')
			f.write(random_reverse_translation(input_seq))
			f.write('\n')