'''
This code was written in Python 3.7

Purpose of this code: 
-Given a user-entered .txt file containing FASTA-formatted DNA sequences of identical length, 
calculate the pairwise percent codon identity of each sequence (Percent ID being how many codons are identical at a given codon position)

EXAMPLE INPUT FORMAT (file.txt):
>1
ATGAGTGCCTGGCACATCTAC
>2
ATGTCGGCATGGCATATCTAC
>3
ATGTCCGCGTGGCATATTTAC
>4
ATGAGCGCGTGGCACATATAC

EXAMPLE OUTPUT FORMAT:
file.txt (input)
Sequence 1, Sequence 2, % Codon ID 
1,2,57.14
1,3,42.86
1,4,57.14
2,3,57.14
2,4,42.86
3,4,57.14

'''

import itertools

class DNA_sequence:
	'''Instances should be NAMED DNA sequences. 
	Ex. DNA_sequence(seq_name=native_DNA, DNA_seq=ATGTATCAGACT)
	'''
	def __init__(self, seq_name, DNA_seq):
		self.seq_name = seq_name
		self.DNA_seq = DNA_seq

	def percent_codon_ID(self,comparison_seq):
		'''Given a comparison sequence, calculates codon identity'''
		assert (len(self.DNA_seq) == len(comparison_seq)), "Sequences are not identical length!"

		identity = 0
		number_of_codons = len(self.DNA_seq)/3

		#Coutnts how many codon positions have identical codons between the two sequences
		for bp_position in range(0, len(self.DNA_seq),3):
			seq_1_codon = self.DNA_seq[bp_position:bp_position+3]
			seq_2_codon = comparison_seq[bp_position:bp_position+3]
			if seq_1_codon == seq_2_codon:
				identity += 1
			else:
				continue
		percent_ID = round(identity/number_of_codons * 100,2)
		return percent_ID

if __name__ =='__main__':
	#Ask user for input and output filenames
	while True:
		filename = input('Input FASTA-formatted filename (.txt)? : ')
		if filename.endswith('.txt'):
			break
		else:
			print('Invalid filename. Input filenames must end with .txt')
	while True:
		output_filename = input('Output filename (.csv)? : ')
		if output_filename.endswith('.csv'):
			break
		else:
			print('Invalid filename. Output filenames must end with .csv')

	#Generates dictionary of {Seq_name:sequence} from file
	with open(filename,'r') as file:
		names = []
		sequences = []
		sequence_holder = ''
		for line in file:
			line = line.strip()
			if line.startswith('>'):
				names.append(line[1:])
				sequences.append(sequence_holder)
				sequence_holder = ''
			else:
				sequence_holder += line
		sequences.append(sequence_holder)
		del sequences[0]
		seq_dict = dict(zip(names,sequences))

	#Generates list of DNA sequences instances from the file
	list_of_seqs = []
	for name, sequence in seq_dict.items():
		list_of_seqs.append(DNA_sequence(seq_name=name, DNA_seq=sequence))

	print(filename)

	#Prints out all possible unique combinations of sequences
	for seq1, seq2 in itertools.combinations(list_of_seqs,2):
		print("Percent ID of %s and %s : %.2f" % (seq1.seq_name, seq2.seq_name, seq1.percent_codon_ID(seq2.DNA_seq)))
	print('\n')

	#Write codon percent ID to file
	print("Writing to .csv output file: %s" % output_filename)
	with open(output_filename,'a') as output_file:
		output_file.write(filename + '\n')
		output_file.write("Sequence 1, Sequence 2, % Codon ID \n")
		for seq1, seq2 in itertools.combinations(list_of_seqs,2):
			output_file.write("%s,%s,%.2f" % (seq1.seq_name, seq2.seq_name, seq1.percent_codon_ID(seq2.DNA_seq)))
			output_file.write("\n")
		output_file.write('\n')