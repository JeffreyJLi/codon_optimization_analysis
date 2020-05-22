'''
This code was written in Python 3.7
Purpose of this code: 
-Generate and manage a .txt file containing HIVE-CUT database codon usage tables for different species
-This code is not utilized by itself; it is imported and used in other scripts (ex. violin_plot_csv_generator.py)

Overview:
The "species" class has access to three methods:
1. get_HIVE_CUT_codon_frequency_table()
-This method first opens the "HIVE_CUT_codon_frequency_database.txt" and checks if the codon usage table for the given species is present
-If the table is not found, it executes the import_HIVE_CUT_codon_frequency_table() method to allow the user to import the codon usage table

2. import_HIVE_CUT_codon_frequency_table()
-This method reads a .txt file containing a single HIVE-CUT codon usage table (must be in the format shown below)
TTT 17.08 (1341983) TCT 16.88 (1326337) TAT 12.06 ( 947642) TGT 10.48 ( 823726)
TTC 17.53 (1377246) TCC 17.35 (1363614) TAC 13.46 (1057591) TGC 10.90 ( 856646)
TTA  8.78 ( 690119) TCA 14.15 (1111829) TAA  0.43 (  33559) TGA  0.77 (  60621)
TTG 13.42 (1054214) TCG  4.08 ( 320758) TAG  0.34 (  26407) TGG 11.67 ( 916734)
                       
CTT 14.11 (1108931) CCT 19.07 (1498271) CAT 11.87 ( 932958) CGT  4.56 ( 358155)
CTC 17.85 (1402940) CCC 18.98 (1491270) CAC 14.65 (1151545) CGC  8.78 ( 690020)
CTA  7.49 ( 588397) CCA 18.67 (1466930) CAA 14.15 (1111835) CGA  6.39 ( 502451)
CTG 36.24 (2848162) CCG  6.19 ( 486595) CAG 35.42 (2783737) CGG 10.70 ( 841125)
                       
ATT 16.50 (1296263) ACT 14.26 (1120650) AAT 18.40 (1445838) AGT 13.98 (1098334)
ATC 18.67 (1467275) ACC 17.80 (1398865) AAC 18.22 (1431977) AGC 19.70 (1547892)
ATA  8.14 ( 639457) ACA 16.54 (1299919) AAA 27.65 (2172879) AGA 13.32 (1046695)
ATG 21.39 (1680759) ACG  5.64 ( 443078) AAG 31.85 (2503111) AGG 12.19 ( 958126)
                       
GTT 11.72 ( 920619) GCT 18.85 (1481387) GAT 23.99 (1885552) GGT 10.77 ( 846406)
GTC 13.44 (1055851) GCC 25.84 (2030499) GAC 24.29 (1908799) GGC 19.83 (1558404)
GTA  7.67 ( 603076) GCA 17.05 (1339676) GAA 33.72 (2649394) GGA 17.07 (1341033)
GTG 25.90 (2035509) GCG  6.03 ( 473616) GAG 39.75 (3123257) GGG 15.33 (1204755)

3. get_HIVE_CUT_frequency_python_dictionary()
-Parses the species's HIVE-CUT codon usage table into a Python dictionary
'''

import re
import os
codon_translation ={ 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    }

#Dictionary containing amino acids and their synonymous codons
amino_acid_codon_dict = {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'], 
          'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
          'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'], 
          'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'], 
          'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'], 
          'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'], 
          'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}

class species:
  '''Given a species name, pulls codon frequency table from the pre-existing database file, or if not present, asks the user to input a .txt file containing the codon frequency table for the given organism'''
  def __init__(self, species_name):
    '''Records input record name from Biopython-parsed .fasta file record.id as the name of the species name used in the codon frequency database'''
    self.name = species_name.capitalize().replace(' ','_')
    self.HIVE_CUT_codon_frequency_table = self.get_HIVE_CUT_codon_frequency_table()

  def get_HIVE_CUT_codon_frequency_table(self):
    ''' Function to check if the codon frequency table from HIVE-CUT of the given species is present in the local database file HIVE_CUT_codon_frequency_database.txt'''
    species_database_name = '>' + self.name #Codon tables in the database are labeled in .fasta format with >example_protein_name (ex. >Homo_sapiens)
    species_record_found = False
    species_database_table = '' #Holder variable for HIVE_cut codon frequency table

    try:
      #This loop opens HIVE_CUT_codon_frequency_database.txt file and if the input species name is present, pulls out its associated database
      with open('HIVE_CUT_codon_frequency_database.txt','r') as HIVE_cut_database:
        for line in HIVE_cut_database:
          line = line.strip()
          if line == species_database_name:
            species_record_found = True
            continue
          if species_record_found:
            if line.startswith('>'):
              break
            else:
              species_database_table+=line.rstrip('\n') + '\n'
      
      #Returns the input species codon frequency table, if found. Else, asks if user would like to import one and then recalls the method to return it. If not, returns None            
      if species_database_table != '':
        return species_database_table
      else:
        print("HIVE CUT codon frequency table not found for species: %s." % self.name)
        while True:
          import_decision = input('Would you like to import the table from another .txt file? (y/n): ')
          if import_decision.lower() == 'y':
            self.import_HIVE_CUT_codon_frequency_table()
            return self.get_HIVE_CUT_codon_frequency_table()
          elif import_decision.lower() == 'n':
            return None

    except FileNotFoundError:
      print("HIVE_CUT_codon_frequency_database.txt file not found.")
      while True:
        import_decision = input('Would you like to import the table from another .txt file? (y/n): ')
        if import_decision.lower() == 'y':
          self.import_HIVE_CUT_codon_frequency_table()
          return self.get_HIVE_CUT_codon_frequency_table()
        elif import_decision.lower() == 'n':
          return None
      return None

  def import_HIVE_CUT_codon_frequency_table(self):
    '''Appends a .txt file containing a given codon frequency table from HIVE CUT'''
    while True:
      input_filename = input("name of .txt file containing a singular HIVE cut codon frequency table corresponding to %s? (exit to exit): " % self.name)
      if input_filename == 'exit':
        return "Exiting table import"
      if os.path.isfile(input_filename):
        break
      else:
        print("Invalid file. Please enter valid file name.")

    '''Appends input HIVE CUT codon frequency table to the HIVE_CUT_codon_frequency_database.txt database file'''
    with open('HIVE_CUT_codon_frequency_database.txt', 'a') as HIVE_cut_database:
      with open(input_filename, 'r') as input_file:
        
        firstline = input_file.readline().strip()
        assert (firstline.startswith('TTT')), "File in incorrect format"

        HIVE_cut_database.write('\n')
        HIVE_cut_database.write(">" + self.name + "\n")
        HIVE_cut_database.write(firstline + "\n")
        for line in input_file:
          HIVE_cut_database.write(line.strip() + "\n")
    print("Appended %s HIVE CUT frequency table to HIVE_CUT_codon_frequency_database.txt" % self.name)
    return None

  def get_HIVE_CUT_frequency_python_dictionary(self):
    #Will only run this code if there is a HIVE CUT codon frequency table for the input species (or else, it doesn't make sense)
    assert self.HIVE_CUT_codon_frequency_table != None, "No HIVE CUT codon frequency table for species %s" % self.name

    residue_frequencyperthousandsum_dict = {}
    codon_fractional_frequency = {}

    #From the HIVE CUT codon frequency table, uses a regex to parse out each codon and their associated codon frequency per thousand into a separate holder dictionary
    for codon, frequency_per_thousand in re.findall(r"([ATCG]{3})\s+(\d+.\d*)\s+\(",self.HIVE_CUT_codon_frequency_table):
      codon_residue = codon_translation[codon]
      residue_frequencyperthousandsum_dict[codon_residue] = residue_frequencyperthousandsum_dict.get(codon_residue,0) + float(frequency_per_thousand)

    #Sums up each residue's total frequency per thousand, and then uses that to calculate each codon's fractional frequency (fractional_frequency = codon's frequency per thousand/residue's total frequency per thousand)
    for codon, frequency_per_thousand in re.findall(r"([ATCG]{3})\s+(\d+.\d*)\s+\(",self.HIVE_CUT_codon_frequency_table):
      codon_residue = codon_translation[codon]
      fractional_frequency = float(frequency_per_thousand)/residue_frequencyperthousandsum_dict[codon_residue]
      codon_fractional_frequency[codon] = fractional_frequency

    #Final check to see if all codons have been successfully calculated. If not, there has been an error.
    assert (len(codon_fractional_frequency.keys()) == 64), "Table does not include 64 codons!"

    return codon_fractional_frequency

