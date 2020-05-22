# Project title

These scripts are associated with the Ranaghan et al. manuscript analyzing codon optimization algorithms. These files are broken into two sections, pairwise_alignment and violin_plot_csv_generator. The pairwise alignment directory contains a script used to generate random reverse translations of a input protein sequence, and a script used to calculate the pairwise codon identity of equal-length DNA sequences. The violin plot csv generator directory contains a script to create/manage/parse a .txt file containing HIVE-CUT codon usage tables. The other script takes in a .txt file containing FASTA-formatted DNA sequences and calculates the codon frequency at each codon position of each sequence. This data is outputted into .csv format for easy access by Excel/Graphpad for plotting as a violin plot. 

## Prerequisites

Python files were written in Python 3.7

numpy
pandas
Biopython

## Quickstart

Installation:
'''
pip install -r requirements.txt
'''

Usage:
'''
python path/script_to_use.py
'''

See individual script files for additional details on how to use.