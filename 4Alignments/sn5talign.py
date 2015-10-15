#!/usr/bin/env python

# use Mafft alignments (from the cutalign.py script) and run transAlign.pl to align the amino acids. 
# cutalign.py used default -sn threshold, here change to -sn5

# ----------------
# Import libraries
# ----------------
import os, sys
from collections import defaultdict

# -------------------
# get the input files
# -------------------
thefiles = [f for f in os.listdir('./mafft_alignments/')]
#print "thefiles: ", thefiles

for myfile in thefiles:
#	print "myfiles: ", myfile

	# this takes all files in the folder and only selects one with .aligned
	if '.aligned' in myfile:
		# this takes specified files from if statement and removes the extension so you just have the file name
		mygene = myfile.split('_')[0]
#		print "mygene: ", mygene
		
		# make tAlign_fastaouts/ directory first. tranAlign alignment
		# align with transAlign.pl. Allow 5 stop codons. alpha=5% in one-tailed t-test.
		# transAlign.pl -h for help
		os.system('transAlign.pl -if ' + '-d' + 'mafft_alignments/' + mygene + '_cut.aligned' + ' -ope -sn5 -b0.05') # align with transAlign
		os.system('mv ' + mygene + '_cut_tAlign.fasta' + ' ./tAlign_fastaouts/')
	else:
		os.system('mv ' +mygene + '* failed/')
			
			
		
							
print "All done!"
