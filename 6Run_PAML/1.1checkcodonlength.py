#!/usr/bin/env python

# Code to check if alignments are multiples of 3 (otherwise, get error in codeml)
# Then move alignments that need re-editing to Fix_sequence directory

import os, shutil

thefiles = [f for f in os.listdir('./phylipforpaml')]
print "Check codon length:"

for myfile in thefiles:
	if '.phylip' in myfile:
		mygene = myfile.split('.')[0]
		myphylip = open('./phylipforpaml/' + myfile, 'r')
		firstline = myphylip.readline()
		seqlen = firstline.split()[1]
		codonlen = float(seqlen)/3
		
#		print seqlen
#		print codonlen
#		print float(seqlen) % 3 != 0
#		print float(seqlen) % 3
	
		if float(seqlen) % 3 != 0:
			out1 = open('fix_codons.txt', 'a')
			out1.write(mygene + '\n')
	
		else:
			continue
out1.close()




