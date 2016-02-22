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
	
		if float(seqlen) % 3 == 0:
			continue
	
		else:
			out1 = open('fix_codons.txt', 'a')
			out1.write(mygene + '\n')
out1.close() #i think this part doesn't work if there are no sequences in out1 

print "All done!"




