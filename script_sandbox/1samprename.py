#!/usr/bin/env python

# code to rename sequence ids in fasta files

import os, re

thefiles = [f for f in os.listdir('./Filestorename')]
#print thefiles

for myfile in thefiles:
	if '.fasta' in myfile:
		mygene = myfile.split('_')[0]
		myfasta = open('./Filestorename/' + myfile, 'r')
		out1 = open(mygene + '_sub.fasta', 'w')
		for aline in myfasta:
			if '>' in aline:
				mysampname = aline.strip().split('.')[0] 
				aline = re.sub(aline, mysampname, mysampname) + '\n'	
				if '_R_' in mysampname:
					aline = re.sub('_R_', '', mysampname) + '\n'
			seq1 = next(myfasta).strip() + '\n'
			out1.write(aline)
			out1.write(seq1)
		out1.close()
		os.system('mv ' + mygene + '_sub.fasta' + ' ./Renamed_fastas')	
