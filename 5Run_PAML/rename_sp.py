#!/usr/bin/env python

# Code to rename species names in fasta files for phylip conversion

import os, re

thefiles = [f for f in os.listdir('./')]

for myfile in thefiles:
	if '.fasta' in myfile:
		mygene = myfile.split('.')[0]
		myfasta = open('./' + myfile, 'r')	
		outfile = open(mygene + '_sub.fasta', 'w')
		for aline in myfasta:
			if '>' in aline:
				renamethis = aline.strip().split('.')[0] 
				aline = re.sub(aline, renamethis, aline)
				if 'a_Phlogophilus_harterti' in renamethis:
					aline = re.sub(aline, '>a_Phart', aline)
				if 'b_Chaetocercus_mulsant' in renamethis:
					aline = re.sub(aline, '>b_Cmuls', aline)
				if 'c_Phaethornis_malaris' in renamethis:
					aline = re.sub(aline, '>c_Pmala', aline)
				if 'd_Aglaeactis_castelnaudii' in renamethis:
					aline = re.sub(aline, '>d_Acast', aline)
				if 'e_Colibri_coruscans' in renamethis:
					aline = re.sub(aline, '>e_Ccoru', aline)
				if 'f_Coeligena_coeligena' in renamethis:
					aline = re.sub(aline, '>f_Ccoel', aline)
				if 'g_Coeligena_violifer' in renamethis:
					aline = re.sub(aline, '>g_Cviol', aline)
				if 'h_Amazilia_amazilia' in renamethis:
					aline = re.sub(aline, '>h_Aamaz', aline)
				if 'i_Metallura_phoebe' in renamethis:
					aline = re.sub(aline, '>i_Mphoe', aline)
				if 'j_Amazilia_viridicauda' in renamethis:
					aline = re.sub(aline, '>j_Aviri', aline)
				if 'k_Adelomyia_melanogenys' in renamethis:
					aline = re.sub(aline, '>k_Amela', aline)
				if 'l_Patagona_gigas' in renamethis:
					aline = re.sub(aline, '>l_Pgiga', aline)				
			theseq = next(myfasta).strip()
			#print '%s\n%s\n' % (aline, theseq)
			outfile.write('%s\n%s\n' % (aline, theseq))
		outfile.close()
		os.system('mv ' + mygene + '_sub.fasta ' + './Renamed_fastas')	
		



