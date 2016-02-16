#!/usr/bin/env python

# Code to rename species names in fasta files for alignment
# CHECK THIS FIRST: For either the nuDNA or mtDNA files, but remember to change file paths

import os, re

#thefiles = [f for f in os.listdir('/crucible/bi4iflp/mlim/Annotation/CDS_files')]
thefiles = [f for f in os.listdir('/crucible/bi4iflp/mlim/Annotation/CDS_files/Mito_CDS_files')]

for myfile in thefiles:
	if '.fa' in myfile:
		mygene = myfile.split('.')[0]
#		myfasta = open('/crucible/bi4iflp/mlim/Annotation/CDS_files/' + myfile, 'r')	
		myfasta = open('/crucible/bi4iflp/mlim/Annotation/CDS_files/Mito_CDS_files/' + myfile, 'r')
		outfile = open(mygene + '_sub.fa', 'w')
		for aline in myfasta:
			if '>' in aline:
				renamethis = aline.strip().split('.')[0] 
				aline = re.sub(aline, renamethis, aline)
				if 'Trinity_A' in renamethis:
					aline = re.sub(aline, '>a_Phlogophilus_harterti.' + mygene, aline)
				if 'Trinity_B' in renamethis:
					aline = re.sub(aline, '>b_Chaetocercus_mulsant.' + mygene, aline)
				if 'Trinity_C' in renamethis:
					aline = re.sub(aline, '>c_Phaethornis_malaris.' + mygene, aline)
				if 'Trinity_D' in renamethis:
					aline = re.sub(aline, '>d_Aglaeactis_castelnaudii.' + mygene, aline)
				if 'Trinity_E' in renamethis:
					aline = re.sub(aline, '>e_Colibri_coruscans.' + mygene, aline)
				if 'Trinity_F' in renamethis:
					aline = re.sub(aline, '>f_Coeligena_coeligena.' + mygene, aline)
				if 'Trinity_G' in renamethis:
					aline = re.sub(aline, '>g_Coeligena_violifer.' + mygene, aline)
				if 'Trinity_H' in renamethis:
					aline = re.sub(aline, '>h_Amazilia_amazilia.' + mygene, aline)
				if 'Trinity_I' in renamethis:
					aline = re.sub(aline, '>i_Metallura_phoebe.' + mygene, aline)
				if 'Trinity_J' in renamethis:
					aline = re.sub(aline, '>j_Amazilia_viridicauda.' + mygene, aline)
				if 'Trinity_K' in renamethis:
					aline = re.sub(aline, '>k_Adelomyia_melanogenys.' + mygene, aline)
				if 'Trinity_L' in renamethis:
					aline = re.sub(aline, '>l_Patagona_gigas.' + mygene, aline)				
			theseq = next(myfasta).strip()
			#print '%s\n%s\n' % (aline, theseq)
			outfile.write('%s\n%s\n' % (aline, theseq))
		outfile.close()
		os.system('mv ' + mygene + '_sub.fa ' + '/crucible/bi4iflp/mlim/Annotation/CDS_files/Renamed_fastas_foralignment')	
		


