#!/usr/bin/env python

import os

thefiles = [f for f in os.listdir('/crucible/bi4iflp/mlim/Annotation/CDS_files/Renamed_fastas_foralignment')]

for myfile in thefiles:
	if '.fa' in myfile:
		mygene = myfile.split('_')[0]

		# make the mafft_alignments/ directory first. in > out
		os.system('mafft --localpair --adjustdirection ' + '/crucible/bi4iflp/mlim/Annotation/CDS_files/Renamed_fastas_foralignment/' + myfile + ' > ' + '/crucible/bi4iflp/mlim/Annotation/mafft_alignments/' + mygene + '.aligned')
