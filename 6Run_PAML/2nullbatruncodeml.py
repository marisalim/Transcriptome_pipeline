#!/usr/bin/env python

#PAML analyses - null model

## To do: add functions to run at different starting values...

import os
from Bio.Phylo.PAML import codeml

def ModelA_nullcodeml(element):
	genename = element.split(".")[0]
	cml = codeml.Codeml()
	cml.alignment = "./phylipforpaml/parttwo/" + element
	cml.tree = "transcriptometree_root.phy"
	cml.out_file = genename + "_modAnull.out"
	cml.working_dir = "./"

	cml.set_options(verbose = 1, 
	runmode = 0,
	seqtype = 1,
	CodonFreq = 2,
	clock = 0,
	aaDist = 0,
	model = 2, 
	NSsites = [2],
	icode = 0,
	Mgene = 0,
	fix_kappa = 0, # kappa estimated
	kappa = 2, # initial or fixed kappa
	fix_omega = 1,# fix omega
	omega = 1, # initial or fixed omega
	cleandata = 1)
	
	cml.run(verbose = True)
	
alignmentfiles = [f for f in os.listdir('./phylipforpaml/parttwo/')]
#alignmentfiles = [f for f in os.listdir('.')]

for afile in alignmentfiles:
	if '.phylip' in afile:
		ModelA_nullcodeml(afile)
		# move outs to the outputs directory	
		os.system('mv *null.out ./Codeml_outputs')

print 'Null models done!'		
print 'Time to find "check convergence" files...'
		
# search output for "check convergence.." so that I can pick these files out by looking at rerun.txt, the results file will be in the Run_codeml directory
outputs = [f for f in os.listdir('./Codeml_outputs')]
out1 = open('rerun_null.txt', 'w')
for outs in outputs:
        if 'null.out' in outs:
                myfile = open('./Codeml_outputs/' + outs, 'r')
                filename = outs.split('.')[0]
                for aline in myfile:
                        if 'check convergence..' in aline:
#                               print filename
                                out1.write(filename + '\n')
out1.close()

print 'Last step - grab the lnL values'

# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml directory
os.system('grep lnL ./Codeml_outputs/*null.out | awk \'{print $1"\t"$5}\' > lnL_null.txt')


		
