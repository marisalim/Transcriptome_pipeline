#!/usr/bin/env python

import os
from Bio.Phylo.PAML import codeml

# can search files for 'check convergence', if found:
# write file name to rerun.txt
# use this rerun.txt list to search the phylipforpaml files
# grab the files to rerun and put back into paml loop
# def codeml function with input 'element' for the file name and starting values...
# in while loop that runs the function, set the different starting values
# eventually, should build this into the original script (nullbatrun or posbatrun)

def ModelA_nullcodeml(element, kappastart):
	cml = codeml.Codeml()
	cml.alignment = "./phylipforpaml/" + element + '.phylip'
	cml.tree = "transcriptometree_root.phy"
	cml.out_file = element + "_modAnull.out"
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
	kappa = kappastart, # initial kappa
	fix_omega = 1,# fix omega
	omega = 1, # omega fixed to 1
	cleandata = 1)
	
	cml.run(verbose = True)
	
def ModelA_poscodeml(element, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = "./phylipforpaml/" + element + '.phylip'
	cml.tree = "transcriptometree_root.phy"
	cml.out_file = element + "_modApos.out"
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
	kappa = kappastart, # initial kappa
	fix_omega = 0,# omega estimated
	omega = omegastart, # initial omega
	cleandata = 1)

	cml.set_options(cleandata = 1)
	cml.run(verbose = True)
	
reruns = open('reruntest.txt', 'r')
kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
#kappa_starts = [0,0.5] # use to test code
omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
#omega_starts = [0,0.5] # use to test code

for aname in reruns:
	if 'null' in aname:
		file_rerun = aname.split('_')[0]
		print '-----------------------------------------------------'
		print 'File to rerun null model on: ', file_rerun
		print '-----------------------------------------------------'
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print 'Ok, all start values have been tried.'
				print 'If you see this message, then the file has still NOT converged. Sigh...'
				break
				
			ModelA_nullcodeml(file_rerun, kappa_starts[counter])
	
			if 'check convergence..' in open(file_rerun + '_modAnull.out').read():
				print 'Nope, try again...', file_rerun
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_rerun + '_modAnull.out ./thereruns')
				print 'Success!!!', file_rerun, ' moved to thereruns directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1
			
	if 'pos' in aname:
		file_rerun = aname.split('_')[0]
		print '-----------------------------------------------------'
		print 'File to rerun pos model on: ', file_rerun
		print '-----------------------------------------------------'
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print 'Ok, all start values have been tried.'
				print 'If the file is not in thereruns dir, then the file has still NOT converged. Sigh...'
				break
				
			ModelA_poscodeml(file_rerun, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_rerun + '_modApos.out').read():
				print 'Nope, try again...', file_rerun
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_rerun + '_modApos.out ./thereruns')
				print 'Success!!!', file_rerun, ' moved to thereruns directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1






		

				
		





				
		
		






		
		
		
		
		
