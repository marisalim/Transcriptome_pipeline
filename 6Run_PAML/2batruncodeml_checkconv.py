#!/usr/bin/env python

# goal here is to combine the nullbatruncodeml.py, postbatruncodeml.py, and rerun_codeml.py scripts
# it will run codeml for both models with different starting values, checking for convergence, and grab lnL values
# TODO: test that it works

# import modules
import os
from Bio.Phylo.PAML import codeml

def ModelA_nullcodeml(element, kappastart):
	cml = codeml.Codeml()
	cml.alignment = './alignments_for_codeml/phylipforpaml/' +  element + '.phylip'
	cml.tree = 'transcriptometree_root.phy'
	cml.out_file = element + '_modAnull.out'
	cml.working_dir = './'

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
	cml.alignment = './alignments_for_codeml/phylipforpaml/' + element + '.phylip'
	cml.tree = 'transcriptometree_root.phy'
	cml.out_file = element + '_modApos.out'
	cml.working_dir = './'

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
	
# Input files
alignmentfiles = [f for f in os.listdir('./alignments_for_codeml/phylipforpaml/')] ## TODO: check that this is the correct directory for inputs

# Choose start values
#kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
kappa_starts = [0,0.5] # use to test code
#omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
omega_starts = [0,0.5] # use to test code

# Run codeml (null model) with different starting values
for aname in alignmentfiles: ## TODO: check that input names correct for file_null and file_pos
	if '.phylip' in aname:
		file_null = aname.split('.')[0]
		print '-----------------------------------------------------'
		print 'File to rerun null model on: ', file_null
		print '-----------------------------------------------------'
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print 'Ok, all start values have been tried.'
				print 'If you see this message, then the file has still NOT converged. Sigh...'
				break
				
			ModelA_nullcodeml(file_null, kappa_starts[counter])
	
			if 'check convergence..' in open(file_null + '_modAnull.out').read(): 
				print 'Nope, try again...', file_null
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_null + '_modAnull.out ./Codeml_outputs')
				print 'Success!!!', file_null, ' moved to Codeml_outputs directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1
			
# Run codeml (positive selection model) with different starting values			
#for aname in alignmentfiles:
#	if '.phylip' in aname:
#		file_pos = aname.split('.')[0]
#		print '-----------------------------------------------------'
#		print 'File to rerun pos model on: ', file_pos
#		print '-----------------------------------------------------'
		
#		notconverged = True #set this to initial state for while loop
#		counter = 0 # to move consecutively through the values
		
#		while notconverged:
#			if counter > 4:
#				print 'Ok, all start values have been tried.'
#				print 'If you see this message, then the file has still NOT converged. Sigh...'
#				break
				
#			ModelA_poscodeml(file_pos, kappa_starts[counter], omega_starts[counter])
	
#			if 'check convergence..' in open(file_pos + '_modApos.out').read():
#				print 'Nope, try again...', file_pos
#				print '---------------------------------------'
				
#			else: 
#				os.system('mv ' + file_pos + '_modApos.out ./Codeml_outputs')
#				print 'Success!!!', file_pos, ' moved to Codeml_outputs directory'
#				print '----------------------------------------------------------------'
#				break
				
#			counter += 1

# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml directory (or whichever directory this script is in)
#os.system('grep lnL ./Codeml_outputs/*null.out | awk \'{print $1"\t"$5}\' > lnL_null.txt')
#os.system('grep lnL ./Codeml_outputs/*pos.out | awk \'{print $1"\t"$5}\' > lnL_pos.txt')




		

				
		





				
		
		






		
		
		
		
		
