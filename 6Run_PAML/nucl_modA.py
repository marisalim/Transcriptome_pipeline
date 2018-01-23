#!/usr/bin/env python

# goal here is to combine the nullbatruncodeml.py, postbatruncodeml.py, and rerun_codeml.py scripts
# it will run codeml for both models with different starting values, checking for convergence, and grab lnL values

# NOTE: this is for nuclear dna sequences

# import modules
import os
import sys
from Bio.Phylo.PAML import codeml

# this is M1A model, the null to ModelA 
def ModelA_nullcodeml(element, gene_batch, tree, treename, kappastart):
	cml = codeml.Codeml()
	cml.alignment = '/gpfs/scratch/mclim/Codeml_nucleargenes/' + gene_batch + element + '.phylip'
	cml.tree = tree #Note: using an unrooted tree
	cml.out_file = element + '_' + treename + '_nucl_modAnull.out'
	cml.working_dir = './'

	cml.set_options(verbose = 1, 
	runmode = 0, # 0: user tree
	seqtype = 1, # 1: codons
	CodonFreq = 2, # 2: F3X4
	clock = 0, # 0 if unrooted tree
	aaDist = 0, # 0: equal
	model = 2, # 2:2 or more dN/dS ratios for branches
	NSsites = [2], # 2: Positive selection
	icode = 0, # 0: standard genetic code
	Mgene = 0, # 0: rates
	fix_kappa = 0, # 0: kappa estimated
	kappa = kappastart, # initial kappa
	fix_omega = 1,# 1: fix omega
	omega = 1, # 1: omega fixed to 1
	cleandata = 1) #1: yes, remove sites with ambiguity data
	
	cml.run(verbose = True)

# this is ModelA, allows positive selection
def ModelA_poscodeml(element, gene_batch, tree, treename, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = '/gpfs/scratch/mclim/Codeml_nucleargenes/' + gene_batch + element + '.phylip'
	cml.tree = tree #Note: using an unrooted tree
	cml.out_file = element + '_' + treename + '_nucl_modApos.out'
	cml.working_dir = './'

	cml.set_options(verbose = 1, 
	runmode = 0,
	seqtype = 1,
	CodonFreq = 2,
	clock = 0, # 0 if unrooted
	aaDist = 0,
	model = 2, 
	NSsites = [2],
	icode = 0,
	Mgene = 0,
	fix_kappa = 0, # 0: kappa estimated
	kappa = kappastart, # initial kappa
	fix_omega = 0,# 0: omega estimated
	omega = omegastart, # initial omega
	cleandata = 1)

	cml.run(verbose = True)

# Input files
output_dir = sys.argv[-1] #output dir name
tree = sys.argv[-2] #this is the tree file
treename = sys.argv[-3] # this is the tree name
gene_batch = sys.argv[-4]
alignmentfiles = [f for f in os.listdir('/gpfs/scratch/mclim/Codeml_nucleargenes/' + gene_batch)]

# Choose start values
kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
#kappa_starts = [0,0.5] # use to test code
omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
#omega_starts = [0,0.5] # use to test code

# Run codeml (null model) with different starting values
for aname in alignmentfiles: 
	if '.phylip' in aname:
		file_null = aname.split('.')[0]
		print('-----------------------------------------------------')
		print('File to run null model on: ', file_null)
		print('-----------------------------------------------------')
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print('Ok, all start values have been tried.')
				print('If you see this message, then the file has still NOT converged. Sigh...')
				break
				
			ModelA_nullcodeml(file_null, gene_batch, tree, treename, kappa_starts[counter])
	
			if 'check convergence..' in open(file_null + '_' + treename + '_nucl_modAnull.out').read(): 
				print('Nope, try again...', file_null)
				print('---------------------------------------')
				
			else: 
				os.system('mv ' + file_null + '_' + treename + '_nucl_modAnull.out ./' + output_dir)
				print('Success!!!', file_null, ' moved to output directory: ', output_dir)
				print('----------------------------------------------------------------')
				break
				
			counter += 1

# Run codeml (positive selection model) with different starting values		
for aname in alignmentfiles:
	if '.phylip' in aname:
		file_pos = aname.split('.')[0]
		print('-----------------------------------------------------')
		print('File to run pos model on: ', file_pos)
		print('-----------------------------------------------------')
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print('Ok, all start values have been tried.')
				print('If you see this message, then the file has still NOT converged. Sigh...')
				break
				
			ModelA_poscodeml(file_pos, gene_batch, tree, treename, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_pos + '_' + treename + '_nucl_modApos.out').read():
				print('Nope, try again...', file_pos)
				print('---------------------------------------')
				
			else: 
				os.system('mv ' + file_pos + '_' + treename + '_nucl_modApos.out ./' + output_dir)
				print('Success!!!', file_pos, ' moved to output directory: ', output_dir)
				print('----------------------------------------------------------------')
				break
				
			counter += 1
# Model A for tree 2

# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml directory (or whichever directory this script is in)

os.system('grep lnL ./' + output_dir +'/*_nucl_modAnull.out | awk \'{print $1"\t"$5}\' > ./' + output_dir + '/' + treename +'lnL_nuclnull.txt')
os.system('grep lnL ./' + output_dir +'/*_nucl_modApos.out | awk \'{print $1"\t"$5}\' > ./' + output_dir + '/' + treename +'lnL_nuclpos.txt')







		



				
		
		






		
		
		
		
		
