#!/usr/bin/env python

# it will run codeml for null Model 1A with different starting values, checking for convergence, and grab lnL values

# NOTE: this is for nuclear dna sequences

# import modules
import os
from Bio.Phylo.PAML import codeml

# this is M1A model, the null to ModelA 
def ModelA_nullcodeml_Mphoe(element, kappastart):
	cml = codeml.Codeml()
	cml.alignment = '/pylon1/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/pylon1/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot_Mphoe.tree' #Note: using an unrooted tree
	cml.out_file = element + '_Mphoe_modAnull.out'
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

# Input files
alignmentfiles = [f for f in os.listdir('/pylon1/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/')]

# Choose start values
kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
#kappa_starts = [0,0.5] # use to test code
omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
#omega_starts = [0,0.5] # use to test code

# Run codeml (null model) with different starting values
# M1A
for aname in alignmentfiles: 
	if '.phylip' in aname:
		file_null = aname.split('.')[0]
		print '-----------------------------------------------------'
		print 'File to run null model on: ', file_null
		print '-----------------------------------------------------'
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print 'Ok, all start values have been tried.'
				print 'If you see this message, then the file has still NOT converged. Sigh...'
				break
				
			ModelA_nullcodeml_Mphoe(file_null, kappa_starts[counter])
	
			if 'check convergence..' in open(file_null + '_Mphoe_modAnull.out').read(): 
				print 'Nope, try again...', file_null
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_null + '_Mphoe_modAnull.out /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Mphoe/Mphoeout_nu')
				print 'Success!!!', file_null, ' moved to Mphoeout_nu directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1

# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml directory (or whichever directory this script is in)
os.system('grep lnL /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Mphoe/Mphoeout_nu/*modAnull.out | awk \'{print $1"\t"$5}\' > /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Mphoe_lnL_null.txt')
		
