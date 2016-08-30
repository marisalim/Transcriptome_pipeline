#!/usr/bin/env python

# it will run codeml for Model A with different starting values, checking for convergence, and grab lnL values

# NOTE: this is for nuclear dna sequences

# import modules
import os
from Bio.Phylo.PAML import codeml

# this is ModelA, allows positive selection
def ModelA_poscodeml_Ccoru(element, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = '/pylon1/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/pylon1/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot_Ccoru.tree' #Note: using an unrooted tree
	cml.out_file = element + '_Ccoru_modApos.out'
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
alignmentfiles = [f for f in os.listdir('/pylon1/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/')]

# Choose start values
kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
#kappa_starts = [0,0.5] # use to test code
omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
#omega_starts = [0,0.5] # use to test code

# Run codeml (positive selection model) with different starting values		
# Model A 
for aname in alignmentfiles:
	if '.phylip' in aname:
		file_pos = aname.split('.')[0]
		print '-----------------------------------------------------'
		print 'File to run pos model on: ', file_pos
		print '-----------------------------------------------------'
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print 'Ok, all start values have been tried.'
				print 'If you see this message, then the file has still NOT converged. Sigh...'
				break
				
			ModelA_poscodeml_Ccoru(file_pos, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_pos + '_Ccoru_modApos.out').read():
				print 'Nope, try again...', file_pos
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_pos + '_Ccoru_modApos.out /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Ccoru/Ccoru_pos/Ccoruout_nu_pos')
				print 'Success!!!', file_pos, ' moved to Ccoruout_nu_pos directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1

# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml directory (or whichever directory this script is in)
os.system('grep lnL /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Ccoru/Ccoru_pos/Ccoruout_nu_pos/*pos.out | awk \'{print $1"\t"$5}\' > /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Ccoru_lnL_pos.txt')








		



				
		
		






		
		
		
		
		
