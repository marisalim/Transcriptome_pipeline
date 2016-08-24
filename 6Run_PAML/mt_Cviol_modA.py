# use transcriptometree_noroot_Cviol.tree

#!/usr/bin/env python

# it will run codeml for both models with different starting values, checking for convergence, and grab lnL values

# NOTE: this is for mitochondrial dna sequences

# import modules
import os
from Bio.Phylo.PAML import codeml

def ModelA_nullcodeml_Cviol(element, kappastart):
	cml = codeml.Codeml()
	cml.alignment = '/pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_phylipfiles/' +  element + '.phylip'
	cml.tree = '/pylon1/bi4iflp/mlim/Run_codeml_mtdna/transcriptometree_noroot_Cviol.tree'
	cml.out_file = element + '_Cviol_modAnull.out'
	cml.working_dir = './'

	cml.set_options(verbose = 1, 
	runmode = 0,
	seqtype = 1,
	CodonFreq = 2,
	clock = 0,
	aaDist = 0,
	model = 2, 
	NSsites = [2],
	icode = 1, #mammalian mtDNA code
	Mgene = 0,
	fix_kappa = 0, # kappa estimated
	kappa = kappastart, # initial kappa
	fix_omega = 1,# fix omega
	omega = 1, # omega fixed to 1
	cleandata = 1)
	
	cml.run(verbose = True)
	
def ModelA_poscodeml_Cviol(element, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = '/pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_phylipfiles/' + element + '.phylip'
	cml.tree = '/pylon1/bi4iflp/mlim/Run_codeml_mtdna/transcriptometree_noroot_Cviol.tree'
	cml.out_file = element + '_Cviol_modApos.out'
	cml.working_dir = './'

	cml.set_options(verbose = 1, 
	runmode = 0,
	seqtype = 1,
	CodonFreq = 2,
	clock = 0,
	aaDist = 0,
	model = 2, 
	NSsites = [2],
	icode = 1, #mammalian mtDNA code
	Mgene = 0,
	fix_kappa = 0, # kappa estimated
	kappa = kappastart, # initial kappa
	fix_omega = 0,# omega estimated
	omega = omegastart, # initial omega
	cleandata = 1)

	cml.run(verbose = True)
	
# Input files
alignmentfiles = [f for f in os.listdir('/pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_phylipfiles/')]

# Choose start values
kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
#kappa_starts = [0,0.5] # use to test code
omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
#omega_starts = [0,0.5] # use to test code

# Run codeml (null model) with different starting values
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
				
			ModelA_nullcodeml_Cviol(file_null, kappa_starts[counter])
	
			if 'check convergence..' in open(file_null + '_Cviol_modAnull.out').read(): 
				print 'Nope, try again...', file_null
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_null + '_Cviol_modAnull.out /pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_codemlout_Cviol')
				print 'Success!!!', file_null, ' moved to mtDNA_codemlout_Cviol directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1
		
# Run codeml (positive selection model) with different starting values			
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
				
			ModelA_poscodeml_Cviol(file_pos, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_pos + '_Cviol_modApos.out').read():
				print 'Nope, try again...', file_pos
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_pos + '_Cviol_modApos.out /pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_codemlout_Cviol')
				print 'Success!!!', file_pos, ' moved to mtDNA_codemlout_Cviol directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1

# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml_mtdna directory (or whichever directory this script is in)
## this runs, but saves file to home folder home/mlim
## or you can just run this separately in python after codeml analysis finishes
os.system('grep lnL /pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_codemlout_Cviol/*modAnull.out | awk \'{print $1"\t"$5}\' > /pylon1/bi4iflp/mlim/Run_codeml_mtdna/Cviol_lnL_null.txt')
os.system('grep lnL /pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_codemlout_Cviol/*pos.out | awk \'{print $1"\t"$5}\' > /pylon1/bi4iflp/mlim/Run_codeml_mtdna/Cviol_lnL_pos.txt')
