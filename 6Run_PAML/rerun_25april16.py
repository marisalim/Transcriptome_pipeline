#!/usr/bin/env python

import os
from Bio.Phylo.PAML import codeml

# for each model, testing 3 different branch tagging schemes
# tree1: < or > 3000m
# tree2: < or > 2000m
# tree3: by HB beta-a genotype
# this is M1A model, the null to ModelA 
def ModelA_nullcodeml_tree1(element, kappastart):
	cml = codeml.Codeml()
	cml.alignment = '/crucible/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/crucible/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot.tree' #Note: using an unrooted tree
	cml.out_file = element + '_t1_modAnull.out'
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

def ModelA_nullcodeml_tree2(element, kappastart):
	cml = codeml.Codeml()
	cml.alignment = '/crucible/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/crucible/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot_sampalt.tree' #Note: using an unrooted tree
	cml.out_file = element + '_t2_modAnull.out'
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

def ModelA_nullcodeml_tree3(element, kappastart):
	cml = codeml.Codeml()
	cml.alignment = '/crucible/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/crucible/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot_hbbgeno.tree' #Note: using an unrooted tree
	cml.out_file = element + '_t3_modAnull.out'
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
def ModelA_poscodeml_tree1(element, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = '/crucible/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/crucible/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot.tree' #Note: using an unrooted tree
	cml.out_file = element + '_t1_modApos.out'
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

def ModelA_poscodeml_tree2(element, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = '/crucible/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/crucible/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot_sampalt.tree' #Note: using an unrooted tree
	cml.out_file = element + '_t2_modApos.out'
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
	
def ModelA_poscodeml_tree3(element, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = '/crucible/bi4iflp/mlim/Run_codeml_nudna/phylipforpaml/' +  element + '.phylip'
	cml.tree = '/crucible/bi4iflp/mlim/Run_codeml_nudna/transcriptometree_noroot_hbbgeno.tree' #Note: using an unrooted tree
	cml.out_file = element + '_t3_modApos.out'
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
	
# From April 2016 codeml run, there are 3 alignments that need to be rerun
# The first 2 for Model A, the last one for M1A: ENSTGUP00000005605.phylip, ENSTGUP00000001486.phylip, ENSTGUP00000016470.phylip
# using this script for other redos too
alignmentfiles = [f for f in os.listdir('/crucible/bi4iflp/mlim/Run_codeml_nudna/RERUN_phylipforpaml/')]
#print alignmentfiles

# Choose start values
kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
#kappa_starts = [0,0.5] # use to test code
omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
#omega_starts = [0,0.5] # use to test code

# Run codeml (null model) with different starting values
# M1A for tree 1
# for aname in alignmentfiles[1]:
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
				
			ModelA_nullcodeml_tree1(file_null, kappa_starts[counter])
	
			if 'check convergence..' in open(file_null + '_t1_modAnull.out').read(): 
				print 'Nope, try again...', file_null
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_null + '_t1_modAnull.out /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs')
				print 'Success!!!', file_null, ' moved to nudna_codeml_outputs directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1
# M1A for tree 2			
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
				
			ModelA_nullcodeml_tree2(file_null, kappa_starts[counter])
	
			if 'check convergence..' in open(file_null + '_t2_modAnull.out').read(): 
				print 'Nope, try again...', file_null
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_null + '_t2_modAnull.out /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_sampalt')
				print 'Success!!!', file_null, ' moved to nudna_codeml_outputs_sampalt directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1			
# M1A for tree 3
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
				
			ModelA_nullcodeml_tree3(file_null, kappa_starts[counter])
	
			if 'check convergence..' in open(file_null + '_t3_modAnull.out').read(): 
				print 'Nope, try again...', file_null
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_null + '_t3_modAnull.out /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_hbbgeno')
				print 'Success!!!', file_null, ' moved to nudna_codeml_outputs_hbbgeno directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1			
			
# take out index 1 file (only run null model), the index 0 and 2 files should be run for Model A
# with the .pop(1) statement, alignmentfiles will now just contain original index 0 and 2 files
#alignmentfiles.pop(1)
			
# Run codeml (positive selection model) with different starting values		
# Model A for tree 1	
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
				
			ModelA_poscodeml_tree1(file_pos, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_pos + '_t1_modApos.out').read():
				print 'Nope, try again...', file_pos
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_pos + '_t1_modApos.out /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs')
				print 'Success!!!', file_pos, ' moved to nudna_codeml_outputs directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1
# Model A for tree 2
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
				
			ModelA_poscodeml_tree2(file_pos, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_pos + '_t2_modApos.out').read():
				print 'Nope, try again...', file_pos
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_pos + '_t2_modApos.out /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_sampalt')
				print 'Success!!!', file_pos, ' moved to nudna_codeml_outputs_sampalt directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1
# Model A for tree 3			
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
				
			ModelA_poscodeml_tree3(file_pos, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_pos + '_t3_modApos.out').read():
				print 'Nope, try again...', file_pos
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_pos + '_t3_modApos.out /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_hbbgeno')
				print 'Success!!!', file_pos, ' moved to nudna_codeml_outputs_hbbgeno directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1			
			
# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml directory (or whichever directory this script is in)
# for tree 1
os.system('grep lnL /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs/*modAnull.out | awk \'{print $1"\t"$5}\' > /crucible/bi4iflp/mlim/Run_codeml_nudna/nu_lnL_null.txt')
os.system('grep lnL /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs/*pos.out | awk \'{print $1"\t"$5}\' > /crucible/bi4iflp/mlim/Run_codeml_nudna/nu_lnL_pos.txt')

# for tree 2
os.system('grep lnL /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_sampalt/*modAnull.out | awk \'{print $1"\t"$5}\' > /crucible/bi4iflp/mlim/Run_codeml_nudna/nu_lnL_null.txt')
os.system('grep lnL /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_sampalt/*pos.out | awk \'{print $1"\t"$5}\' > /crucible/bi4iflp/mlim/Run_codeml_nudna/nu_lnL_pos.txt')

# for tree 3
os.system('grep lnL /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_hbbgeno/*modAnull.out | awk \'{print $1"\t"$5}\' > /crucible/bi4iflp/mlim/Run_codeml_nudna/nu_lnL_null.txt')
os.system('grep lnL /crucible/bi4iflp/mlim/Run_codeml_nudna/nudna_codeml_outputs_hbbgeno/*pos.out | awk \'{print $1"\t"$5}\' > /crucible/bi4iflp/mlim/Run_codeml_nudna/nu_lnL_pos.txt')







