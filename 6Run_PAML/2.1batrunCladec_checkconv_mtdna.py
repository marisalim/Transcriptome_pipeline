#!/usr/bin/env python

# goal here is to combine the nullbatruncodeml.py, postbatruncodeml.py, and rerun_codeml.py scripts
# it will run codeml for both models with different starting values, checking for convergence, and grab lnL values

# NOTE: this is for mitochondrial dna sequences

# import modules
import os
from Bio.Phylo.PAML import codeml

def M2arelcodeml(element, kappastart, omegastart):	
	#cml = codeml.Codeml()
	#cml.alignment = '../For_codeml/mtDNA_phylipfiles/' +  element + '.phylip'
	#cml.tree = 'transcriptometree_noroot.tree'
	#cml.out_file = element + '_M2arelnull.out'
	#cml.working_dir = './'

	#cml.set_options(verbose = 1, 
#	runmode = 0,
#	seqtype = 1,
#	CodonFreq = 2,
#	clock = 0,
#	aaDist = 0,
#	model = 0, 
#	NSsites = [22], ##FIXME! this is not a siteclass model in Bio.Phylo.PAML _parse_codeml.py 
#	icode = 1, #mammalian mtDNA code
#	Mgene = 0,
#	fix_kappa = 0, # kappa estimated
#	kappa = kappastart, 
#	fix_omega = 0,# omega estimated
#	omega = omegastart, 
#	cleandata = 1)
	
#	cml.run(verbose = True)

	ctl_file = open('codemlM2arel.ctl', 'w')
	
	ctl_file.write(
	'seqfile = ../For_codeml/mtDNA_phylipfiles/' +  element + '.phylip' + '\n' +
	'treefile = transcriptometree_noroot.tree' + '\n' +
	'outfile =' + element + '_M2arelnull.out' + '\n' +
	'verbose = 1' + '\n' +
	'runmode = 0' + '\n' +
	'seqtype = 1' + '\n' +
	'CodonFreq = 2' + '\n' +
	'clock = 0' + '\n' +
	'aaDist = 0' + '\n' +
	'model = 0' + '\n' +
	'NSsites = 22' + '\n' +
	'icode = 1' + '\n' + #mammalian mtDNA code
	'Mgene = 0' + '\n' +
	'fix_kappa = 0' + '\n' + # kappa estimated
	'kappa = ' + str(kappastart) + '\n' + 
	'fix_omega = 0' + '\n' + # omega estimated
	'omega = ' + str(omegastart) + '\n' +
	'cleandata = 1')
	
	ctl_file.close()
		
def cladeCcodeml(element, kappastart, omegastart):
	cml = codeml.Codeml()
	cml.alignment = '../For_codeml/mtDNA_phylipfiles/' + element + '.phylip'
	cml.tree = 'transcriptometree_noroot.tree'
	cml.out_file = element + '_cladeC.out'
	cml.working_dir = './'

	cml.set_options(verbose = 1, 
	runmode = 0,
	seqtype = 1,
	CodonFreq = 2,
	clock = 0,
	aaDist = 0,
	model = 3, 
	NSsites = [2],
	icode = 1, #mammalian mtDNA code
	Mgene = 0,
	fix_kappa = 0, # kappa estimated
	kappa = kappastart, 
	fix_omega = 0,# omega estimated
	omega = omegastart, 
	cleandata = 1)
	
	cml.run(verbose = True)

# Input files
alignmentfiles = [f for f in os.listdir('../For_codeml/mtDNA_phylipfiles/')]

# Choose start values
#kappa_starts = [0, 0.5, 1, 3.5, 15] # range of start values for kappa
kappa_starts = [0,0.5] # use to test code
#omega_starts = [0, 0.5, 1, 3.5, 15] # range of start values for omega
omega_starts = [0,0.5] # use to test code

# Run codeml (null model) with different starting values
for aname in alignmentfiles: 
	if '.phylip' in aname:
		file_M2arel = aname.split('.')[0]
		print '-----------------------------------------------------'
		print 'File to run M2arel model on: ', file_M2arel
		print '-----------------------------------------------------'
		
		notconverged = True #set this to initial state for while loop
		counter = 0 # to move consecutively through the values
		
		while notconverged:
			if counter > 4:
				print 'Ok, all start values have been tried.'
				print 'If you see this message, then the file has still NOT converged. Sigh...'
				break
				
			M2arelcodeml(file_M2arel, kappa_starts[counter], omega_starts[counter])
			os.system('codeml /Volumes/Trochilidae/TestingDat/Feb16_makephylipforcodeml/Run_codeml/codemlM2arel.ctl')
	
			if 'check convergence..' in open(file_M2arel + '_M2arelnull.out').read(): 
				print 'Nope, try again...', file_M2arel
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_M2arel + '_M2arelnull.out ./mtDNA_codeml_outputs')
				print 'Success!!!', file_M2arel, ' moved to mtDNA_codeml_outputs directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1
		
# Run codeml (positive selection model) with different starting values			
#for aname in alignmentfiles:
#	if '.phylip' in aname:
#		file_cladeC = aname.split('.')[0]
#		print '-----------------------------------------------------'
#		print 'File to run cladeC model on: ', file_cladeC
#		print '-----------------------------------------------------'
		
#		notconverged = True #set this to initial state for while loop
#		counter = 0 # to move consecutively through the values
		
#		while notconverged:
#			if counter > 4:
#				print 'Ok, all start values have been tried.'
#				print 'If you see this message, then the file has still NOT converged. Sigh...'
#				break
				
#			cladeCcodeml(file_cladeC, kappa_starts[counter], omega_starts[counter])
	
#			if 'check convergence..' in open(file_cladeC + '_cladeC.out').read():
#				print 'Nope, try again...', file_cladeC
#				print '---------------------------------------'
				
#			else: 
#				os.system('mv ' + file_cladeC + '_cladeC.out ./mtDNA_codeml_outputs')
#				print 'Success!!!', file_cladeC, ' moved to mtDNA_codeml_outputs directory'
#				print '----------------------------------------------------------------'
#				break
				
#			counter += 1

# grab the lnL outputs to check if results significant, the results file will be in the Run_codeml_mtdna directory (or whichever directory this script is in)
## this runs, but saves file to home folder home/mlim
## or you can just run this separately in python after codeml analysis finishes
#os.system('grep lnL ./mtDNA_codeml_outputs/*M2arelnull.out | awk \'{print $1"\t"$5}\' > lnL_M2a_relnull.txt')
#os.system('grep lnL ./mtDNA_codeml_outputs/*cladeC.out | awk \'{print $1"\t"$5}\' > lnL_cladeC.txt')
