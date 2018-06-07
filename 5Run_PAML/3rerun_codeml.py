#!/usr/bin/env python

# can search files for 'check convergence', if found:
# write file name to rerun.txt
# use this rerun.txt list to search the phylipforpaml files
# grab the files to rerun and put back into paml loop
# def codeml function with input 'element' for the file name and starting values...
# in while loop that runs the function, set the different starting values
# eventually, should build this into the original script (nullbatrun or posbatrun)

# import modules
import os
from Bio.Phylo.PAML import codeml

def M2arelcodeml(element, kappastart, omegastart):	
	ctl_file = open('codemlM2arel.ctl', 'w')
	
	ctl_file.write(
	'seqfile = /crucible/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_phylipfiles/' +  element + '.phylip' + '\n' +
	'treefile = /crucible/bi4iflp/mlim/Run_codeml_mtdna/transcriptometree_noroot.tree' + '\n' +
	'outfile = ' + element + '_M2arelnull.out' + '\n' +
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
	cml.alignment = '/crucible/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_phylipfiles/' + element + '.phylip'
	cml.tree = '/crucible/bi4iflp/mlim/Run_codeml_mtdna/transcriptometree_noroot.tree'
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

#write rerun text files
outputs = [f for f in os.listdir('/home/mlim/')]
out1 = open('rerun_M2arelnull.txt', 'w')
out2 = open('rerun_cladeC.txt', 'w')
for outs in outputs:
        if 'M2arelnull.out' in outs:
            myfile = open('/home/mlim/' + outs, 'r')
            filename = outs.split('.')[0]
            for aline in myfile:
                if 'check convergence..' in aline:
#                   print filename
                    out1.write(filename + '\n')
		if 'cladeC.out' in outs:
			myfile = open('/home/mlim/' + outs, 'r')
			filename = outs.split('.')[0]
			for aline in myfile:
				if 'check convergence..' in aline:
					out2.write(filename + '\n')
out1.close()	
out2.close()

# TODO: should i make this run separately for each model or together??	
reruns = open('rerun.txt', 'r')
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
				os.system('mv ' + file_rerun + '_modAnull.out ./Codeml_outputs')
				print 'Successful rerun!!!', file_rerun, ' moved to Codeml_outputs directory'
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
				print 'If you see this message, then the file has still NOT converged. Sigh...'
				break
				
			ModelA_poscodeml(file_rerun, kappa_starts[counter], omega_starts[counter])
	
			if 'check convergence..' in open(file_rerun + '_modApos.out').read():
				print 'Nope, try again...', file_rerun
				print '---------------------------------------'
				
			else: 
				os.system('mv ' + file_rerun + '_modApos.out ./Codeml_outputs')
				print 'Successful reruns!!!', file_rerun, ' moved to Codeml_outputs directory'
				print '----------------------------------------------------------------'
				break
				
			counter += 1






		

				
		





				
		
		






		
		
		
		
		
