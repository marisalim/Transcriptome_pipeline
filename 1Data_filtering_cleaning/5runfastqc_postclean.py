#!/usr/bin/env python

# Code to run fastqc post trimming and post merging

# import modules
import os, sys, multiprocessing

# set path
path = '/Volumes/Trochilidae/TestingDat/Cleanreads/CleanreadsOUT'

# set an ID common to all files
ID = 'final'

# Define run FastQC function
def runfastqc_postclean(element):
	variables = dict(sample = element)
	
	# This command takes input file and runs it through FastQC, output directory is FastQC_postclean
	commands = """
	echo "Processing {sample}"
	fastqc -o /Volumes/Trochilidae/TestingDat/Cleanreads/FastQC_postclean {sample}
	""".format(**variables)
	
	command_list = commands.split('\n')
	for cmd in command_list:
		os.system(cmd)

# This loop creates file path list as inputs for the runfastqc_postclean function
def thesamps():
	myfiles = []
	for root, dirs, files in os.walk(path):
		for filenames in files:
			filepath = os.path.join(root, filenames)
			if ID in filepath:
				myfiles.append(filepath)
			else:
				continue
	
	pool = multiprocessing.Pool()
	pool.map(runfastqc_postclean, myfiles)
	
if __name__ == "__main__":
	thesamps()

