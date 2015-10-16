#!/usr/bin/env python

# This code should be run from the folder containing the data that is to be concatenated
# The results will be moved to the Concatreads folder 

# import modules
import os, sys, multiprocessing

# set path - to create the list of indexes
path = '/Volumes/Trochilidae/TestingDat/Passfilterreads'
#path = '/Volumes/Trochilidae/Pythonscripts/Testpycode/testingmerge/samp4'

# set an ID common to all files
ID = 'CGML'

# Define concatenation function
def concat(element):
#	print element
	
	variables = dict(index = element)
#	print variables
	
	# This command will concatenate read 1 files and read 2 files, then move them to a new output directory called Concatreads
	commands = """
	echo "Processing {index}"
	cat *{index}_*_R1_*.fastq.gz > {index}_R1.fq.gz
	cat *{index}_*_R2_*.fastq.gz > {index}_R2.fq.gz
	echo "Moving concatenated files to Concatreads directory: {index}"
	mv {index}*.fq.gz /Volumes/Trochilidae/TestingDat/Concatreads
	""".format(**variables)
	
	command_list = commands.split('\n')
	for cmd in command_list:
		os.system(cmd)

# This loop creates file path list as inputs for the concat function
mysamps = []
for root, dirs, files in os.walk(path):
	for filenames in files:
		if ID in filenames:
			sampname = filenames.split('_')
			sampname = sampname[1] #change to 1 for my data, 0 for the samp4 data
			mysamps.append(sampname)
	mysamps = set(mysamps)
#	print mysamps

pool = multiprocessing.Pool()
pool.map(concat, mysamps)
