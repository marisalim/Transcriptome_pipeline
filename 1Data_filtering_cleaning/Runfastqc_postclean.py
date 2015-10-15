#!/usr/bin/env python

# Code to run fastqc post trimming and post merging

# import modules
import os, sys, multiprocessing

# set path
path = '/Volumes/Trochilidae/TestingDat/Cleanreads/CleanreadsOUT'

ID = 'final'

def runfastqc_postclean(element):
	variables = dict(sample = element)
	
	commands = """
	echo "Processing {sample}"
	fastqc -o /Volumes/Trochilidae/TestingDat/Cleanreads/FastQC_postclean {sample}
	""".format(**variables)
	
	command_list = commands.split('\n')
	for cmd in command_list:
		os.system(cmd)
		
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

