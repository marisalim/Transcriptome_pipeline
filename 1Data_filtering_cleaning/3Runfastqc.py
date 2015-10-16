#!/usr/bin/env python

# Code to run FastQC - this currently runs fastqc on the files that have been filtered and 
# concatenated, but not cleaned (not trimmed, still has adapters)

# import modules
import os, sys

# set path
path = '/Volumes/Trochilidae/TestingDat/Concatreads'

# set an ID common to all files
ID  = 'CGML'

# Define run FastQC function
def runfastqc():
	myfiles = []
	for root, dirs, files in os.walk(path):
		for filenames in files:
			if ID in filenames:
				myfiles.append(filenames)
		
	for samps in myfiles:
		variables = dict(sample = samps)
		
		# This command takes input file and runs it through FastQC, output directory is FastQC_pass1
		commands = """
		echo "Processing {sample}"
		fastqc -o /Volumes/Trochilidae/TestingDat/Cleanreads/FastQC_pass1 {sample}
		""".format(**variables)
		
		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)

runfastqc()