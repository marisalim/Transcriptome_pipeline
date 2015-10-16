#!/usr/bin/env python

# Code to concatenate read 1 files and read 2 files
# File structure: All files in one folder called Passfilterreads/ 

# Set path - where the data and script are
path = '/Volumes/Trochilidae/TestingDat/Passfilterreads/'
#path = '/Volumes/Trochilidae/Pythonscripts/Testpycode/testingmerge/samp4' #test data

# import modules 
import os, sys, glob, multiprocessing

ID = 'CGML'

# one method 
def concatfiles(elements):
		myfolders = [root for root, dirs, files in os.walk(path)]
		del myfolders[0]
#		print myfolders
		
		for dir in myfolders: 
			currentpath = dir
			
			# Create the names for the output files 
			head, tail = os.path.split(currentpath)
			outread1 = os.path.join(currentpath, tail + "_R1.fq.gz")
			outread2 = os.path.join(currentpath, tail + "_R2.fq.gz")
#			print outread1
#			print outread2
		
			# Create list of R1 and R2 input files for each sample 
			fileread1list = [fileread1 for fileread1 in glob.glob(os.path.join(currentpath, "*R1*.fastq.gz")) if os.path.isfile(fileread1)]
			fileread2list = [fileread2 for fileread2 in glob.glob(os.path.join(currentpath, "*R2*.fastq.gz")) if os.path.isfile(fileread2)]
#			print fileread1list
#			print fileread2list
						
			# Concatenate R1 files
			with open(outread1, "w") as outfile:
				for fname in fileread1list:
					with open(fname) as infile:
#						print infile
						for line in infile:	
							outfile.write(line)	   
										
			# Concatenate R2 files
			with open(outread2, "w") as outfile:
				for fname in fileread2list:
					with open(fname) as infile:
#						print infile		
						for line in infile:
							outfile.write(line) 

def samplist():
	myfolders = [root for root, dirs, files in os.walk(path)]
	del myfolders[0]
	print myfolders
	
	pool = multiprocessing.Pool() #can specify number of processes (process = ) or leave undefined 
	pool.map(concatfiles, myfolders)
	
if __name__ == "__main__":
	samplist()
	
# another method - without multiprocessing
def concat(element):
	for samps in element:	
		variables = dict(index=samps)
		
		commands = """
		echo "Processing {index}"
		cat *{index}_*_R1_*.fastq.gz > {index}_R1.fq.gz
		cat *{index}_*_R2_*.fastq.gz > {index}_R2.fq.gz
		echo "Moving concatenated files to Concatreads directory"
		mv {index}*.fq.gz /Volumes/Trochilidae/TestingDat/Concatreads
		""".format(**variables)
		
		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)
			
mysamps = []
for root, dirs, files in os.walk(path):
	for filenames in files:
		if ID in filenames:
			sampname = filenames.split('_')
			sampname = sampname[0] # 0 for test data, 1 for real data
			mysamps.append(sampname)
	mysamps = set(mysamps) # remove duplicates

concat(mysamps)

# another method - with multiprocessing
def concat(element):
	
	variables = dict(index=element)
		
	commands = """
	echo "Processing {index}"
	cat *{index}_*_R1_*.fastq.gz > {index}_R1.fq.gz
	cat *{index}_*_R2_*.fastq.gz > {index}_R2.fq.gz
	echo "Moving concatenated files to Concatreads directory"
	mv {index}*.fq.gz /Volumes/Trochilidae/TestingDat/Concatreads
	""".format(**variables)
		
	command_list = commands.split('\n')
	for cmd in command_list:
		os.system(cmd)
			
mysamps = []
for root, dirs, files in os.walk(path):
	for filenames in files:
		if ID in filenames:
			sampname = filenames.split('_')
			sampname = sampname[1] # 0 for test data, 1 for real data
			mysamps.append(sampname)
	mysamps = set(mysamps) # remove duplicates

pool = multiprocessing.Pool()
pool.map(concat, mysamps)













	
