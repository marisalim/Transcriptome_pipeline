#!/usr/bin/env python

# Code to merge read 1 files and read 2 files
# File structure: For each sample, there is a folder that contains the read files. 

# import modules 
import os, sys, glob, multiprocessing

# Set working directory (folder that contains the sample directory)
# Real data directory 
#path = "/Volumes/Trochilidae/TestingDat"

# Test data directory
path = "/Volumes/Trochilidae/Pythonscripts/Testpycode/testingmerge"

def mergefiles(elements):
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
	pool.map(mergefiles, myfolders)
	
if __name__ == "__main__":
	samplist()
