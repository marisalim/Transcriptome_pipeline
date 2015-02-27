#!/usr/bin/env python

# Code to clean reads: remove low quality reads/PCR dups/optical dups (Trimmomatic), trim adapters (Trimmomatic),
# merge overlapping PE reads (FLASH, COPE), filter for contamination (use Bowtie2 to align PE reads to possible contamination source genomes)

# set path
# test data path
path = '/Volumes/Trochilidae/TestingDat/Cleanreads/Testreads'

# import modules
import os, sys, multiprocessing

# Code testing in progress ---
# let's start without multiprocessing
def trim():
	myfiles = []
	for root, dirs, files in os.walk(path):
		for filenames in files:
			if ID in filenames:
				myfiles.append(filenames)
			else:
				continue
	
	mysamps = [i.split('_', 1)[0] for i in myfiles]
	mysamps = set(mysamps)
	
	for asample in mysamps:
		variables = dict(
		adapterfile = '/Volumes/Trochilidae/TestingDat/Cleanreads/MLadapters.fa',
		trimmomatic = '/usr/local/bin/Trimmomatic-0.33/trimmomatic-0.33.jar',
		read1in = asample + '_R1.fq.gz',
		read2in = asample + '_R2.fq.gz',
		read1out = asample + '_R1_trimmed.fq.gz',
		read2out = asample + '_R2_trimmed.fq.gz',
		read1out_unpaired = asample + '_R1_trimmedunpaired.fq.gz',
		read2out_unpaired = asample + '_R2_trimmedunpaired.fq.gz',
		sampleID = asample)
		
		commands = """
		echo "These are the input files: {read1in} and {read2in}"
		echo "These are the read 1 output files: {read1out} and {read1out_unpaired}"
		echo "These are the read 2 output files: {read2out} and {read2out_unpaired}"
		
		java -classpath {trimmomatic} org.usadellab.trimmomatic.TrimmomaticPE -phred33 {read1in} {read2in} {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} ILLUMINACLIP:{adfile}:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36 LEADING:3 TRAILING:3

		""".format(**variables)
		
		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)

trim()




