#!/usr/bin/env python

import os, sys

path = '/Volumes/Trochilidae/TestingDat/CleanreadsOUT'

ID = 'CGML'

def assemblereads():

	myfiles = []
	for root, dirs, files in os.walk(path):
		for filenames in files:
			if ID in filenames:
				myfiles.append(filenames)
			else:
				continue
	
	mysamps = [i.split('_', 1)[0] for i in myfiles]
	mysamps = set(mysamps)
	
	for element in mysamps:
		variables = dict(
		sample = element,
		read1 = element + '_trinity_left.fq',
		read2 = element + '_trinity_right.fq'
		) 

		# Trinity commands
		# seqType = type of reads, fa or fq
		# JM = Jellyfish Memory, number of GB of system memory to use for k-mer counting by jellyfish
		# left = (for paired reads) left reads
		# right = (for paired reads) right reads 
		# output = name of directory for output (will be created if it doesn't already exist) default: "/usr/uo/7/mlim/trinity_out_dir"
		# CPU = number of CPUs to use, default is 2
		# full_cleanup = only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
	
		commands = """
		Trinity --seqType fq --JM 10G --left {read1} --right {read2} --output {sample} --CPU 12 --full_cleanup
		""".format(**variables)

		command_list = commands.split("\n")
		for cmd in command_list:
			os.system(cmd)

	

	

