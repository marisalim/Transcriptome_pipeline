#!/usr/bin/env python

import os, sys, multiprocessing

path = '/Volumes/Trochilidae/TestingDat/CleanreadsOUT'

def assemblereads(element):

	variables = dict(
	sample = element,
	read1 = element + '_subset_trinity_left.fq',
	read2 = element + '_subset_trinity_right.fq'
	) #name your output

	## change directory
	commands = """
	/home/analysis/Downloads/trinityrnaseq_r20140413p1/Trinity --seqType fq --JM 10G --left {read1} --right {read2} --output {sample} --CPU 12 --full_cleanup
	""".format(**variables)

	command_list = commands.split("\n")
	for cmd in command_list:
		os.system(cmd)

	
# something like this ... but will have to modify!
def samplist():
	myfiles = []
	for root, dirs, files in os.walk(path):
		for filenames in files:
			if ID in filenames:
				myfiles.append(filenames)
			else:
				continue
	
	mysamps = [i.split('_', 1)[0] for i in myfiles]
	mysamps = set(mysamps)

	pool = multiprocessing.Pool()
	pool.map(assemblereads, mysamps)

if __name__ == "__main__":
	samplist()
