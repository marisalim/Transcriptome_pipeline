#!/usr/bin/env python

import os
import sys
import argparse
import multiprocessing
from Bio.Nexus import Nexus

# this part defines the a function to run the fasta to nexus file conversion script
# first, get the gene name to use as output nexus file name
# then run the script to do the file conversion
def tonexusphy(element):
#	print element # this is the fasta input file

	genename = element.split("_")
	variables = dict(
	thefile = element, 
	geneID = genename[0]
	) 

#	print variables

	commands = """
	python fasta2Nexus.py {thefile} {geneID}.nexus
	python nexus2phylip.py {geneID}.nexus {geneID}.phylip
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)


# get your list of input fasta files
thefiles = [f for f in os.listdir('.')]

# make dictionary to put newly converted nexus files into
file_list = []

# here, the fasta to nexus file conversion script gets run on all the input files grabbed by thedir (should rename to thefile)
for afile in thefiles:
	if '_sub.fasta' in afile:
#		print "afile: ", afile

		tonexusphy(afile) # do the conversion
		file_list.append(afile.split('_')[0]+'.phylip') # take the gene name only for nexus file name

os.system('mv ' + '*.phylip' + ' ./phylipforpaml/')
os.system('rm ' + '*.nexus')
