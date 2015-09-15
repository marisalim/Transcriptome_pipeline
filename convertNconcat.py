#!/usr/bin/env python

import os
import sys
import argparse
import multiprocessing
from Bio.Nexus import Nexus

# this part defines the a function to run the fasta to nexus file conversion script
# first, get the gene name to use as output nexus file name
# then run the script to do the file conversion
def tonexus(element):
#	print element # this is the fasta input file

	genename = element.split("_")
	variables = dict(
#	thefile = './test_renamed/' + element,
	thefile = './Renamed_fastas/' + element, 
	geneID = genename[0]
	) 

#	print variables

	commands = """
	python fasta2Nexus.py {thefile} {geneID}.nexus
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)


# get your list of input fasta files
#thefiles = [f for f in os.listdir('./test_renamed')]
thefiles = [f for f in os.listdir('./Renamed_fastas')]

# make dictionary to put newly converted nexus files into
file_list = []

# here, the fasta to nexus file conversion script gets run on all the input files grabbed by thedir (should rename to thefile)
for afile in thefiles:
	if '_sub.fasta' in afile:
#		print "afile: ", afile

		tonexus(afile) # do the conversion
		file_list.append(afile.split('_')[0]+'.nexus') # take the gene name only for nexus file name


# grab the nexus files
nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list] 
#print "nexi: ", nexi

# this part concatenates all the sequences for each species
combined = Nexus.combine(nexi)
combined.write_nexus_data(filename=open('hbird_phylogeny.nexus', 'w'))

os.system('python' + ' nexus2phylip.py ' + 'hbird_phylogeny.nexus ' + 'hbird_phylogeny.phylip')


