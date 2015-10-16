#!/usr/bin/env python

# Code to remove reads that failed Illumina's quality control filter and to rename/simplify the read names

# import modules
import os, sys, multiprocessing

# set path
path = '/Volumes/Trochilidae/TestingDat'

# set an ID common to all files
ID = 'CGML'

# Define the passfilter function
def passfilter(element):
#	print element
	
	outname = element.split('/')
	outname = outname[-1]
	outname = outname[:-3]
	outname = 'filtered_' + outname
#	print outname

	variables = dict(index = str(element), outfile = str(outname))
#	print variables

	# This command will find all the reads that passed Illumina's filter (filtering out ones that failed), and rename the sequence identifiers
	commands = """
	echo "Processing {index}"
	echo "This is the outfile: {outfile}"
	gzcat {index} | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v '^--$' | sed 's/ 1:N:0:.*/\\/1/g' | sed 's/ 2:N:0:.*/\\/2/g' > {outfile}
	echo "Compressing {outfile}"
	gzip {outfile}
	""".format(**variables)

	command_list = commands.split('\n')
	for cmd in command_list:
		os.system(cmd)

# This loop creates file path list as inputs for the passfilter function
myfilepaths = []
for root, dirs, files in os.walk(path):
	for filename in files:
		path = os.path.join(root, filename)
		if ID in path:
			myfilepaths.append(path)
		else:
			continue

pool = multiprocessing.Pool()
pool.map(passfilter, myfilepaths)
