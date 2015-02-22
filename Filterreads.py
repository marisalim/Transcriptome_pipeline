#!/usr/bin/env python

# Code for removing reads that failed Illumina's quality control filters and modify sequence identifiers

# set path
path = '/Volumes/Trochilidae/Pythonscripts/Testpycode/passfilter'

# import modules
import os, sys, multiprocessing 

# an ID common to all files
ID = 'CGML'

# without multiprocessing
def passfilter(element):
	for filenames in element:
		head, tail = os.path.split(filenames)
		outname = os.path.join("test" + tail)
		outname = outname.split('.')
		outname = '.'.join(outname[:-1])
		
		variables = dict(
		index = str(filenames), 
		outfile = str(outname))
		
		commands = """
		echo "Processing {index}"
		echo "This is the outfile: {outfile}"
		gzcat {index} | grep -A 3 '^@.* [^:]*:N:[^:]:' | grep -v '^--$' | sed 's/ 1:N:0:.*/\\/1/g' | sed 's/ 2:N:0.*/\\/2/g' > {outfile}
		echo "Compressing {outfile}
		gzip {outfile}
		""".format(**variables)
		
		cmd_list = commands.split('\n')
		for cmd in cmd_list:
			os.system(cmd)
	
myfilepaths = []
for root, dirs, files in os.walk(path):
	for filename in files:
		path = os.path.join(root, filename)
		if ID in path:
			myfilepaths.append(path)
		else:
			continue
				
passfilter(myfilepaths)

# with multiprocessing
def passfilter(element):
	
	outname = element.split("/")
	outname = outname[-1]
	outname = outname[:-3]
	
	variables = dict(
	index = str(filenames), 
	outfile = str(outname))
		
	commands = """
	echo "Processing {index}"
	echo "This is the outfile: {outfile}"
	gzcat {index} | grep -A 3 '^@.* [^:]*:N:[^:]:' | grep -v '^--$' | sed 's/ 1:N:0:.*/\\/1/g' | sed 's/ 2:N:0.*/\\/2/g' > {outfile}
	echo "Compressing {outfile}
	gzip {outfile}
	""".format(**variables)
		
	cmd_list = commands.split('\n')
	for cmd in cmd_list:
		os.system(cmd)

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




