#!/usr/bin/env python 

# this code is modified from: https://github.com/hongqin/Simple-reciprocal-best-blast-hit-pairs/blob/master/RBH-v1.py
# to do: would like to make it choose RBH by specifically by e-value or percent identity. it seems to choose best hit by the bit score right now.

Usage = """RBH BLASTOUTPUT1 BLASTOUTPUT2 RBH-list-outfile """

# import modules
import sys, re

if len(sys.argv) < 3:
	print(Usage)
 
debug = 9

# input files
infl1 = sys.argv[1] #Blast result from A to B
infl2 = sys.argv[2] #Blast result from B to A
outfile = sys.argv[3] #Output file name

#parse first BLAST results
FL1 = open(infl1, 'r')
D1 = {} #dictionary for BLAST file ONE
for Line in FL1:
	if ( Line[0] != '#' ):
		Line.strip()
		Elements = re.split('\t', Line)
		queryId = Elements[0]
		subjectId = Elements[1]
		querystart = Elements[6]
		queryend = Elements[7]
		targetstart = Elements[8]
		targetend = Elements[9]
		evalue = Elements[10]
		bitscore = Elements[11]
		if ( not ( queryId in D1.keys() ) ):
			D1[queryId] = subjectId, querystart, queryend, targetstart, targetend, evalue, bitscore  #pick the first hit
#			D1[queryId] = subjectId  #pick the first hit	

if (debug): D1.keys() 

#parse second BLAST results
FL2 = open(infl2, 'r')
D2 = {}
for Line in FL2:
	if ( Line[0] != '#' ):
		Line.strip()
		Elements = re.split('\t', Line)
		queryId = Elements[0]
		subjectId = Elements[1]
		querystart = Elements[6]
		queryend = Elements[7]
		startbp = Elements[8]
		endbp = Elements[9]
		evalue = Elements[10]
		bitscore = Elements[11]
		if ( not ( queryId in D2.keys() ) ):
			D2[queryId] = subjectId, querystart, queryend, targetstart, targetend, evalue, bitscore  #pick the first hit
#			D2[queryId] = subjectId  #pick the first hit

if (debug): D2.keys() 

#print "D1: ", D1
#print "D2: ", D2

#Now, pick the shared pairs
SharedPairs={}
for id1 in D1.keys():
	value1 = D1[id1][0]
	for id2 in D2.keys():
	
#		print "if value1 == id2: ", (value1, id2)
#		print "if id1 == D2[value2]: ", (id1, D2[value1][0])
		
		if ( value1 in id2 ):
			if ( id1 == D2[value1][0] ) : #a shared best reciprocal pair
				SharedPairs[id1] = value1, D1[id1][1], D1[id1][2], D1[id1][3], D1[id1][4], D1[id1][5], D1[id1][6]
#				print "SharedPairs: ", SharedPairs
#				SharedPairs[value1] = id1
		
if (debug): SharedPairs 

# now, save the RBH results in output file
outfl = open( outfile, 'w')

for k1 in SharedPairs.keys():
	line = k1 + '\t' + SharedPairs[k1][0] + '\t' + SharedPairs[k1][1] + '\t' + SharedPairs[k1][2] + '\t' + SharedPairs[k1][3] + '\t' + SharedPairs[k1][4] + '\t' + SharedPairs[k1][5] + '\t' + SharedPairs[k1][6]
	outfl.write(line)
	
outfl.close()

print "Done. RBH from", sys.argv[1], "and", sys.argv[2], "are in", sys.argv[3]
