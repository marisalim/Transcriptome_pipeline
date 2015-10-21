#!/usr/bin/env python

# code modified from Mark to run self blast on the gene alignments files (prior to using Mafft) and trimming off unalignable parts
# then put Mafft alignments into transAlign.pl to align the amino acids
# NOTE: two settings for transAlign - 1) default and 2) lower threshold for number of allowed stop codons
# Decide which settings you want for transAlign!

# ----------------
# Import libraries
# ----------------
import os, sys
from collections import defaultdict

# -------------------
# get the input files
# -------------------
thefiles = [f for f in os.listdir('.') if os.path.isfile(f)]
#print "thefiles: ", thefiles

for myfile in thefiles:
	#print "myfiles: ", myfile

	# this takes all files in the folder and only selects one with .fasta
	if '.fasta' in myfile:
		# this takes specified files from if statement and removes the extension so you just have the file name
		mygene = myfile.split('.')[0]
#		print "mygene: ", mygene

		# -----------------
		# do the self-blast
		# -----------------
		# this assigns key = gene, value = name of mygene file
		variables = dict(gene = mygene) #names the output
#		print "variables: ", variables

		commands = """
		makeblastdb -dbtype nucl -in {gene}.fasta
		blastn -query {gene}.fasta -db {gene}.fasta -evalue 1e-10 -outfmt 6 > {gene}.out
		""".format(**variables)
		
		cmd_list = commands.split("\n")
		for cmd in cmd_list:
			os.system(cmd)

		# open self blast for next parts
		myselfblastfile = open(mygene + '.out', 'r')
		myselfblast = defaultdict(dict)
#		print myselfblast
#		print mygene

		# ---------------------------------
		# get the start and end coordinates
		# ---------------------------------
		# for each sample - find the lowest starting point and highest end point that was locally aligned (via blast) to the other samples
		# the idea is to trim away the regions that don't align as a pre-filter to using muscle or mafft (data cleaning)
		for line in myselfblastfile:
			info = line.strip().split('\t')
#			print "info: ", info
			query = info[0].split('.')[0]
			db = info[1].split('.')[0]
			start = int(info[6])
			end = int(info[7])
#			print "query = %s and db = %s ; startbp = %d and endbp = %d" % (query, db, startbp, endbp)
		
			if query == db:
				# if the query id is the same as the db id, skip and move to next line
				continue
			elif query in myselfblast.keys():
				# now we are going to assign start and end values?
				if myselfblast[query]['start'] > start:
					myselfblast[query]['start'] = start	
				else:
					myselfblast[query]['start'] = myselfblast[query]['start']
				if myselfblast[query]['end'] < end:
					myselfblast[query]['end'] = end
				else:
					myselfblast[query]['end'] = myselfblast[query]['end']
			else:
				myselfblast[query]['end'] = end
				myselfblast[query]['start'] = start

#		print "myselfblast: ", myselfblast
			
		# ---------------------------------------------
		# trim contigs if needed					
		# ---------------------------------------------
		myseqcounter = 12
		if myseqcounter == len(myselfblast.keys()): # should be 12
			thefasta = open(mygene + '.fasta', 'r')

			# make outfile
			out1 = open(mygene + '_cut.fasta', 'w') # trimmed
			
			for line in thefasta:
				# get the sequence sample id
				if '>' in line: 
					myindex = line.strip().split('.')[0][1:]
#					print "myindex: ", myindex
					seq1 = next(thefasta).strip() #next line is the sequence
					
#					print myselfblast
#					print myselfblast[myindex]['start']
#					print myselfblast[myindex]['end'] 

					seq1 = seq1[myselfblast[myindex]['start']-1:myselfblast[myindex]['end']]+ '\n'
					
					out1.write(line)
					out1.write(seq1)
			out1.close()

			# make the mafft_alignments/ directory first. in > out
			os.system('mafft --localpair --adjustdirection ' + mygene +'_cut.fasta' + ' > ' + 'mafft_alignments/' + mygene + '_cut.aligned')
			
			# make tAlign_fastaouts/ directory first. tranAlign alignment
			# align with transAlign, default settings
			os.system('transAlign.pl -if ' + '-d' + 'mafft_alignments/' + mygene + '_cut.aligned' + ' -ope') 
			# align with transAlign but change thresholds for allowed number of stop codons
			#os.system('transAlign.pl -if ' + '-d' + 'mafft_alignments/' + mygene + '_cut.aligned' + ' -ope -sn5 -b0.05')
			os.system('mv ' + mygene + '_cut_tAlign.fasta' + ' ./tAlign_fastaouts/')
		else:
			os.system('mv ' +mygene + '* failed/')
			
			
		
							
print "All done!"
