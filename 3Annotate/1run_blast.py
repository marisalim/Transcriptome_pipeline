#!/usr/bin/env python

# Code to run blastx and tblastn on nuclear DNA and mtDNA
# TODO: test run this code! I just wrote it combining the commands from the 4 blast .sh batch files, but haven't tested it yet

# import modules
import os, sys

# set path
path = '/brashear/$USER/Blast_searches'

def blastnuc():
	myfiles = []
	for root, dirs, files in os.walk(path):
		for filenames in files:
			if ID in filenames:
				myfiles.append(filenames)
			else:
				continue
	mysamps = [i.split('.', 1)[0] for i in myfiles] #TO DO: want to get the 'Trinity_L' label
	mysamps = set(mysamps)
	print mysamps
	
	for asample in mysamps:
	
		samplett = [asample.split('_', 1)[1]
		samplett = set(samplett) ## do i need to do this?
	
		# Enter blast query, blast database, and output directory
		variables = dict(
		blastthreads = 16,
		blastquery = './' + asample + '.fasta.transdecoder.cds',
		blastdb = './Taeniopygia_guttata.taeGut3.2.4.pep.all.fa',
		outdir = './Blast_' + samplett,
		blastx_out = './outdir/metodb_' + samplett + '.out',
		blastx_log = './outdir/metodb_' + samplett + '.log',
		tblastn_out = './outdir/dbtome_' + samplett + '.out',
		tblastn_log = './outdir/dbtome_' + samplett + '.log')
#		print variables

		commands = """
		blastx -query blastquery -num_threads blastthreads -db blastdb -evalue 1e-10 -query_gencode 1 -outfmt 6 -out blastx_out > blastx_log 2>&1
		tblastn -query blastdb -num_threads blastthreads -db blastquery -evalue 1e-10 -db_gencode 1 -outfmt 6 -out tblastn_out > tblastn_log 2>&1
		""".format(**variables)
		
		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)
blastnuc()

def blastmtdna():
	myfiles = []
	for root, dirs, files in os.walk(path):
		for filenames in files:
			if ID in filenames:
				myfiles.append(filenames)
			else:
				continue
	mysamps = [i.split('.', 1)[0] for i in myfiles] #TO DO: want to get the 'Trinity_L' label
	mysamps = set(mysamps)
	print mysamps
	
	for asample in mysamps:
	
		samplett = [asample.split('_', 1)[1]
		samplett = set(samplett) ## do i need to do this?
	
		# Enter blast query, blast database, and output directory
		variables = dict(
		blastthreads = 16,
		blastquery = './mtdnablasts/' + asample + '.fasta',
		blastdb = './mtDNA_Aversicolor.pep.fa',
		outdir = './Blast_' + samplett,
		blastx_out = './outdir/mtdna_metodb_' + samplett + '1.out',
		blastx_log = './outdir/mtdna_metodb_' + samplett + '1.log',
		tblastn_out = './outdir/mtdna_dbtome_' + samplett + '1.out',
		tblastn_log = './outdir/mtdna_dbtome_' + samplett + '1.log')
#		print variables

		commands = """
		blastx -query blastquery -num_threads blastthreads -db blastdb -query_gencode 2 -outfmt 6 -out blastx_out > blastx_log 2>&1
		tblastn -query blastdb -num_threads blastthreads -db blastquery -db_gencode 2 -outfmt 6 -out tblastn_out > tblastn_log 2>&1		
		""".format(**variables)
		
		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)
blastmtdna()