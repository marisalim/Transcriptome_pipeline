#!/usr/bin/env python

# Run trimmomatic (remove adapters, trim low quality bp, remove reads that are too short), run flash (merge PE reads), format outputs for trinity

# import modules
import os, sys, multiprocessing

# set path
path = '/Volumes/Trochilidae/TestingDat/Concatreads'

# set an ID common to all files
ID = 'CGML001E' # use when running full dataset
# let's try 1 sample of real dataset
#ID = 'CGML001A'

 print 'REMINDER: you need to source programs!'
# -------------------------------------------------------------

# Define trim and merge function, withOUT multiprocessing
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
	print mysamps

	for asample in mysamps:
	
		# Enter adapter file, set directory for trimmomatic, and create all the input and output names
		variables = dict(
		adapterfile = './MLadapters.fa', #note: because the current directory is TestingDat/Cleanreads, don't need to put full path (in fact i get errors if i do, because it adds the full path on top of the current dir path)
		trimmomatic = '/usr/local/bin/Trimmomatic-0.33/trimmomatic-0.33.jar',
		read1in = '/Volumes/Trochilidae/TestingDat/Concatreads/' + asample + '_R1.fq.gz',
		read2in = '/Volumes/Trochilidae/TestingDat/Concatreads/' + asample + '_R2.fq.gz',
		read1out = './CleanreadsOUT/' + asample + '_R1_trimmed.fq.gz',
		read2out = './CleanreadsOUT/' + asample + '_R2_trimmed.fq.gz',
		read1out_unpaired = './CleanreadsOUT/' + asample + '_R1_trimmedunpaired.fq.gz',
		read2out_unpaired = './CleanreadsOUT/' + asample + '_R2_trimmedunpaired.fq.gz', 
		sampleID = './CleanreadsOUT/' + asample)
#		print variables

		# These commands trim, merge, and format for Trinity
		commands = """
		echo "These are the input files: {read1in} {read2in}"
		echo "These are the read1 output files: {read1out} {read1out_unpaired}"
		echo "These are the read2 output files: {read2out} {read2out_unpaired}"
		echo "This is the adapter file: {adapterfile}"
		
		# Trimmomatic: trim reads
		# PE = paired end mode, -phred33 = use phred + 33 quality scores
		java -classpath {trimmomatic} org.usadellab.trimmomatic.TrimmomaticPE -phred33 {read1in} {read2in} {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} ILLUMINACLIP:{adapterfile}:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36 LEADING:3 TRAILING:3
		
		# FLASH: merge reads
		# M = maxOverlap, m = minOverlap, x = mismatchRatio, f = averageFragment Length, o = prefixofOutputFiles
		flash {read1out} {read2out} -M 100 -m 5 -x 0.05 -f 300 -o {sampleID}

		echo "Ok, FLASH is done running! Now, we will format the output files."

		# Create these files - get used downstream in sequence alignment (e.g., Bowtie2 input files)
		cat {sampleID}.notCombined_1.fastq > {sampleID}_final1E.fq
		cat {sampleID}.notCombined_2.fastq > {sampleID}_final2E.fq
		# put all unpaired reads in one file
		cat {read1out_unpaired} {read2out_unpaired} {sampleID}.extendedFrags.fastq > {sampleID}_finalunpairedE.fq
		
		# Rename sequence identifiers and format outputs for Trinity
		# rename the .extendedFrags.fastq file by adding /1 to the read headers (this file will be added to the other unpaired read 1 outputs) 
		sed -e '/\@HW/s=$=/1=' {sampleID}.extendedFrags.fastq > {sampleID}.extendedFrags_left.fastq
		# attach the unpaired reads to left/read1s (or the right/read2s, doesn't matter, but trinity can't take an unpaired reads file separately)
		cat {sampleID}.notCombined_1.fastq {read1out_unpaired} {sampleID}.extendedFrags_left.fastq > {sampleID}_trinity_leftE.fq
		cat {sampleID}.notCombined_2.fastq {read2out_unpaired} > {sampleID}_trinity_rightE.fq
		
		# Remove files that aren't needed, but keep trim files so that I can run through FastQC
		echo "File outputs are done. Time to clean up!"
		
		# USE THIS IF YOU WANT TO RUN POSTCLEAN FASTQC
#		rm {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*
		
		# USE THIS IF YOU WANT TO REMOVE ALL UNNECESSARY OUTPUTS
#		rm {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*
		
		""".format(**variables)

		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)		
trim()
