#!/usr/bin/env python

# Run trimmomatic (remove adapters, trim low quality bp, remove reads that are too short), run flash (merge PE reads), format outputs for trinity

# set path
path = '/Volumes/Trochilidae/TestingDat/Concatreads'

# import modules
import os, sys, multiprocessing

# set an ID common to all files
ID = 'CGML'

 print 'REMINDER: you need to source programs!'
# -------------------------------------------------------------

# Define trim and merge function, with multiprocessing
def trim(element):

	# Enter adapter file, set directory for trimmomatic, and create all the input and output names
	variables = dict(
	adapterfile = './MLadapters.fa',
	trimmomatic = '/usr/local/bin/trimmomatic-0.33/trimmomatic-0.33.jar',
	read1in = '/Volumes/Trochilidae/TestingDat/Concatreads/' + element + '_R1.fq.gz',
	read2in = '/Volumes/Trochilidae/TestingDat/Concatreads/' + element + '_R2.fq.gz',
	read1out = './CleanreadsOUT/' + element + '_R1_trimmed.fq.gz',
	read2out = './CleanreadsOUT/' + element + '_R2_trimmed.fq.gz',
	read1out_unpaired = './CleanreadsOUT/' + element + '_R1_trimmedunpaired.fq.gz',
	read2out_unpaired = './CleanreadsOUT/' + element + '_R2_trimmedunpaired.fq.gz',
	sampleID = './CleanreadsOUT/' + element)

	# These commands trim, merge, and format for Trinity
	commands = """

	echo "These are the input files: {read1in} and {read2in}"
	echo "These are the read 1 output files: {read1out} and {read1out_unpaired}"
	echo "These are the read 2 output files: {read2out} and {read2out_unpaired}"
	
	# Trimmomatic: trim reads
	java -classpath {trimmomatic} org.usadellab.trimmomatic.TrimmomaticPE -phred33 {read1in} {read2in} {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} ILLUMINACLIP:{adapterfile}:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36 LEADING:3 TRAILING:3

	# FLASH: merge reads
	# M = maxOverlap, m = minOverlap, x = mismatchRatio, f = averageFragment Length, o = prefixofOutputFiles
	flash {read1out} {read2out} -M 100 -m 5 -x 0.05 -f 300 -o {sampleID}
	
	echo "Ok, FLASH is done running! Now, we will format the output files."
	
	# Create these files - get used downstream in sequence alignment (e.g., Bowtie2 input files)
	cat {sampleID}.notCombined_1.fastq > {sampleID}_final1.fq
	cat {sampleID}.notCombined_2.fastq > {sampleID}_final2.fq
	cat {read1out_unpaired} {read2out_unpaired} {sampleID}.extendedFrags.fastq > {sampleID}_finalunpaired.fq
	
	# Rename sequence identifiers and format outputs for Trinity
	sed -e '/\@HW/s=$=/1=' {sampleID}.extendedFrags.fastq > {sampleID}.extendedFrags_left.fastq
	# formatting for Trinity - Trinity just takes these two files as inputs:
	cat {sampleID}.notCombined_1.fastq {read1out_unpaired} {sampleID}.extendedFrags_left.fastq > {sampleID}_trinity_left.fq
	cat {sampleID}.notCombined_2.fastq {read2out_unpaired} > {sampleID}_trinity_right.fq
	
	echo "File outputs are done. Time to clean up!"
	
	# USE THIS IF YOU WANT TO RUN POSTCLEAN FASTQC
#	rm {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*
	
	# USE THIS IF YOU WANT TO REMOVE ALL UNNECESSARY OUTPUTS
#	rm {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*

	""".format(**variables)
	
	command_list = commands.split('\n')
	for cmd in command_list:
		os.system(cmd)

# This loop creates file path list as inputs for the trim function
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
	pool.map(trim, mysamps)

if __name__ == "__main__":
	samplist()
