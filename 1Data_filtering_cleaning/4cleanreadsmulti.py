#!/usr/bin/env python

# Run trimmomatic (remove adapters, trim low quality bp, remove reads that are too short), run flash (merge PE reads), format outputs for trinity

# set path
path = '/Volumes/Trochilidae/TestingDat/Concatreads'

# import modules
import os, sys, multiprocessing

ID = 'CGML'

# REMINDER: you need to source programs!

# -------------------------------------------------------------

# with multiprocessing
def trim(element):

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

	commands = """

	echo "These are the input files: {read1in} and {read2in}"
	echo "These are the read 1 output files: {read1out} and {read1out_unpaired}"
	echo "These are the read 2 output files: {read2out} and {read2out_unpaired}"

	java -classpath {trimmomatic} org.usadellab.trimmomatic.TrimmomaticPE -phred33 {read1in} {read2in} {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} ILLUMINACLIP:{adapterfile}:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36 LEADING:3 TRAILING:3

	flash {read1out} {read2out} -M 100 -m 5 -x 0.05 -f 300 -o {sampleID}
	
	echo "Ok, FLASH is done running! Now, we will format the output files."
	
	# these files are for alignment step (e.g., Bowtie2)
	cat {sampleID}.notCombined_1.fastq > {sampleID}_final1.fq
	cat {sampleID}.notCombined_2.fastq > {sampleID}_final2.fq
	gunzip {read1out_unpaired} 
	gunzip {read2out_unpaired}
	cat {sampleID}_R1_trimmedunpaired.fq {sampleID}_R2_trimmedunpaired.fq {sampleID}.extendedFrags.fastq > {sampleID}_finalunpaired.fq
	
	# renaming for trinity
	sed -e '/\@HW/s=$=/1=' {sampleID}.extendedFrags.fastq > {sampleID}.extendedFrags_left.fastq
	# formatting for trinity - trinity just takes these two files as inputs:
	cat {sampleID}.notCombined_1.fastq {sampleID}_R1_trimmedunpaired.fq {sampleID}.extendedFrags_left.fastq > {sampleID}_trinity_left.fq
	cat {sampleID}.notCombined_2.fastq {sampleID}_R2_trimmedunpaired.fq > {sampleID}_trinity_right.fq
	
	echo "File outputs are done. Time to clean up!"
	
	rm {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*

	echo "Zip files. Save space!"
	gzip {sampleID}_trinity_left.fq
	gzip {sampleID}_trinity_right.fq
	gzip {sampleID}_final1.fq
	gzip {sampleID}_final2.fq
	gzip {sampleID}_finalunpaired.fq
	
	""".format(**variables)
	
	command_list = commands.split('\n')
	for cmd in command_list:
		os.system(cmd)

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
