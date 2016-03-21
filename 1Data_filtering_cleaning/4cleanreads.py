#!/usr/bin/env python

# Run trimmomatic (remove adapters, trim low quality bp, remove reads that are too short), run flash (merge PE reads), format outputs for trinity

# set path
path = '/Volumes/Trochilidae/TestingDat/Concatreads'

# import modules
import os, sys, multiprocessing

ID = 'CGML001B' # use when running full dataset

#let's try 1 sample of real dataset
#ID = 'CGML001A'

# let's start without multiprocessing
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
#	print mysamps

	for asample in mysamps:
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

		commands = """
		echo "These are the input files: {read1in} {read2in}"
		echo "These are the read1 output files: {read1out} {read1out_unpaired}"
		echo "These are the read2 output files: {read2out} {read2out_unpaired}"
		echo "This is the adapter file: {adapterfile}"
		
		# run trimmomatic
		# PE = paired end mode, -phred33 = use phred + 33 quality scores
		java -classpath {trimmomatic} org.usadellab.trimmomatic.TrimmomaticPE -phred33 {read1in} {read2in} {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} ILLUMINACLIP:{adapterfile}:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36 LEADING:3 TRAILING:3
		
		# run flash on the paired outputs from trimmomatic
		# M = maxOverlap, m = minOverlap, x = mismatchRatio, f = averageFragment Length, o = prefixofOutputFiles
		flash {read1out} {read2out} -M 100 -m 5 -x 0.05 -f 300 -o {sampleID}

		echo "Ok, FLASH is done running! Now, we will format the output files."

		# formatting flash output names
		cat {sampleID}.notCombined_1.fastq > {sampleID}_final1.fq
		cat {sampleID}.notCombined_2.fastq > {sampleID}_final2.fq
		
		# put all unpaired reads in one file
		gunzip {read1out_unpaired} 
		gunzip {read2out_unpaired}
		cat {sampleID}_R1_trimmedunpaired.fq {sampleID}_R2_trimmedunpaired.fq {sampleID}.extendedFrags.fastq > {sampleID}_finalunpaired.fq
		
		# rename the .extendedFrags.fastq file by adding /1 to the read headers (this file will be added to the other unpaired read 1 outputs for trinity_left) 
		sed -e '/\@HW/s=$=/1=' {sampleID}.extendedFrags.fastq > {sampleID}.extendedFrags_left.fastq

		# attach the unpaired reads to left/read1s (or the right/read2s, doesn't matter, but trinity can't take an unpaired reads file separately)
		echo "Make input files for Trinity!"
		cat {sampleID}.notCombined_1.fastq {sampleID}_R1_trimmedunpaired.fq {sampleID}.extendedFrags_left.fastq > {sampleID}_trinity_left.fq
		cat {sampleID}.notCombined_2.fastq {sampleID}_R2_trimmedunpaired.fq > {sampleID}_trinity_right.fq
		
		# remove files that aren't needed, but keep trim files so that I can run through FastQC - for testing purposes
		echo "File outputs are done. Time to clean up!"
#		rm {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*
		
		# remove files that aren't needed downstream - use this after testing phase done
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
trim()
