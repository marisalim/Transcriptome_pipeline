#!/usr/bin/env python

# Code to clean reads: remove low quality reads/PCR dups/optical dups (Trimmomatic), trim adapters (Trimmomatic),
# merge overlapping PE reads (FLASH, COPE), filter for contamination (use Bowtie2 to align PE reads to possible contamination source genomes)

# set path
# test data path
#path = '/Volumes/Trochilidae/TestingDat/Cleanreads/Testreads'

# import modules
import os, sys, multiprocessing

ID = 'CGML'

# Code testing in progress ---
# let's start without multiprocessing
# def trim():
	# myfiles = []
	# for root, dirs, files in os.walk(path):
		# for filenames in files:
			# if ID in filenames:
				# myfiles.append(filenames)
			# else:
				# continue
	
	# mysamps = [i.split('_', 1)[0] for i in myfiles]
	# mysamps = set(mysamps)
	
	# for asample in mysamps:
		# variables = dict(
		# adapterfile = '/Volumes/Trochilidae/TestingDat/Cleanreads/MLadapters.fa',
		# trimmomatic = '/usr/local/bin/Trimmomatic-0.33/trimmomatic-0.33.jar',
		# read1in = asample + '_testclean_R1.fq.gz',
		# read2in = asample + '_testclean_R2.fq.gz',
		# read1out = './TestcleanOUT/' + asample + '_R1_trimmed.fq.gz',
		# read2out = './TestcleanOUT/' + asample + '_R2_trimmed.fq.gz',
		# read1out_unpaired = './TestcleanOUT/' + asample + '_R1_trimmedunpaired.fq.gz',
		# read2out_unpaired = './TestcleanOUT/' + asample + '_R2_trimmedunpaired.fq.gz',
		# sampleID = './TestcleanOUT/' + asample)
		
		# commands = """
		# echo "These are the input files: {read1in} and {read2in}"
		# echo "These are the read 1 output files: {read1out} and {read1out_unpaired}"
		# echo "These are the read 2 output files: {read2out} and {read2out_unpaired}"
		
		# java -classpath {trimmomatic} org.usadellab.trimmomatic.TrimmomaticPE -phred33 {read1in} {read2in} {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} ILLUMINACLIP:{adfile}:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36 LEADING:3 TRAILING:3

		# flash {read1out} {read2out} -M 100 -m 5 -x 0.05 -f 300 -o {sampleID}
	
		# cat {sampleID}.notCombined_1.fastq > {sampleID}_final1.fq
		# cat {sampleID}.notCombined_2.fastq > {sampleID}_final2.fq

		# cat {read1out_unpaired} {read2out_unpaired} {sampleID}.extendedFrags.fastq > {sampleID}_finalunpaired.fq
		
		# sed -e '/\@HW/s=$=/1=' {sampleID}.extendedFrags.fastq > {sampleID}.extendedFrags_left.fastq
	
		# cat {sampleID}.notCombined_1.fastq {read1out_unpaired} {sampleID}.extendedFrags_left.fastq > {sampleID}_trinity_left.fq
		# cat {sampleID}.notCombined_2.fastq {read2out_unpaired} > {sampleID}_trinity_right.fq
		# rm {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*

		# """.format(**variables)
		
		# command_list = commands.split('\n')
		# for cmd in command_list:
			# os.system(cmd)

# trim()

# Code for my data - just need to change the output files
path = '/Volumes/Trochilidae/TestingDat/Concatreads'

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
	
	for asample in mysamps:
		variables = dict(
        adapterfile = './MLadapters.fa', #note: because the current directory is TestingDat/Cleanreads, don't need to put full path (in fact i get e$
        trimmomatic = '/usr/local/bin/Trimmomatic-0.33/trimmomatic-0.33.jar',
        read1in = '/Volumes/Trochilidae/TestingDat/Concatreads/' + asample + '_R1.fq.gz',
        read2in = '/Volumes/Trochilidae/TestingDat/Concatreads/' + asample + '_R2.fq.gz',
        read1out = './CleanreadsOUT/' + asample + '_R1_trimmed.fq.gz',
        read2out = './CleanreadsOUT/' + asample + '_R2_trimmed.fq.gz',
        read1out_unpaired = './CleanreadsOUT/' + asample + '_R1_trimmedunpaired.fq.gz',
        read2out_unpaired = './CleanreadsOUT/' + asample + '_R2_trimmedunpaired.fq.gz',
        sampleID = './CleanreadsOUT/' + asample)
		
		commands = """
		echo "These are the input files: {read1in} and {read2in}"
		echo "These are the read 1 output files: {read1out} and {read1out_unpaired}"
		echo "These are the read 2 output files: {read2out} and {read2out_unpaired}"
		
		java -classpath {trimmomatic} org.usadellab.trimmomatic.TrimmomaticPE -phred33 {read1in} {read2in} {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} ILLUMINACLIP:{adfile}:2:40:15 SLIDINGWINDOW:4:20 MINLEN:36 LEADING:3 TRAILING:3

		flash {read1out} {read2out} -M 100 -m 5 -x 0.05 -f 300 -o {sampleID}
	
		echo "Ok, FLASH is done running! Now, we will format the output files."
	
		cat {sampleID}.notCombined_1.fastq > {sampleID}_final1.fq
		cat {sampleID}.notCombined_2.fastq > {sampleID}_final2.fq

		cat {read1out_unpaired} {read2out_unpaired} {sampleID}.extendedFrags.fastq > {sampleID}_finalunpaired.fq
		
		sed -e '/\@HW/s=$=/1=' {sampleID}.extendedFrags.fastq > {sampleID}.extendedFrags_left.fastq
	
		cat {sampleID}.notCombined_1.fastq {read1out_unpaired} {sampleID}.extendedFrags_left.fastq > {sampleID}_trinity_left.fq
		cat {sampleID}.notCombined_2.fastq {read2out_unpaired} > {sampleID}_trinity_right.fq
		
		echo "File outputs are done. Time to clean up!"
		rm {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*

		""".format(**variables)
		
		command_list = commands.split('\n')
		for cmd in command_list:
			os.system(cmd)

trim()

# with multiprocessing

def trim(element):
	
	variables = dict(
	adapterfile = './MLadapters.fa', #note: because the current directory is TestingDat/Cleanreads, don't need to put full path (in fact i get e$
    trimmomatic = '/usr/local/bin/Trimmomatic-0.33/trimmomatic-0.33.jar',
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
	cat {sampleID}.notCombined_1.fastq > {sampleID}_final1.fq
	cat {sampleID}.notCombined_2.fastq > {sampleID}_final2.fq

	cat {read1out_unpaired} {read2out_unpaired} {sampleID}.extendedFrags.fastq > {sampleID}_finalunpaired.fq
	
	sed -e '/\@HW/s=$=/1=' {sampleID}.extendedFrags.fastq > {sampleID}.extendedFrags_left.fastq
	
	cat {sampleID}.notCombined_1.fastq {read1out_unpaired} {sampleID}.extendedFrags_left.fastq > {sampleID}_trinity_left.fq
	cat {sampleID}.notCombined_2.fastq {read2out_unpaired} > {sampleID}_trinity_right.fq
	echo "File outputs are done. Time to clean up!"
	rm {read1out} {read1out_unpaired} {read2out} {read2out_unpaired} {sampleID}.extendedFrags_left.fastq {sampleID}.extendedFrags.fastq {sampleID}.notCombined_1.fastq {sampleID}.notCombined_2.fastq {sampleID}*hist*

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


