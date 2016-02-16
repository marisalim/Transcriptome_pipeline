#!/usr/bin/env python

# This script is modified from Sonal Singhal: It will extract the CDS, translate sequence and check reading frame, and then write resulting sequences to file for nuclear DNA sequences and mtDNA sequences
# First, load modules on Greenfield first: module load python/2.7.10-mkl; module load blast/2.2.31-all

import re
import subprocess
import glob
import random
import os
import numpy as np 			

specieses = ['Trinity_A', 'Trinity_B', 'Trinity_C', 'Trinity_D', 'Trinity_E', 'Trinity_F', 'Trinity_G', 'Trinity_H', 'Trinity_I', 'Trinity_J', 'Trinity_K', 'Trinity_L']
#specieses = ['Trinity_A'] # test with this
#myspeciesnames = ['a_Phlogophilus_harterti'] # test with this
genes = {}
out_dir = '/crucible/bi4iflp/mlim/Annotation/CDS_files/'
mito_out_dir = '/crucible/bi4iflp/mlim/Annotation/CDS_files/Mito_CDS_files/'
prot_db = '/crucible/bi4iflp/mlim/Annotation/Taeniopygia_guttata.taeGut3.2.4.pep.all.fa' 

def get_protein(prot_db):
	protein = {}
	id = ''
	f = open(prot_db, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			protein[id] = ''
		else:
			seq = l.rstrip()
			seq = seq.replace('*', '')
			protein[id] += seq
	f.close()
	return protein

def get_seqs(species, genes):
	f = open('/crucible/bi4iflp/mlim/Trinity_assembly_files/%s/%s.fasta.annotated' % (species, species), 'r')
	print "This is the annotated file: ", f
	for l in f:
		if re.search('>', l):
			d = re.split('\t', l.rstrip())
			#print "This is the header line: ", d # looks like this: this is header line:  ['>contig9974', '5u150_gs151_ge435_3u436', 'ENSTGUP00000000263', 'HIG1 hypoxia inducible domain family, member 1A', '5e-49']

			if len(d) > 1:
				gene = d[2] # this is the gene start and end info part
				seq = f.next().rstrip()
				cds_start = int(re.search('gs(\d+)', d[1]).group(1))
				cds_end = int(re.search('ge(\d+)', d[1]).group(1))
				#print "This is cds_start: ", cds_start # looks like this: This is cds_start:  151
				#print "This is cds_end: ", cds_end # looks like this:  This is cds_end:  435

				cds = seq[(cds_start - 1):cds_end]
				#print "This is resulting cds: ", cds # looks like this (this is just random example, not same as above d, start, or end): This is resulting cds:  GTTGCAGGGGCAGGTGTGCTGACAACATTGCTGCTGATGCTGATTTTGCTGGTCAGGCTGCCCTTCATCAAGGACAAGGAGAAGAAAAGCCCCCTGGGAATGCACTTCCTCTTTCTCTTGGGGACTCTGGGACTGTTTGGGCTGACGTTTGCTTTTATCATCCAGGAAGATGAGACGGTGTGCTCTGCCCGAAGGTTTCTCTGGGGAGTTATC

				if gene not in genes:
					genes[gene] = {}
				genes[gene][species] = cds
	return genes

# this is for nuclear DNA
def translate(seq):
	map = {	'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    		'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
  		'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    		'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
    		'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    		'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    		'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    		'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    		'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
   		'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    		'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    		'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
    		'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    		'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    		'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    		'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

	triplets = [seq[i:i+3] for i in range(0, len(seq), 3)]
	aa = ''
	for triplet in triplets:
		if triplet in map:
			aa += map[triplet]
		else:
			aa += 'X'
	return aa

def nu_check_frame(prots, seq, gene):
	for i in range(0,3):
		# translate that frame
		aa = translate(seq[i:])
		# do the blast
		query = 'query.fa'
		q = open(query, 'w')
		q.write('>query\n%s\n' % aa)
		q.close()
		subject = 'subject.fa'
                s = open(subject, 'w')
                s.write('>subject\n%s\n' % prots[gene])
                s.close()

		blast = subprocess.Popen('blastp -query %s -subject %s -outfmt 6' % (query, subject), shell=True, stdout=subprocess.PIPE)
		top_match = 0 
		for l in blast.stdout:
			d = re.split('\t', l.rstrip())
			if float(d[2]) > top_match:
				top_match = float(d[2])
		
		if (top_match > 30):
			return i
	return None	
	
def write_files(prots, out_dir, genes):
	to_align_files = []
	for gene in genes:
		# Note: gene = the Ensembl peptide ID, genes = Trinity_*: the sequence
		# are the focal species in this gene?
		# Note: if you write one species here, you can get genes for just that species, but if you want to align all 12 species, you have to write all of the names here in the if statement
		# BUT: this statement grabs anything with any of these labels...not just the ones with all 12...
		if 'Trinity_A' in genes[gene] and 'Trinity_B' in genes[gene] and 'Trinity_C' in genes[gene] and 'Trinity_D' in genes[gene] and 'Trinity_E' in genes[gene] and 'Trinity_F' in genes[gene] and 'Trinity_G' in genes[gene] and 'Trinity_H' in genes[gene] and 'Trinity_I' in genes[gene] and 'Trinity_J' in genes[gene] and 'Trinity_K' in genes[gene] and 'Trinity_L' in genes[gene]: 
			print 'Yes, 12 species for %s' % gene
			
			max_length = np.max([len(x) for x in genes[gene].values()])
			out_file = '%s/%s.fa' % (out_dir, gene)
			o = open(out_file, 'w')

			for species in genes[gene]:
				# check frame
				frame = nu_check_frame(prots, genes[gene][species], gene)
				# check length
				if frame != None:
					#print "This is frame: ", frame # what does 0 mean?
					aa = translate(genes[gene][species][frame:])
					len_seq = len(re.search('^([^\*]+)', aa).group(1)) * 3
					if (len_seq / float(max_length)) > 0.8:				
						seq = genes[gene][species][frame:(frame+len_seq)]
						o.write('>%s.%s\n%s\n' % (species, gene, seq))
			o.close()
			to_align_files.append(out_file) # this appends results, so if you rerun, delete old results first
			
		else:
			print 'No, fewer than 12 species for %s' % gene
			continue
			
	return to_align_files		

# this is for mtDNA	
def mito_translate(seq):
	mito_map = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    		'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
  		'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
    		'TGT':'C', 'TGC':'C', 'TGA':'W', 'TGG':'W',
    		'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    		'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    		'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
    		'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
    		'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M',
   		'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    		'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
    		'AGT':'S', 'AGC':'S', 'AGA':'*', 'AGG':'*',
    		'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    		'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    		'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
    		'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
	
	mito_triplets = [seq[i:i+3] for i in range(0, len(seq), 3)]
	mito_aa = ''
	for mito_triplet in mito_triplets:
		if mito_triplet in mito_map:
			mito_aa += mito_map[mito_triplet]
		else:
			mito_aa += 'X'
	return mito_aa
		
def mito_check_frame(prots, seq, gene):
	for i in range(0,3):
		# translate that frame
		mito_aa = mito_translate(seq[i:])
		# do the blast
		query = 'query.fa'
		q = open(query, 'w')
		q.write('>query\n%s\n' % mito_aa)
		q.close()
		subject = 'subject.fa'
                s = open(subject, 'w')
                s.write('>subject\n%s\n' % prots[gene])
                s.close()

		blast = subprocess.Popen('blastp -query %s -subject %s -outfmt 6' % (query, subject), shell=True, stdout=subprocess.PIPE)
		top_match = 0 
		for l in blast.stdout:
			d = re.split('\t', l.rstrip())
			if float(d[2]) > top_match:
				top_match = float(d[2])
		
		if (top_match > 30):
			return i
	return None
				
def mito_write_files(prots, mito_out_dir, genes):
	mito_to_align_files = []
	for gene in genes:
		# are the focal species in this gene?
		if 'Trinity_A' in genes[gene] and 'Trinity_B' in genes[gene] and 'Trinity_C' in genes[gene] and 'Trinity_D' in genes[gene] and 'Trinity_E' in genes[gene] and 'Trinity_F' in genes[gene] and 'Trinity_G' in genes[gene] and 'Trinity_H' in genes[gene] and 'Trinity_I' in genes[gene] and 'Trinity_J' in genes[gene] and 'Trinity_K' in genes[gene] and 'Trinity_L' in genes[gene]: 
#			print 'Yes, 12 species for %s' % gene
			
			max_length = np.max([len(x) for x in genes[gene].values()])
			out_file = '%s%s.fa' % (mito_out_dir, gene)
			o = open(out_file, 'w')

			for species in genes[gene]:
				# check frame
				frame = mito_check_frame(prots, genes[gene][species], gene)
				# check length
				if frame != None:
					#print "This is frame: ", frame # what does 0 mean?
					mito_aa = mito_translate(genes[gene][species][frame:])

					seq_search = re.search('^([^\*]+)', mito_aa)
					# need to add this because sometimes the search comes up with no result (NoneType)
					if seq_search != None:
						len_seq = len(seq_search.group(1))*3
					
						if (len_seq / float(max_length)) > 0.8:	
							seq = genes[gene][species][frame:(frame+len_seq)]
							o.write('>%s.%s\n%s\n' % (species, gene, seq))
			o.close()
			mito_to_align_files.append(out_file) # this appends results, so if you rerun, delete old results first
			
		else:
#			print 'No, fewer than 12 species for %s' % gene
			continue
			

	return mito_to_align_files				
			
# Generate inputs for final write_files()	
prots = get_protein(prot_db)
for species in specieses:
	genes = get_seqs(species, genes)

# Write files for nuclear genes
to_align_files = write_files(prots, out_dir, genes)

# Write files for mitochondrial genes
mito_to_align_files = mito_write_files(prots, mito_out_dir, genes)

print "All functions completed yay yay yay! Ok, now let's rename the Trinity_* to species names!"
