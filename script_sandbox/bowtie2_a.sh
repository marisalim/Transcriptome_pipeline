#!/bin/bash
#PBS -l nodes=1:ppn=15
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N bt2_A
#
set echo
module load bowtie2/2.2.5
module load samtools/1.2
module load bcftools/1.2
#
# Data
#
OUTDIR=/crucible/bi4iflp/mlim/A_out
INDEXFILES=/crucible/bi4iflp/mlim/index_A/assemAindex
READ1=../Readfiles_align/CGML001A_final1.fq.gz
READ2=../Readfiles_align/CGML001A_final2.fq.gz
UNPAIR=../Readfiles_align/CGML001A_finalunpaired.fq.gz
#
mkdir -p $OUTDIR
cd $OUTDIR
#
# Run Bowtie2
#
bowtie2 -x $INDEXFILES -1 $READ1 -2 $READ2 -k 2 -S A_paired.sam
bowtie2 -x $INDEXFILES -U $UNPAIR -S A_unpaired.sam
#
# Run SAMtools and bcftools
samtools view -bS A_paired.sam > A_paired.bam
samtools view -bS A_unpaired.sam > A_unpaired.bam
samtools merge -f A_both.bam A_paired.bam A_unpaired.bam
samtools sort A_both.bam A_sorted
samtools index A_sorted.bam
samtools mpileup -d 1000000 -u -I -t DP,SP -B -f ../index_A/assemAindex.fasta A_sorted.bam | bcftools call -c - > A_snps.vcf 
