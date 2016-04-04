#!/bin/bash 
#PBS -l nodes=1:ppn=45
#PBS -l walltime=24:00:00 
#PBS -j oe 
#PBS -q batch 
#PBS -m e 
#PBS -M marisa.lim@stonybrook.edu 
#PBS -N trinityG
#
set echo
module load trinity/2.1.1-all
module load bowtie/1.1.1
module load samtools/0.1.19
ulimit -s unlimited
# run trinity
# note: bflyGCThreads has to be <32
Trinity --seqType fq --left /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/CGML001G_trinity_left.fq.gz --right /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/CGML001G_trinity_right.fq.gz --CPU 45 --bflyGCThreads 15 --min_kmer_cov 2 --max_memory 500G --group_pairs_distance 999 --output /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/trinity_results/trinity_CGML001G > /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/trinity_results/trinity_CGML001G/trinity_G.log
