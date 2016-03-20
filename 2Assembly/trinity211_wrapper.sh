#!/bin/bash 
#PBS -l nodes=1:ppn=15 
#PBS -l walltime=24:00:00 
#PBS -j oe 
#PBS -q batch 
#PBS -m e 
#PBS -M marisa.lim@stonybrook.edu 
#PBS -N trinityB
#
set echo
module load trinity/2.1.1-all
module load bowtie/1.1.1
module load samtools/0.1.19

# run trinity
Trinity --seqType fq --left /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/CGML001B_trinity_left.fq --right /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/CGML001B_trinity_right.fq --CPU 15 --bflyGCThreads 15 --min_kmer_cov 2 --max_memory 100G --group_pairs_distance 999 --output /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/trinity_results/trinity_CGML001B > /crucible/bi4iflp/mlim/Trinity_v2.1.1_test/trinity_results/trinity_CGML001B/trinity_B.log
