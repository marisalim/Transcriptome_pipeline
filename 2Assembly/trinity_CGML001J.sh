#!/bin/bash
#PBS -l ncpus=48
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N trinity_CGML001J
#-W group_list=IBN100014
set -x
source /usr/share/modules/init/bash
module load trinity/r2014-07-17
module load bowtie/1.1.1
module load samtools/0.1.19
ulimit -u unlimited
ja
cd /brashear/mlim/Trinityfiles/CGML001Jfiles2/trinityResults/CGML001J
export OMP_NUM_THREADS=32
Trinity --seqType fq --left /brashear/mlim/Trinityfiles/CGML001Jfiles2/CGML001J_trinity_left.fq --right /brashear/mlim/Trinityfiles/CGML001Jfiles2/CGML001J_trinity_right.fq --bflyGCThreads 32 --min_kmer_cov 2 --JM 345G --group_pairs_distance 999 --output /brashear/mlim/Trinityfiles/CGML001Jfiles2/trinityResults/CGML001J --CPU 48 > /brashear/mlim/Trinityfiles/CGML001Jfiles2/trinityResults/CGML001J/trinity.log
ja -chlst
