#!/bin/bash
#PBS -l ncpus=16
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N blastx_trinityL_mtdna
#
set -x
source /usr/share/modules/init/bash
module load ncbi-blast/2.2.30
ja
#
# Data
#
OUTDIR=/brashear/$USER/Blast_searches/Blast_L
PROT=/brashear/$USER/Blast_searches/mtdnablasts/Trinity_L.fasta
PROTDB=/brashear/$USER/Blast_searches/mtDNA_Aversicolor.pep.fa
BLASTTHREADS=16
#
mkdir -p $OUTDIR
cd $OUTDIR
#
# Take protein file and run Blast against database, vice versa
#
date
blastx -query $PROT -num_threads $BLASTTHREADS -db $PROTDB -query_gencode 2 -outfmt 6 -out mtdna_metodb_L1.out > mtdna_metodb_L1.log 2>&1
ja -chlst

