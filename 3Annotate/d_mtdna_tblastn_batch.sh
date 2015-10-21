#!/bin/bash
#PBS -l ncpus=16
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N tblastn_trinityC_mtdna
#
set -x
source /usr/share/modules/init/bash
module load ncbi-blast/2.2.30
ja
#
# Data
#
OUTDIR=/brashear/$USER/Blast_searches/Blast_C
PROT=/brashear/$USER/Blast_searches/mtdnablasts/Trinity_C.fasta
PROTDB=/brashear/$USER/Blast_searches/mtDNA_Aversicolor.pep.fa
#
mkdir -p $OUTDIR
cd $OUTDIR
#
# Take protein file and run Blast against database, vice versa
#
date
tblastn -query $PROTDB -num_threads 16 -db $PROT -db_gencode 2 -outfmt 6 -out mtdna_dbtome_C1.out > mtdna_dbtome_C1.log 2>&1
ja -chlst

