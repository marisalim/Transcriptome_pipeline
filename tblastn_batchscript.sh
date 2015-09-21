#!/bin/bash
#PBS -l ncpus=16
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N tblastn_trinityL
#
set -x
source /usr/share/modules/init/bash
module load ncbi-blast/2.2.30
ja
#
# Data
#
OUTDIR=/brashear/$USER/Blast_searches/Blast_L
PROT=/brashear/$USER/Blast_searches/transdecoder_cds/Trinity_L.fasta.transdecoder.cds
PROTDB=/brashear/$USER/Blast_searches/Taeniopygia_guttata.taeGut3.2.4.pep.all.fa
#
mkdir -p $OUTDIR
cd $OUTDIR
#
# Take protein file and run Blast against database, vice versa
#
date
tblastn -query $PROTDB -num_threads 16 -db $PROT -evalue 1e-10 -db_gencode 1 -outfmt 6 -out dbtome_L.out > dbtome_L.log 2>&1
ja -chlst

