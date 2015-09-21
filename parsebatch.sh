#!/bin/bash
#PBS -l nodes=1:ppn=15
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N parse_k
#
set echo
#
module load python/2.7.10-mkl
#
# run parseBlast.py
# note: the paths need to be explicit otherwise python can't find the files
METODB=/crucible/bi4iflp/mlim/Blast_searches/Blast_K/metodb_K.out
DBTOME=/crucible/bi4iflp/mlim/Blast_searches/Blast_K/dbtome_K.out
OUTPUT=/crucible/bi4iflp/mlim/Blast_searches/Blast_K/infoK_rbh.out
# 
python /crucible/bi4iflp/mlim/Blast_searches/theblastparser.py $METODB $DBTOME $OUTPUT


