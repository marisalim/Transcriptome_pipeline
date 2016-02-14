#!/bin/bash
#PBS -l nodes=1:ppn=15
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N mk_align
#
set echo
module load python/2.7.10-mkl
module load mafft/7.245-all
#
# Create alignments using Mafft
python /crucible/bi4iflp/mlim/Annotation/mafft_align.py

