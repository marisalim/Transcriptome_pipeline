#!/bin/bash 
#PBS -l nodes=1:ppn=15
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N get_cds
#
set echo
module load python/2.7.10-mkl
module load blast/2.2.31-all
#
# Get CDS - writes a file for each gene for each assembly
python /crucible/bi4iflp/mlim/Annotation/Get_CDS.py
