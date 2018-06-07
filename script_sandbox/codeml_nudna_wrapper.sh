#!/bin/bash 

#PBS -l nodes=1:ppn=15 
#PBS -l walltime=24:00:00 
#PBS -j oe 
#PBS -q batch 
#PBS -m e 
#PBS -M marisa.lim@stonybrook.edu 
#PBS -N nudna_modA
#
set echo 
module load python/2.7.10-mkl 
module load paml/4.9a
# run codeml
python /crucible/bi4iflp/mlim/Run_codeml_nudna/2batruncodeml_checkconv.py
