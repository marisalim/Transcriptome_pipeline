#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=128GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x

module load python/2.7.11_gcc
module load paml/4.9a

# run codeml script
cd /pylon1/bi4iflp/mlim/Run_codeml_mtdna/mtDNA_codemlout_Aviri
python /pylon1/bi4iflp/mlim/Run_codeml_mtdna/mt_Aviri_modA.py
