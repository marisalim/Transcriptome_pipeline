#!/bin/bash
#PBS -l nodes=4:ppn=28,walltime=4:00:00
#PBS -N mitoB_paml
#PBS -q short

set echo

cd /gpfs/scratch/mclim
module load shared
module load anaconda/3 #note that this is using python 3
source activate biopython #to use Biopython

python mt_modA.py B codeml_trees/transcriptometree_noroot_B.tree mito_outsB/ >> mitoB_codeml.log 2>&1

