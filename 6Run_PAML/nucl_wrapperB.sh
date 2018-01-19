#!/bin/bash
#PBS -l nodes=4:ppn=28,walltime=48:00:00
#PBS -N nuclB_paml
#PBS -q long

set echo

cd /gpfs/scratch/mclim
module load shared
module load anaconda/3 #note that this is using python 3
source activate biopython #to use Biopython

mkdir nucl_outsB_batch1
python nucl_modA.py batch1/ B codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch1/ >> nuclB_batch1_codeml.log 2>&1

#mkdir nucl_outsB_batch2
#python nucl_modA.py batch2/ B codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch2/ >> nuclB_batch2_codeml.log 2>&1
#mkdir nucl_outsB_batch3
#python nucl_modA.py batch3/ B codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch3/ >> nuclB_batch3_codeml.log 2>&1
#mkdir nucl_outsB_batch4
#python nucl_modA.py batch4/ B codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch4/ >> nuclB_batch4_codeml.log 2>&1
#mkdir nucl_outsB_batch5
#python nucl_modA.py batch5/ B codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch5/ >> nuclB_batch5_codeml.log 2>&1
