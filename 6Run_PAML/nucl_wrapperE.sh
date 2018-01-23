#!/bin/bash
#PBS -l nodes=4:ppn=28,walltime=24:00:00
#PBS -N nuclE_paml4_5
#PBS -q long

set echo

cd /gpfs/scratch/mclim/nucl_E
module load shared
module load anaconda/3 #note that this is using python 3
source activate biopython #to use Biopython

#mkdir nucl_outsE_batch1
#python /gpfs/scratch/mclim/nucl_modA.py batch1/ E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_batch1/ >> nuclE_batch1_codeml.log 2>&1

#mkdir nucl_outsE_batch2
#python /gpfs/scratch/mclim/nucl_modA.py batch2/ E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_batch2/ >> nuclE_batch2_codeml.log 2>&1
#mkdir nucl_outsE_batch3
#python /gpfs/scratch/mclim/nucl_modA.py batch3/ E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_batch3/ >> nuclE_batch3_codeml.log 2>&1
mkdir nucl_outsE_batch4
python /gpfs/scratch/mclim/nucl_modA.py batch4/ E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_batch4/ >> nuclE_batch4_codeml.log 2>&1
mkdir nucl_outsE_batch5
python /gpfs/scratch/mclim/nucl_modA.py batch5/ E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_batch5/ >> nuclE_batch5_codeml.log 2>&1
