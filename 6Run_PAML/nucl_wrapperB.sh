#!/bin/bash
#PBS -l nodes=4:ppn=28,walltime=1:00:00
#PBS -N nuclB_paml_redo
#PBS -q short

# long, 24 hours, 4 nodes for the full batches
# less time needed for debugging and running redos

set echo

cd /gpfs/scratch/mclim/nucl_B
module load shared
module load anaconda/3 #note that this is using python 3
source activate biopython #to use Biopython

#mkdir nucl_outsB_batch1
#python /gpfs/scratch/mclim/nucl_modA.py batch1/ B /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch1/ >> nuclB_batch1_codeml.log 2>&1

#mkdir nucl_outsB_batch2
#python /gpfs/scratch/mclim/nucl_modA.py batch2/ B /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch2/ >> nuclB_batch2_codeml.log 2>&1
#mkdir nucl_outsB_batch3
#python /gpfs/scratch/mclim/nucl_modA.py batch3/ B /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch3/ >> nuclB_batch3_codeml.log 2>&1
#mkdir nucl_outsB_batch4
#python /gpfs/scratch/mclim/nucl_modA.py batch4/ B /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch4/ >> nuclB_batch4_codeml.log 2>&1
#mkdir nucl_outsB_batch5
#python /gpfs/scratch/mclim/nucl_modA.py batch5/ B /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_batch5/ >> nuclB_batch5_codeml.log 2>&1

#mkdir nucl_outsB_redos
# copy these to Codeml_nucleargenes/convergence_redo_batch/
# /pos_redos/: ENSTGUP00000001031.phylip (batch5)
# /null_redos/: ENSTGUP00000013596.phylip (batch1) -- DONE
python /gpfs/scratch/mclim/nucl_redo_modA.py B /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_B.tree nucl_outsB_redos/ >> nuclB_redos_codeml.log 2>&1
