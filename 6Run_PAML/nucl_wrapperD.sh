#!/bin/bash
#PBS -l nodes=4:ppn=28,walltime=2:00:00
#PBS -N nuclD_paml_redos
#PBS -q short

# long, 24 hours, 4 nodes for the full batches
# less time needed for debugging and running redos

set echo

cd /gpfs/scratch/mclim/nucl_D
module load shared
module load anaconda/3 #note that this is using python 3
source activate biopython #to use Biopython

#mkdir nucl_outsD_batch1
#python /gpfs/scratch/mclim/nucl_modA.py batch1/ D /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_D.tree nucl_outsD_batch1/ >> nuclD_batch1_codeml.log 2>&1

#mkdir nucl_outsD_batch2
#python /gpfs/scratch/mclim/nucl_modA.py batch2/ D /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_D.tree nucl_outsD_batch2/ >> nuclD_batch2_codeml.log 2>&1
#mkdir nucl_outsD_batch3
#python /gpfs/scratch/mclim/nucl_modA.py batch3/ D /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_D.tree nucl_outsD_batch3/ >> nuclD_batch3_codeml.log 2>&1
#mkdir nucl_outsD_batch4
#python /gpfs/scratch/mclim/nucl_modA.py batch4/ D /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_D.tree nucl_outsD_batch4/ >> nuclD_batch4_codeml.log 2>&1
#mkdir nucl_outsD_batch5
#python /gpfs/scratch/mclim/nucl_modA.py batch5/ D /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_D.tree nucl_outsD_batch5/ >> nuclD_batch5_codeml.log 2>&1

mkdir nucl_outsD_redos
# copy these to Codeml_nucleargenes/convergence_redo_batch/
# /pos_redos/: ENSTGUP00000018138.phylip (batch1), ENSTGUP00000009446.phylip (batch2), ENSTGUP00000010247.phylip (batch2), ENSTGUP00000008478.phylip (batch3)
# /null_redos/:ENSTGUP00000010255.phylip (batch2), ENSTGUP00000006084.phylip (batch3)
python /gpfs/scratch/mclim/nucl_redo_modA.py D /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_D.tree nucl_outsD_redos/ >> nuclD_redos_codeml.log 2>&1