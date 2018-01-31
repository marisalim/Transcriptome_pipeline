#!/bin/bash
#PBS -l nodes=4:ppn=28,walltime=2:00:00
#PBS -N nuclE_paml_redos
#PBS -q short

# long, 24 hours, 4 nodes for the full batches
# less time needed for debugging and running redos

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
#mkdir nucl_outsE_batch4
#python /gpfs/scratch/mclim/nucl_modA.py batch4/ E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_batch4/ >> nuclE_batch4_codeml.log 2>&1
#mkdir nucl_outsE_batch5
#python /gpfs/scratch/mclim/nucl_modA.py batch5/ E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_batch5/ >> nuclE_batch5_codeml.log 2>&1

mkdir nucl_outsE_redos
# copy these to Codeml_nucleargenes/convergence_redo_batch/
# /pos_redos/: ENSTGUP00000009446.phylip (batch2), ENSTGUP00000005605.phylip (batch4), ENSTGUP00000002255.phylip (batch4)
# /null_redos/: ENSTGUP00000013071.phylip (batch1), ENSTGUP00000013188.phylip (batch1), ENSTGUP00000002255.phylip (batch4), ENSTGUP00000000654.phylip (batch5)
python /gpfs/scratch/mclim/nucl_redo_modA.py E /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot_E.tree nucl_outsE_redos/ >> nuclE_redos_codeml.log 2>&1