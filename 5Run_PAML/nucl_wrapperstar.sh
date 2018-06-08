#!/bin/bash
#PBS -l nodes=4:ppn=28,walltime=1:00:00
#PBS -N nuclstar_paml
#PBS -q short

# long, 24 hours, 4 nodes for the full batches
# less time needed for debugging and running redos

set echo

cd /gpfs/scratch/mclim/nucl_star
module load shared
module load anaconda/3 #note that this is using python 3
source activate biopython #to use Biopython

#mkdir nucl_outsstar_batch1
#python /gpfs/scratch/mclim/nucl_modA.py batch1/ star /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot.tree nucl_outsstar_batch1/ >> nuclstar_batch1_codeml.log 2>&1

#mkdir nucl_outsstar_batch2
#python /gpfs/scratch/mclim/nucl_modA.py batch2/ star /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot.tree nucl_outsstar_batch2/ >> nuclstar_batch2_codeml.log 2>&1
#mkdir nucl_outsstar_batch3
#python /gpfs/scratch/mclim/nucl_modA.py batch3/ star /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot.tree nucl_outsstar_batch3/ >> nuclstar_batch3_codeml.log 2>&1
#mkdir nucl_outsstar_batch4
#python /gpfs/scratch/mclim/nucl_modA.py batch4/ star /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot.tree nucl_outsstar_batch4/ >> nuclstar_batch4_codeml.log 2>&1
#mkdir nucl_outsstar_batch5
#python /gpfs/scratch/mclim/nucl_modA.py batch5/ star /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot.tree nucl_outsstar_batch5/ >> nuclstar_batch5_codeml.log 2>&1

#mkdir nucl_outsstar_redos
python /gpfs/scratch/mclim/nucl_redo_modA.py star /gpfs/scratch/mclim/codeml_trees/transcriptometree_noroot.tree nucl_outsstar_redos/ >> nuclstar_redos_codeml.log 2>&1
