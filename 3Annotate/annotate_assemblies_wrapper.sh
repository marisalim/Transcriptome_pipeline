#!/bin/bash
#PBS -l nodes=1:ppn=15
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -q batch
#PBS -m e
#PBS -M marisa.lim@stonybrook.edu
#PBS -N annotate_assemblies
#
set echo
module load blast/2.2.31-all
module load exonerate/2.2.0
module load cd-hit/4.6.4-all
module load perl/5.18.4-threads
#
# Take Trinity assemblies and annotate with zebra finch
#
perl /crucible/bi4iflp/mlim/Annotation/Annotate_assemblies_singhal.pl -a /crucible/bi4iflp/mlim/Trinity_assembly_files -b /crucible/bi4iflp/mlim/Annotation/Taeniopygia_guttata.taeGut3.2.4.pep.all.fa -c 1e-10 -e 15 -f /crucible/bi4iflp/mlim/Annotation/Taeniopygia_guttata_gene_name.txt
