#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=2:00:00
#PBS -N test_etetoolkit
#PBS -q short

module load shared
module load anaconda/3
source activate etetoolkit

#-t tree file
#--alg input fasta file
#-o output file directory name with all outputs from PAML and control file
#--model lists codeml models
#--tests models used in LRT, first model is the null,second model is alt
#--clear_tree removes all tags on tree before doing next set
#--clear_all overwrites any data in output directory
#--mark tells etetoolkit which branch to tag on tree

# trees
# c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz
# e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri
# c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart
# k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri
# c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga
# f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri
# c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls
# d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

cd /gpfs/scratch/mclim/etetoolkit_pamlsims/

#echo 'ENSTGUP00000008523_MRPL42'
#ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000008523_sub.fasta -o ENSTGUP00000008523_MRPL42 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000013977_CLPB'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000013977_sub.fasta -o ENSTGUP00000013977_CLPB --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

