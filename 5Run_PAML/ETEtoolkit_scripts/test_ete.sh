#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=10:00:00
#PBS -N test_etetoolkit
#PBS -q long

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

echo 'ENSTGUP00000002414_AGTRAP'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000002414_sub.fasta -o ENSTGUP00000002414_AGTRAP --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000004086_EPAS1'
ete3 evol -t transcriptometree_noroot_notag.nwk --alg Renamed_fastas/ENSTGUP00000004086_sub.fasta -o ENSTGUP00000004086_EPAS1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000005988_ATG9A'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000005988_sub.fasta -o ENSTGUP00000005988_ATG9A --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000013829_ATP1B1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000013829_sub.fasta -o ENSTGUP00000013829_ATP1B1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000001031_EEF2'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000001031_sub.fasta -o ENSTGUP00000001031_EEF2 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000018304_COX1'
ete3 evol -t transcriptometree_noroot_notag.nwk --alg Renamed_fastas/ENSTGUP00000018304_sub.fasta -o ENSTGUP00000018304_COX1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000018302_ND1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000018302_sub.fasta -o ENSTGUP00000018302_ND1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000007661_NDUFS8'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000007661_sub.fasta -o ENSTGUP00000007661_NDUFS8 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000006887_NDUFB10'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000006887_sub.fasta -o ENSTGUP00000006887_NDUFB10 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000013770_NDUFB4'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000013770_sub.fasta -o ENSTGUP00000013770_NDUFB4 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000000167_CACTIN'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000000167_sub.fasta -o ENSTGUP00000000167_CACTIN --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000000194_UQCRQ'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000000194_sub.fasta -o ENSTGUP00000000194_UQCRQ --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000001468_TIMM17A'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000001468_sub.fasta -o ENSTGUP00000001468_TIMM17A --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000001492_unknown'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000001492_sub.fasta -o ENSTGUP00000001492_unknown --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000001501_DNAJA1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000001501_sub.fasta -o ENSTGUP00000001501_DNAJA1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000001776_CCNI'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000001776_sub.fasta -o ENSTGUP00000001776_CCNI --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000002317_MDH1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000002317_sub.fasta -o ENSTGUP00000002317_MDH1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000003366_MRPS14'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000003366_sub.fasta -o ENSTGUP00000003366_MRPS14 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000003786_TSPAN3'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000003786_sub.fasta -o ENSTGUP00000003786_TSPAN3 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000004374_C1QBP'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000004374_sub.fasta -o ENSTGUP00000004374_C1QBP --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000004383_CCT3'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000004383_sub.fasta -o ENSTGUP00000004383_CCT3 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000005075_MEA1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000005075_sub.fasta -o ENSTGUP00000005075_MEA1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000006121_DNAJB2'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000006121_sub.fasta -o ENSTGUP00000006121_DNAJB2 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000006200_CSTL1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000006200_sub.fasta -o ENSTGUP00000006200_CSTL1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000007318_RTN4'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000007318_sub.fasta -o ENSTGUP00000007318_RTN4 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000008371_NOB1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000008371_sub.fasta -o ENSTGUP00000008371_NOB1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000008523_MRPL42'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000008523_sub.fasta -o ENSTGUP00000008523_MRPL42 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000008535_MFGE8'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000008535_sub.fasta -o ENSTGUP00000008535_MFGE8 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000008539_SDHA'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000008539_sub.fasta -o ENSTGUP00000008539_SDHA --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000009466_PDHB'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000009466_sub.fasta -o ENSTGUP00000009466_PDHB --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000009734_unknown'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000009734_sub.fasta -o ENSTGUP00000009734_unknown --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000010265_EEF1B2'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000010265_sub.fasta -o ENSTGUP00000010265_EEF1B2 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000011312_MRPS26'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000011312_sub.fasta -o ENSTGUP00000011312_MRPS26 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000011555_MLF1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000011555_sub.fasta -o ENSTGUP00000011555_MLF1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000011919_ATP6V1D'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000011919_sub.fasta -o ENSTGUP00000011919_ATP6V1D --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000011922_JPH1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000011922_sub.fasta -o ENSTGUP00000011922_JPH1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000013074_CDKN1B'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000013074_sub.fasta -o ENSTGUP00000013074_CDKN1B --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000013412_AAMDC'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000013412_sub.fasta -o ENSTGUP00000013412_AAMDC --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000013977_CLPB'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000013977_sub.fasta -o ENSTGUP00000013977_CLPB --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000014330_unknown'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000014330_sub.fasta -o ENSTGUP00000014330_unknown --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000014844_HADHB'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000014844_sub.fasta -o ENSTGUP00000014844_HADHB --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000015014_unknown'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000015014_sub.fasta -o ENSTGUP00000015014_unknown --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000017602_MGST3'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000017602_sub.fasta -o ENSTGUP00000017602_MGST3 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000018128_SEPP1'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000018128_sub.fasta -o ENSTGUP00000018128_SEPP1 --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri

echo 'ENSTGUP00000018150_FEM1A'
ete3 evol -t /gpfs/scratch/mclim/etetoolkit_pamlsims/transcriptometree_noroot_notag.nwk --alg /gpfs/scratch/mclim/etetoolkit_pamlsims/Renamed_fastas/ENSTGUP00000018150_sub.fasta -o ENSTGUP00000018150_FEM1A --model bsA bsA1 --test bsA1,bsA --clear_tree --mark c_Pmala=d_Acast=g_Cviol=k_Amela=l_Pgiga=h_Aamaz e_Ccoru=f_Ccoel=a_Phart=i_Mphoe=b_Cmuls=j_Aviri c_Pmala=e_Ccoru=d_Acast=f_Ccoel=g_Cviol=a_Phart k_Amela=i_Mphoe=l_Pgiga=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=d_Acast=k_Amela=i_Mphoe=l_Pgiga f_Ccoel=g_Cviol=a_Phart=b_Cmuls=h_Aamaz=j_Aviri c_Pmala=e_Ccoru=g_Cviol=a_Phart=l_Pgiga=b_Cmuls d_Acast=f_Ccoel=k_Amela=i_Mphoe=h_Aamaz=j_Aviri


