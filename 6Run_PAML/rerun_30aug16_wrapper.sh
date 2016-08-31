#!/bin/bash
#SBATCH -N 1
#SBATCH -p LM
#SBATCH --mem=128GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marisa.lim@stonybrook.edu 

#echo commands to stdout
set -x

module load python/2.7.11_gcc
module load paml/4.9a

# run codeml script
# reruns for Acast
cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Acast
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Acast_nullmod.py

cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Acast/Acast_pos
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Acast_modA.py

# reruns for Cviol
cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Cviol
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Cviol_nullmod.py

cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Cviol/Cviol_pos
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Cviol_modA.py

# reruns for Ccoru
cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Ccoru
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Ccoru_nullmod.py

cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Ccoru/Ccoru_pos
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Ccoru_modA.py

# reruns for Mphoe
#cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Mphoe
#python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Mphoe_nullmod.py

cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Mphoe/Mphoe_pos
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Mphoe_modA.py

# reruns for Pgiga
cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Pgiga
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Pgiga_nullmod.py

cd /pylon1/bi4iflp/mlim/Run_codeml_nudna/nuDNA_codemlout_Pgiga/Pgiga_pos
python /pylon1/bi4iflp/mlim/Run_codeml_nudna/nu_Pgiga_modA.py