List of scripts and brief explanation of function. Miscellaneous scripts (for testing, old, etc.) are in the script_sandbox/ directory.
Code contributions from Sonal Singhal and Mark Phuong.

# 1Data_filtering_cleaning: 
1passfilterandmodifyname.py (remove reads that failed Illumina filter, edit sequence ids)
2conc.py (concatenate read files for each sample)
3runfastqc.py (run FastQC pre-cleanup)
4cleanreads.py or 4cleanreadsmulti.py (remove adapters and low quality reads, trim ends)
5runfastqc_postclean.py (run FastQC post-cleanup)

# 2Assembly:
trinity211[sample]_wrapper.sh (batch scripts for Trinity v2.1.1 assembly)
2assemblystats.r (create histograms for assembly contigs)

# 3Annotate:
Annotate_assemblies_singhal.pl (define CDS regions and annotate; use annotate_assemblies_wrapper.sh) 
Get_CDS.py (extract CDS regions for alignment; use Get_CDS_wrapper.sh)
rename_cdsfiles.py (format species names in fasta files)

# 4Alignments:
mafft_align.py (align annotated CDS with Mafft; use mafft_align_wrapper.sh)

# 5Run_PAML:
rename_sp.py (edit fasta file names)
1.0convertnexusphylip.py (convert fasta > nexus > phylip for PAML; use fasta2Nexus.py and nexus2phylip.py)
1.1checkcodonlength.py (make sure sequence lengths are multiples of 3 - otherwise PAML will give error)
  # single-branch foreground analyses:
  nu_[species abbr]_modA.py (run PAML model A for nuclear genes; use nu_[species abbr]_wrapper.sh)
  nu_[species abbr]_nullmod.py (run PAML null model for nuclear genes; use nu_[species abbr]_null_wrapper.sh)
  mt_[species abbr]_modA.py (run PAML for mitochondrial genes; use mt_[species abbr]_wrapper.sh) 
  4aCalc_LRT_highaltsp.r (calculate LRT from codeml results, calculate corrected p-value, analyze positively selected genes)
  
  # multi-branch foreground analyses:
  nucl_modA.py (run PAML for nuclear genes; use nucl_wrapper[B,D,E,star].sh)
  mt_modA.py (run PAML for mitochondrial genes; use mt_wrapper[2,3,B,D,E, star].sh)    
  4bCalc_LRT_Jan2018reanalysis.r (calculate LRT from codeml results, analyze PSGs)







