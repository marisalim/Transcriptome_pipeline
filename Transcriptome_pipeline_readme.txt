List of scripts and brief explanation of function. Code contributions from Sonal Singhal and Mark Phuong are noted.
- Python scripts were written for Python2.7. All scripts were run on Mac or Linux operating systems.
- As currently written, for any scripts that run 1 sample at a time, you need to manually change directory paths for each sample. Have included 1 script as an example.

### 1Data_filtering_cleaning: 
1passfilterandmodifyname.py - remove reads that failed Illumina filter, edit sequence ids
2conc.py - concatenate read (R1 and R2) files for each sample
3runfastqc.py - run FastQC on reads pre-cleanup
4cleanreads.py - for 1 sample at a time, remove adapters and low quality reads, trim ends
4cleanreadsmulti.py - for all dataset samples, remove adapters and low quality reads, trim ends
5runfastqc_postclean.py - run FastQC on reads post-cleanup

### 2Assembly:
trinity211_wrapper.sh - scripts for Trinity v2.1.1 assembly on server, 1 sample at a time
2assemblystats.r - create histograms for assembly contigs

### 3Annotate:
Annotate_assemblies_singhal.pl - define CDS regions and annotate; script from Sonal Singhal
annotate_assemblies_wrapper.sh - wrapper script that runs Annotate_assemblies_singhal.pl on server
Get_CDS.py - extract CDS regions for alignment
Get_CDS_wrapper.sh - wrapper script that runs Get_CDS.py on server
rename_cdsfiles.py - format species names in fasta files

### 4Alignments:
mafft_align.py - align annotated CDS with Mafft
mafft_align_wrapper.sh - wrapper script that runs mafft_align.py on server

### 5Run_PAML:
rename_sp.py - edit fasta file header names 
1.0convertnexusphylip.py - convert fasta > nexus > phylip for PAML 
fasta2Nexus.py - conversion script used by 1.0convergnexusphylip.py for fasta to nexus; script from Mark Phuong
nexus2phylip.py - conversion script used by 1.0convergnexusphylip.py for nexus to phylip; script from Mark Phuong
1.1checkcodonlength.py - make sure sequence lengths are multiples of 3, otherwise PAML will give error 
nu_[species abbr]_modA.py - run PAML model A for nuclear genes, 1 species at a time; use nu_[species abbr]_wrapper.sh (single-branch foreground)
nu_[species abbr]_nullmod.py - run PAML null model for nuclear genes, 1 species at a time; use nu_[species abbr]_null_wrapper.sh (single-branch foreground)
mt_[species abbr]_modA.py - run PAML for mitochondrial genes, 1 species at a time; use mt_[species abbr]_wrapper.sh (single-branch foreground)
nucl_modA.py - run PAML for nuclear genes; use nucl_wrapper [B,D,E,star].sh (multi-branch foreground)
mt_modA.py - run PAML for mitochondrial genes; use mt_wrapper.sh (multi-branch foreground)

### 6Data_analysis:
1aCalc_LRT_highaltsp.r - calculate LRT from codeml results, calculate corrected p-value, analyze positively selected genes (single-branch foreground)
1bCalc_LRT_Jan2018reanalysis.r - calculate LRT from codeml results, analyze PSGs (multi-branch foreground)
1cGOterm_comparison_treeBstarDE.r - compare gene ontology (GO) terms for PAML results from different branch tagging scenarios (multi-branch foreground)
