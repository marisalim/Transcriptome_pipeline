Transcriptome_pipeline
-
Author: Marisa Lim

Contact: marisalim@stonybrook.edu

This repository contains code for assembling and analyzing transcriptome data.

---
List of scripts and brief explanation of function. Miscellaneous scripts (for testing, old, etc.) are in the script_sandbox/ directory.
Code contributions from Sonal Singhal and Mark Phuong.

Python scripts were written for Python2.7. All scripts were run on Mac or Linux operating systems.

### 1Data_filtering_cleaning: 
Script  | Step
--- | --- 
1passfilterandmodifyname.py | remove reads that failed Illumina filter, edit sequence ids
2conc.py | concatenate read (R1 and R2) files for each sample
3runfastqc.py | run FastQC on reads pre-cleanup
4cleanreads.py | for 1 sample at a time, remove adapters and low quality reads, trim ends
4cleanreadsmulti.py | for all dataset samples, remove adapters and low quality reads, trim ends
5runfastqc_postclean.py | run FastQC on reads post-cleanup

### 2Assembly:
Script | Step 
--- | ---
trinity211[sample]_wrapper.sh | batch scripts for Trinity v2.1.1 assembly
2assemblystats.r | create histograms for assembly contigs

### 3Annotate:
Script | Step 
--- | ---
Annotate_assemblies_singhal.pl | define CDS regions and annotate
annotate_assemblies_wrapper.sh | wrapper script that runs Annotate_assemblies_singhal.pl on server
Get_CDS.py | extract CDS regions for alignment
Get_CDS_wrapper.sh | wrapper script that runs Get_CDS.py on server
rename_cdsfiles.py | format species names in fasta files

### 4Alignments:
Script | Step 
--- | ---
mafft_align.py | align annotated CDS with Mafft
mafft_align_wrapper.sh | wrapper script that runs mafft_align.py on server

### 5Run_PAML:
Script | Step | Type 
--- | --- | --- |
rename_sp.py | edit fasta file header names |
1.0convertnexusphylip.py | convert fasta > nexus > phylip for PAML |
fasta2Nexus.py  | conversion script used by 1.0convergnexusphylip.py for fasta to nexus |
nexus2phylip.py | conversion script used by 1.0convergnexusphylip.py for nexus to phylip |
1.1checkcodonlength.py | make sure sequence lengths are multiples of 3, otherwise PAML will give error |
nu_[species abbr]_modA.py | run PAML model A for nuclear genes; use nu_[species abbr]_wrapper.sh | single-branch foreground analyses
nu_[species abbr]_nullmod.py | run PAML null model for nuclear genes; use nu_[species abbr]_null_wrapper.sh | single-branch foreground analyses
mt_[species abbr]_modA.py | run PAML for mitochondrial genes; use mt_[species abbr]_wrapper.sh | single-branch foreground analyses
4aCalc_LRT_highaltsp.r | calculate LRT from codeml results, calculate corrected p-value, analyze positively selected genes | single-branch foreground analyses
nucl_modA.py | run PAML for nuclear genes; use nucl_wrapper[B,D,E,star].sh | multi-branch foreground analyses
mt_modA.py | run PAML for mitochondrial genes; use mt_wrapper[2,3,B,D,E, star].sh | multi-branch foreground analyses
4bCalc_LRT_Jan2018reanalysis.r | calculate LRT from codeml results, analyze PSGs | multi-branch foreground analyses

---


License: This code is available under a BSD 2-Clause License.

Copyright (c) 2015. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
in the documentation and/or other materials provided with the distribution. 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
