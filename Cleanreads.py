#!/usr/bin/env python

# Code to clean reads: remove low quality reads/PCR dups/optical dups (Trimmomatic), trim adapters (Trimmomatic),
# merge overlapping PE reads (FLASH, COPE), filter for contamination (use Bowtie2 to align PE reads to possible contamination source genomes)

# set path


# import modules
import os, sys, multiprocessing

