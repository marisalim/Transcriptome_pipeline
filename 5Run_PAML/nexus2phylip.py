#http://biopython.org/wiki/AlignIO

#!/usr/bin/env python

# Script converts nexus file to phylip file. From Mark Phuong.

# Import modules
from Bio import SeqIO
import dendropy
import os
import sys
from collections import defaultdict
from pprint import pprint
import argparse
import multiprocessing
from Bio.Alphabet import generic_dna, Gapped
from Bio import AlignIO

# input alignment
alignment = AlignIO.read(open(sys.argv[1]), 'nexus',alphabet=Gapped(generic_dna))

# output name to write
output = open(sys.argv[2], 'w')

# write converted output
AlignIO.write(alignment, output, "phylip-sequential")

