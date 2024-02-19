# Use as python3 fa2fq.py read.fasta.gz
# Creates read.fastq.gz

import sys
import gzip
import os

fasta = {}

# Read in the file from arg1

with gzip.open(sys.argv[1],mode="rt") as infile:
    for lines in infile:
        lines = lines.rstrip()
        if lines.startswith('>'):
            lines = lines.replace('>','@')
            fasta[lines] = ''
            current = lines
        else:
            fasta[current] += lines

# Export Fastq with I quality

fastq = open(sys.argv[1].split('.fasta')[0] + '.fastq','w')
for f in fasta:
    fastq.write(f + '\n')
    fastq.write(fasta[f] + '\n')
    fastq.write('+\n')
    fastq.write("I"*len(fasta[f]) + '\n')

fastq.close()

os.system("gzip *fastq")
