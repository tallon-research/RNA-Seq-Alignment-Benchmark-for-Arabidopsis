# Invoke as python3 read_stats input.bam input_1.fasta input_2.fasta
# Second read file is optional.  If 2 arguments are present, assumes single-end.

import sys
import os
readlib = {}
scores = {}

# Function to test if alignment position is inside the gene transcript region.

def inbetween(num1,num2,num3):
    bignum = max(num2,num3)
    lilnum = min(num2,num3)
    if (num1 > lilnum) and (num1 < bignum):
        return True
    else:
        return False

# Input Arguments

bamfile = sys.argv[1]
readfile = sys.argv[2]

# Convert BAM to temporary SAM for quick reading.

command = 'samtools view -h ' + bamfile + ' > temp.sam'
os.system(command)

# Unzip Read File(s)

os.system('gunzip -f ' + sys.argv[2])
zreadfile = readfile.split('.gz')[0]

# Generate readfile information.

readset = readfile.split('errormodel')[1]
readset = 'errormodel' + readset.split('/')[0]
readfile = readfile.split('/')[-1]
readfilename = readfile.split('.')[0]

# Cycle through single-end reads (or paired set 1).

with open(zreadfile) as infile:
    for lines in infile:
        lines = lines.rstrip()
        if lines.startswith('>'):
            readname = lines.split('/')[0]
            readname = readname[1:]
            truechr = lines.split(':')[2]
            truecorg1 = lines.split(':')[3]
            truecorg2 = lines.split(':')[4]
            readlib[readname] = [truechr,int(truecorg1),int(truecorg2)]

# Check for third argument, if present, and append data from paired set 2.

if len(sys.argv) > 3:
    os.system('gunzip -f ' + sys.argv[3])
    zreadfile2 = (sys.argv[3]).split('.gz')[0]
    with open(zreadfile2) as infile:
        for lines in infile:
            lines = lines.rstrip()
            if lines.startswith('>'):
                readname = lines.split('/')[0]
                readname = readname[1:]
                if readname in readlib:
                    continue
                else:
                    truechr = lines.split(':')[2]
                    truecorg1 = lines.split(':')[3]
                    truecorg2 = lines.split(':')[4]
                    readlib[readname] = [truechr,int(truecorg1),int(truecorg2)]
                    
# Read through SAM file and generate statistics grouped by bitwise SAM flag.
                    
with open('temp.sam') as infile:
    for lines in infile:
        if lines.startswith('read'):
            lines = lines.rstrip()
            values = lines.split('\t')
            readname = values[0].split('/')[0]
            bitflag = values[1]
            truechr = values[2]
            truechr = truechr.split(' ')[0] # BBMap does not truncate anything.
            firsali = int(values[3])
            if (inbetween(firsali,readlib[readname][1],readlib[readname][2]) == True) and (truechr == readlib[readname][0]):
                passcheck = 'Correct'
            else:
                passcheck = 'Incorrect'
            if bitflag not in scores:
                scores[bitflag] = {'Correct':0,'Incorrect':0}
            scores[bitflag][passcheck] += 1
            
# Print the output
print('############################################################')
print(sys.argv[1])
print('Total Reads: ' + str(len(readlib)))
print('bitflag\tAlignment\tNumber of Alignments')
for s in scores:
    for m in scores[s]:
        print('{} {} {}'.format(s,m,scores[s][m]))
print('############################################################')
# Remove the temp SAM file.

os.system('rm temp.sam')
os.system('gzip ' + zreadfile)

if (len(sys.argv) > 3):
    os.system('gzip ' + zreadfile2)
