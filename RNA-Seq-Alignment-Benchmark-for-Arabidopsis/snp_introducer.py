# This program will first pull in reads from Polyester (Stranded, with Error) and compare to the source to find the variations.  Keeping a record.
# Next, variants will be introduced at user-defined percentages, also keeping record.
# The ideal CIGAR string based on the source will be determined.  Then, modifications via the variant and error profile will be used to modify the CIGAR string to its actual.

# Run with python3 snp_introducer.py <input.fasta> <output prefix>

import random
import sys

intronLoc = {}
intronLen = {}
readorient = {}

outfilenameprefix = sys.argv[2]
of = open(outfilenameprefix + '.dam','w')
of2 = open(outfilenameprefix + '.fasta','w')


add_perc = int(sys.argv[3])

with open('Arabidopsis_thaliana.TAIR10.46.introns2') as infile:
    for lines in infile:
        lines = lines.rstrip()
        values = lines.split('\t')
        tscript = values[2]
        if tscript not in intronLoc:
            startloc = int(values[5])
            endloc = int(values[6])
            #lentally = endloc - startloc
            lentally = endloc - startloc + 1
            intronLoc[tscript] = [startloc]
            intronLen[tscript] = [lentally]
        else:
            nextstartloc = int(values[5])
            nextendloc = int(values[6])
            #nextlen = nextendloc - nextstartloc
            nextlen = nextendloc - nextstartloc + 1
            intronLoc[tscript].append(nextstartloc - lentally)
            lentally += nextlen
            intronLen[tscript].append(nextlen)
            
# Get transcript lengths

transcriptlen = {}

with open('Arabidopsis_thaliana.TAIR10.cdna.all.fa') as infile:
    for lines in infile:
        if lines.startswith('>'):
            curtran = (lines.split(' ')[0])[1:]
            transcriptlen[curtran] = 0
        else:
            lines = lines.rstrip()
            transcriptlen[curtran] += len(lines)

def cigarUnroller(cigarString):
    try:
        curNum = ''
        curLet = ''
        cigar = ''
        for char in cigarString:
            if char.isdigit() == True:
                curNum += char
            else:
                curLet = char
                cigar += (curLet * int(curNum))
                curNum = ''
        return cigar
    except:
        return ''

def revComp(seq):
    seq = seq.replace('A','1')
    seq = seq.replace('G','2')
    seq = seq.replace('T','A')
    seq = seq.replace('C','G')
    seq = seq.replace('1','T')
    seq = seq.replace('2','C')
    seq = seq[::-1]
    return seq

# Store as variantC[transcript] = [type,start,end,replacement]

def variant_catalogue(varfile):
    variantC = {}
    with open(varfile) as infile:
        for lines in infile:
            values = lines.split('\t')
            transcript = values[0]
            vartype = values[6]
            start = int(values[2])
            end = int(values[3])
            replacement = (values[-1].split('Variant_seq=')[1]).split(';')[0]
            replacement = random.choice(replacement.split(','))
            replacement = replacement.replace('-','')
            if transcript not in variantC:
                variantC[transcript] = [[start,end,replacement]]
            else:
                variantC[transcript].append([start,end,replacement])
    return variantC

# Load in the exon sequence and fasta file to generate ideal.  Load in real.  Add errors from variantC file.  Generate ideal CIGAR.  Generate modified from differences.

def readloader(exonfile,readfile):
    readlib = {}
    if '_1.fa' in readfile:
        mate = 'mate1:'
    else:
        mate = 'mate2:'
    transcripts = {}
    with open(exonfile) as infile:
        for lines in infile:
            lines = lines.rstrip()
            name = lines.split('_')[0]
            if name not in transcripts:
                transcripts[name] = ''
            seq = lines.split('\t')[5]
            transcripts[name] += seq
    with open(readfile) as infile:
        for lines in infile:
            if lines.startswith('>'):
                tx = lines.split(' ')[0]
                tx = tx.split('/')[1]
                re = lines.split('/')[0][1:]
                try:
                    loc = lines.split(mate)[1]
                except:
                    continue
                loc = loc.split(';')[0]
                loc1,loc2 = [int(x) for x in loc.split('-')]
                read_think = transcripts[tx][loc1-1:loc2]
                if mate == 'mate2:':
                    read_think = revComp(read_think)
            else:
                read_is = lines.rstrip()
                readlib[re] = [read_think,read_is,loc1,loc2,tx]
    return readlib

def cigarRoller(readFile,cigLibrary):
    with open(readFile) as infile:
        for lines in infile:
            lines = lines.rstrip()
            readlocs = []
            readlens = []
            if lines.startswith('>'):
                orientation = (lines.split(' gene')[0]).split(':')[-1]
                if '_1.fa' in readFile:
                    mate = 'mate1:'
                elif '_2.fa' in readFile:
                    mate = 'mate2:'
                else:
                    print('Error in :' + lines)
                    continue
                if mate not in lines:
                    continue
                read = lines[1:].split('/')[0]
                if read not in readorient:
                    readorient[read] = orientation
                tscript = (lines[1:].split('/')[1]).split(' ')[0]
                locz = (lines.split(mate)[1]).split(';')[0]
                if orientation == '1':
                    start = int(locz.split('-')[0])
                    end = int(locz.split('-')[1])
                elif orientation == '-1':
                    start = int(locz.split('-')[0])
                    end = int(locz.split('-')[1])
                    fraglen = (end - start)
                    tlengo = transcriptlen[tscript]
                    start = tlengo - end + 1
                    end = start + fraglen
                schr = lines.split(':')[2]
                sstart = int(lines.split(':')[3])
                sstart = sstart + start - 1
                # Cycle through intron locations.
                if tscript in intronLoc:
                    for x in range(len(intronLoc[tscript])):
                        ilocation = intronLoc[tscript][x]
                        if (ilocation > start) and (ilocation < end):
                            readlocs.append(ilocation - start + 1)
                            readlens.append(intronLen[tscript][x])
                        if (ilocation <= start):
                            sstart += intronLen[tscript][x]
                if len(readlocs) == 0:
                    #fullinfo = [schr,sstart,str(end-start+1) + 'M']
                    fullinfo = [schr,sstart,str(end-start + 1) + 'M',orientation]
                    cigLibrary[read] += fullinfo
                    #print(read + '\t' + schr + '\t' + str(sstart) + '\t' + tscript + '\t' + str(end-start+1) + 'M')
                else:
                    cigar = ''
                    runtot = 0
                    for r in range(len(readlocs)):
                        cigar += str(readlocs[r]-runtot) + 'M'
                        runtot += (readlocs[r]-runtot)
                        cigar += str(readlens[r]) + 'N'
                    cigar += str(end-start+1-runtot) + 'M'
                    fullinfo = [schr,sstart,cigar,orientation]
                    cigLibrary[read] += fullinfo
    return cigLibrary

def cuba(cigstring):
    curcount = 1
    curletter = cigstring[0]
    cigar = ''
    for c in cigstring[1:]:
        if c == curletter:
            curcount += 1
        else:
            cigar += str(curcount) + curletter
            curletter = c
            curcount = 1
    cigar += str(curcount) + curletter
    return cigar
            

varC = (variant_catalogue('Exonic_Errors.tsv'))
dicto = (readloader('Culled_Exon_Var_List_plus1.tsv',sys.argv[1]))
dicto = (cigarRoller(sys.argv[1],dicto))

for d in dicto:
    try:
        if dicto[d][-1] == '-1':
            dicto[d][0] = revComp(dicto[d][0])
            dicto[d][1] = revComp(dicto[d][1])
        if len(dicto[d][0]) != len(dicto[d][1]):
            print(dicto[d])
            continue
        modifications = []
        cigar = cigarUnroller(dicto[d][-2])
        actual = ''
        theoretical = ''
        start = int(dicto[d][6])
        i = 0
        for c in cigar:
            if c == 'M':
                theoretical += dicto[d][0][i]
                actual += dicto[d][1][i]
                i += 1
            else:
                theoretical += '9'
                actual += '9'
        i = 0
        cigar = list(cigar)
        actual = list(actual)
        for nuc in theoretical:
            if nuc != actual[i]:
                cigar[i] = 'X'
            i += 1
        end = start + len(cigar)
        if dicto[d][-5] in varC:
            for v in varC[dicto[d][-5]]:
                if (v[0] >= start) and (v[1] <= end):
                    if random.randint(0,100) < add_perc:
                        modifications.append(v)
        else:
            modifications = []
        for m in reversed(modifications):
            adjstart = m[0] - start
            adjend = m[1] - start + 1
            old = actual[adjstart:adjend]
            new = list(m[2])
            if len(old) == len(new):
                cigar[adjstart:adjend] = ['X'] * len(new)
                actual[adjstart:adjend] = new
            elif len(old) > len(new):
                newcig = []
                for o in old:
                    new.append('*')
                for n in range(len(old)):
                    if new[n] == old[n]:
                        newcig.append('M')
                    elif new[n] == '*':
                        newcig.append('D')
                    else:
                        newcig.append('X')
                cigar[adjstart:adjend] = newcig
                actual[adjstart:adjend] = list(m[2])
            else:
                newcig = []
                for n in new:
                    old.append('*')
                for o in range(len(new)):
                    if old[o] == new[o]:
                        newcig.append('M')
                    elif old[o] == '*':
                        newcig.append('I')
                    else:
                        newcig.append('X')
                cigar[adjstart:adjend] = newcig
                actual[adjstart:adjend] = list(m[2])
        cigar = ''.join(cigar)
        cigar = cuba(cigar)
        actual = (''.join(actual))
        actual = actual.replace('9','')
        of2.write('>' + d + '\n')
        if dicto[d][-1] == '1':
            of2.write(actual + '\n')
            of.write(d + '\t' + dicto[d][-4] + '\t' + str(dicto[d][6]) + '\t' + cigar + '\t' + str(len(actual)) + '\t' + dicto[d][-5] + '\n')
        else:
            of2.write(revComp(actual) + '\n')
            of.write(d + '\t' + dicto[d][-4] + '\t' + str(dicto[d][6]) + '\t' + cigar + '\t' + str(len(actual)) + '\t' + dicto[d][-5] + '\n')
        #print(d + '\t' + actual + '\t' + str(dicto[d][6]) + '\t' + cigar + '\t' + dicto[d][-5])
    except:
        # Clause for polyester busted reads (not sure why, bug in their code).
        continue
        
        
    
    
