import sys
import os
import gzip
import tqdm

dam1 = sys.argv[1]
dam2 = sys.argv[2]
sam = sys.argv[3]
corctch = sys.argv[4] #sub will change SAM processing to avoid subread's corrupt output, bbm will convert = to M to match other aligners.
outfile = sys.argv[5]

alignermatch = {}
alignermiss = {}




def samProcessor(samorbamfile):
	filebase = samorbamfile.split('.')[0]
	newout = filebase + '.filtersam'
	if corctch != 'sub':
		os.system('samtools view -F 260 ' + samorbamfile + ' > ' + newout)
	else:
		os.system('cp ' + samorbamfile + ' ' + newout)
	return newout


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

def cigAdjust(actualstring,alignstring,actualstart,alignstart):
    if actualstring == '':
        return actualstring,alignstring
    diff = abs(actualstart - alignstart)
    if actualstart < alignstart:
        alignstring = ('*' * diff) + alignstring
    elif actualstart > alignstart:
        actualstring = ('*' * diff) + actualstring
    return actualstring,alignstring

def stringoverlap(string1,string2):
    # string1 is always actual.
	if len(string1) < len(string2):
		string1 = string1.ljust(len(string2),'*')
	else:
		string2 = string2.ljust(len(string1),'*')
	if corctch == 'bbm':
		string2 = string2.replace('=','M')
	for x in range(max(len(string1),len(string2))):
		if string2[x] != string1[x]:
			if string2[x] not in alignermiss:
				alignermiss[string2[x]] = {}
			if string1[x] not in alignermiss[string2[x]]:
				alignermiss[string2[x]][string1[x]] = 0
			alignermiss[string2[x]][string1[x]] += 1
		else:
			if string2[x] == '*':
				continue
			if string2[x] not in alignermatch:
				alignermatch[string2[x]] = 0
			alignermatch[string2[x]] += 1
	return

def errorpad(damstring):
	for x in range(len(damstring)):
		if damstring[x] not in alignermiss:
			alignermiss[damstring[x]] = {}
		if '*' not in alignermiss[damstring[x]]:
			alignermiss[damstring[x]]['*'] = 0
		alignermiss[damstring[x]]['*'] += 1
	return

def damload(damfile):
    damlist = {}
    with open(damfile,'rt') as infile:
        for lines in infile:
            lines = lines.rstrip()
            values = lines.split('\t')
            cigstring = cigarUnroller(values[3])
            mcount = cigstring.count('M')
            mcount += cigstring.count('=')
            ncount = cigstring.count('N')
            xcount = cigstring.count('X')
            dcount = cigstring.count('D')
            icount = cigstring.count('I')
            damlist[values[0]] = [values[1],int(values[2]),values[3],values[4],mcount,ncount,xcount,dcount,icount]
    return damlist

def samload(samfile):
    sam1list = {}
    sam2list = {}
    with open(samfile, encoding="utf8", errors='ignore') as infile:
        for lines in infile:
            lines = lines.rstrip()
            if lines.startswith('@'):
                continue
            else:
                values = lines.split('\t')
                if values[0] in sam1list:
                    sam2list[values[0]] = [values[2].split(' ')[0],int(values[3]),values[5]]
                else:
                    sam1list[values[0]] = [values[2].split(' ')[0],int(values[3]),values[5]]
    return sam1list,sam2list

def returnlength(cigar):
    if cigar == '':
        return 0
    else:
        numz = ''
        leng = 0
        for c in cigar:
            if c.isdigit():
                numz += c
            else:
                numz += ','
        numz = numz.split(',')
        for n in numz:
            try:
                leng += int(n)
            except:
                pass
        return leng

# Load in dam/sam files, store into lists of lists...

print('Loading DAM files')

dam1dict = damload(dam1)
dam2dict = damload(dam2)

print('Extracting primary alignments from SAM/BAM file')

filtsam = samProcessor(sam)

print('Loading primary alignments')
sam1dict,sam2dict = samload(filtsam)

# Cycle through, compare for overlap.  Use relative position to ensure proper matchup between dam and sam entries.
# Use chromosome check to bypass unecessary overlap for completely wrong alignments (wrong chromo).

uniquekeys = set(list(dam1dict.keys()) + list(dam2dict.keys()))

total = 0
totalm = 0
totaln = 0
totalx = 0
totald = 0
totali = 0

print('Generating Comparison Data')

for u in tqdm.tqdm(uniquekeys):
    # Handle possible missing sam entry.  Important for tophat2 split output (accepted_hits).
    if u not in sam1dict:
        sam1dict[u] = ['-99',0,'']
    if u not in sam2dict:
        sam2dict[u] = ['-99',1,'']
    # Handle possible missing dam entry.  Important for unpaired reads from polyester.
    if u not in dam1dict:
        dam1dict[u] = ['99',0,'']
    if u not in dam2dict:
        dam2dict[u] = ['99',1,'']
    # Relative positioning...
    if sam1dict[u][1] < sam2dict[u][1]:
        lessersam = sam1dict[u]
        greatersam = sam2dict[u]
    else:
        greatersam = sam1dict[u]
        lessersam = sam2dict[u]
    if dam1dict[u][1] < dam2dict[u][1]:
        lesserdam = dam1dict[u]
        greaterdam = dam2dict[u]
    else:
        greaterdam = dam1dict[u]
        lesserdam = dam2dict[u]
    # Adjust total bases.
    total += int(greaterdam[3])
    total += int(lesserdam[3])
    totalm += (greaterdam[4] + lesserdam[4])
    totaln += (greaterdam[5] + lesserdam[5])
    totalx += (greaterdam[6] + lesserdam[6])
    totald += (greaterdam[7] + lesserdam[7])
    totali += (greaterdam[8] + lesserdam[8])
    # Chromo check.
    if greaterdam[0] != greatersam[0]:
        errorpad(cigarUnroller(greaterdam[2]))
    else:
        greaterdamstring,greatersamstring = cigAdjust(cigarUnroller(greaterdam[2]),cigarUnroller(greatersam[-1]),greaterdam[1],greatersam[1])
        stringoverlap(greaterdamstring,greatersamstring)
    if lesserdam[0] != lessersam[0]:
        errorpad(cigarUnroller(lesserdam[2]))
    else:
        lesserdamstring,lessersamstring = cigAdjust(cigarUnroller(lesserdam[2]),cigarUnroller(lessersam[-1]),lesserdam[1],lessersam[1])
        stringoverlap(lesserdamstring,lessersamstring)
    

of = open(outfile,'w')
of.write('Total Read Bases\t' + str(total) + '\n')
of.write('Total Junction Bases\t' + str(totaln) + '\n')
of.write('Total Error Bases\t' + str(totalx) + '\n')
of.write('Total Deletion Bases\t' + str(totald) + '\n')
of.write('Total Insertion Bases\t' + str(totali) + '\n')

of.write('Actual=Prediction\tCount\n')
for a in alignermatch:
	of.write(a + '\t' + str(alignermatch[a]) + '\n')
of.write('Prediction\tActual\tCount\n')
for a in alignermiss:
	for a2 in alignermiss[a]:
		of.write(a + '\t' + a2 + '\t' + str(alignermiss[a][a2]) + '\n')

of.close()

# cleanup

os.system('rm ' + filtsam)
