 # Open GVF, extract percentage of variants.
# Open GTF, build transcripts and add in variants.
gtf_file = 'Arabidopsis_thaliana.TAIR10.46.gtf'
genome_file = 'Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'

genome = {}
exons = {}

def revComp(seq):
    seq = seq.replace('A','1')
    seq = seq.replace('G','2')
    seq = seq.replace('T','A')
    seq = seq.replace('C','G')
    seq = seq.replace('1','T')
    seq = seq.replace('2','C')
    seq = seq[::-1]
    return seq

with open(genome_file) as infile:
    for lines in infile:
        lines = lines.rstrip()
        if lines.startswith('>'):
            chromo = lines[1:].split(' ')[0]
            genome[chromo] = ''
        else:
            genome[chromo] += lines

with open(gtf_file) as infile:
    for lines in infile:
        if lines.startswith('#'):
            continue
        lines = lines.rstrip()
        values = lines.split('\t')
        if values[2] == 'exon':
            transcript = (lines.split('transcript_id "')[1]).split('"')[0]
            exon_num = (lines.split('exon_number "')[1]).split('"')[0]
            strand = values[6]
            chromo = values[0]
            start = int(values[3])
            end = int(values[4])
            if strand == '+':
                print(transcript + '_' + exon_num + '\t' + chromo + '\t' + strand + '\t' + str(start) + '\t' + str(end) + '\t' + genome[chromo][start-1:end])
            elif strand == '-':
                print(transcript + '_' + exon_num + '\t' + chromo + '\t' + strand + '\t' + str(start) + '\t' + str(end) + '\t' + revComp(genome[chromo][start-1:end]))
