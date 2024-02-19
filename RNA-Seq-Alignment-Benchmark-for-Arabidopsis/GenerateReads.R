library(polyester)
library(Biostrings)

# Read the fasta file into a Biostrings dataset
fasta = readDNAStringSet('/home/tallon/RNAseq_benchmark/Arabidopsis_genome/Arabidopsis_thaliana.TAIR10.cdna.all.fa')

# Calculate reads per transcript with normalization for transcript length and desired coverage
coverage = 2
read_length = 100
readspertx = round(coverage * width(fasta) / read_length)
readspertx = readspertx + 1  # Add 1 to avoid 0 count reads which can crash polyester

# Determine fold-change between conditions (setting all to 1 for simplicity)
# Using the length of the fasta file to know how many transcripts we have
fc = sample(c(1), size=length(fasta), replace=TRUE)

# Run the simulation experiment
simulate_experiment(
  dna = '/home/tallon/RNAseq_benchmark/Arabidopsis_genome/Arabidopsis_thaliana.TAIR10.cdna.all.fa',
  reads_per_transcript = readspertx,
  num_reps = c(3, 3),
  fold_changes = fc,
  readlen = 100,
  error_rate = 0.005,
  error_model = 'illumina5',
  paired = TRUE,
  distr = 'normal',
  strand_specific = TRUE,
  outdir = "/home/tallon/RNAseq_benchmark/100_005_Illumina5_Stranded_Paired"
)
