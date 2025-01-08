#adapted from https://stackoverflow.com/questions/48240579/how-to-speed-up-fasta-subsampling-program-for-python
#!/usr/bin/env python3
import random
import sys

# Import arguments
infile = open(sys.argv[1], "r")
outfile = open(sys.argv[2], "w")

# Define list and dictionary
fNames = []
dicfastq_seq = {}
dicfastq_qual = {}

# Extract fastq file into the two lists
for line in infile:
    if line.startswith("@"):
        fNames.append(line.rstrip())
        Id = line.rstrip()
        line_group_count = 1
    elif line_group_count == 1:
        dicfastq_seq[Id] = line.rstrip()
        line_group_count += 1
    elif line_group_count == 2:
        line_group_count += 1
    elif line_group_count == 3:
        dicfastq_qual[Id] = line.rstrip()

# Print total number of sequences in the original file
print("There are "+str(len(fNames))+" in the input file")
# How many random sequences do you want?
num = float(input("What fraction of reads do you want to rarefy to?\n"))

# Create subsamples
subsample_num = float(len(fNames))*num
subsample = []
subsample = random.sample(fNames, int(subsample_num))

# Take random items out of the list for the total number of samples required
for j in subsample:
    print(j, file = outfile)
    print(dicfastq_seq[j], file = outfile)
    print('+', file = outfile)
    print(dicfastq_qual[j], file = outfile)

infile.close()
outfile.close()
print("Done.")