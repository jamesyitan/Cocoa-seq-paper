import os
import sys
import pickle
import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dir_path = os.path.dirname(os.path.realpath(__file__)) #get the path of the current directory
#input to terminal will be the name of the directory to process (should be the same as the one you demultiplexed)
dir_name = sys.argv[1]

#load pickle data files
with open(dir_path+'/'+dir_name+'/barcode_read_counts.pickle', 'rb') as barcode_data:
    barcode_read_counter = pickle.load(barcode_data)
with open(dir_path+'/'+dir_name+'/standard_umi_counts.pickle', 'rb') as std_umi_data:
    standard_umi_counter = pickle.load(std_umi_data)

##Remove barcodes that have mock community 16S amplicons in them, just selecting the "empty" droplets
#for each barcode, determine ratio of reads that are standards
frac_per_bc = dict()

clean_bc_list = list(standard_umi_counter.keys())
for bc in clean_bc_list:
    num_reads = barcode_read_counter[bc]
    num_stds = sum(standard_umi_counter[bc].values())
    frac_std_to_16S = num_stds/(num_reads) #dict for reads includes umi standards too
    frac_per_bc[bc] = frac_std_to_16S

for bc in clean_bc_list:
    if frac_per_bc[bc] < 0.90: #to select for droplets with mostly standards
        del standard_umi_counter[bc]

##UMI clean-up
bc_list = list(standard_umi_counter.keys())
start_num_bc = len(bc_list)

#remove barcode samples that have less than 20 unique UMIs and that have less than 1000 reads/barcode
for bc in bc_list:
    if len(standard_umi_counter[bc].keys()) < 20 or sum(standard_umi_counter[bc].values()) < 1000:
        del standard_umi_counter[bc]

clean_bc_list = list(standard_umi_counter.keys())
clean_num_bc = len(clean_bc_list)
print("Number of total non-noise barcode groups: %d" %(clean_num_bc))

#remove UMIs that are present at less than 2 occurrences
for bc in clean_bc_list:
    mean_umi_abun_per_bc = sum(standard_umi_counter[bc].values())/len(standard_umi_counter[bc])
    for umi in list(standard_umi_counter[bc].keys()):
        if standard_umi_counter[bc][umi] < 2:
            removed_umi = standard_umi_counter[bc].pop(umi)

##For a random subset of the barcodes
# Select a random subset of categories
subset_standard_umi_counter = random.sample(standard_umi_counter.keys(), 10)
# Extract the subset of the nested dictionary
subset_dict = {bc: standard_umi_counter[bc] for bc in subset_standard_umi_counter}

# Flatten the subset of the nested dictionary into a DataFrame
subset_standard_umi_counter = []
for bc, umi in subset_dict.items():
    for umi, count in subset_dict[bc].items():
        subset_standard_umi_counter.append({'barcode': bc, 'umi': umi, 'count': count})

df = pd.DataFrame(subset_standard_umi_counter)

# Create a strip plot
plt.figure(figsize=(10, 5))
sns.swarmplot(data=df, x='barcode', y='count',size=3,color=".2")
plt.xticks([])
plt.title('Random distributions of UMI occurences')
plt.savefig('umi_swarmplot.pdf')