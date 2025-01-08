import os
import sys
import pickle
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

dir_path = os.path.dirname(os.path.realpath(__file__)) #get the path of the current directory
#input to terminal will be the name of the directory to process (should be the same as the one you demultiplexed)
dir_name = sys.argv[1]

#load pickle data files
with open(dir_path+'/'+dir_name+'/barcode_read_counts.pickle', 'rb') as barcode_data:
    barcode_read_counter = pickle.load(barcode_data)

#cumulative read per sample curve
cum_reads = [0]
droplet_count = [0]
droplet_count_sofar = 0
total_reads_sofar = 0

for bc, count in barcode_read_counter.items():
    droplet_count_sofar = droplet_count_sofar + 1
    total_reads_sofar = total_reads_sofar + count
    droplet_count.append(droplet_count_sofar)
    cum_reads.append(total_reads_sofar)

ax = plt.subplot(111)
ax.plot(droplet_count,cum_reads,color='black')
ax.set_xlabel('Number of cumulative barcode samples')
ax.set_ylabel('Total number of reads')
plt.savefig('%s_curve.pdf' % dir_name, transparent=True)