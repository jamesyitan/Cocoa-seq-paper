##This script is to characterize the distribution of reads/barcode and then remove noisy reads

##run in an interactive job on high performance computing cluster
##salloc --account=ninalin1 --mem=16G --time=4:00:00 --partition=standard

import os
import sys
import configparser
import pickle
import matplotlib.pyplot as plt
import math

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def hist(dict, samp_name, bin_size):
    bottom=0.5
    barcode_occurrences = list(dict.values())

    # Create subplots with two rows (for horizontal layout)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))

    # Calculate the number of bins based on the bin size
    num_bins = math.ceil(int(max(barcode_occurrences) / bin_size))

    # Plot the first histogram on the top subplot
    _, bins, _ = axs[0].hist(barcode_occurrences, num_bins, log=True,bottom=bottom,color='black',edgecolor='black',linewidth=1)
    axs[0].set_ylabel('Frequency')

    # Filter values falling into the lowest bin
    #set threshold for the second histogram
    threshold=100
    lowest_values = [value for value in barcode_occurrences if value <= threshold]
    new_num_bins = math.ceil(int(max(lowest_values) / 1))

    # Plot the second histogram on the bottom subplot
    _, bins2, _ = axs[1].hist(lowest_values, bins=new_num_bins, log=True,bottom=bottom,color='black',edgecolor='black',linewidth=1)
    axs[1].set_xlabel('reads per barcode family')
    axs[1].set_title(samp_name)

    # Save the figure
    fig.savefig('%s_hist.pdf' % samp_name, transparent=True)

    num_barcodes_tot = len(list(dict.keys()))
    return num_barcodes_tot

def filter(thres):
    filtered_barcode_counter = {key: value for key, value in barcode_read_counter.items() if value >= thres}
    filtered_barcode = list(filtered_barcode_counter.keys())
    return filtered_barcode

def read_fastq(fwd_fastq):
    #read in fastq file
    read_stream = open(dir_path+'/'+str(fwd_fastq), 'r')

    while True:
        #Read 4 lines from each FastQ
        name = read_stream.readline().rstrip() #Read name
        r1_seq = read_stream.readline().rstrip() #Read seq
        read_stream.readline().rstrip() #+ line
        r1_qual = read_stream.readline().rstrip() #Read qual
        
        # changed to allow for empty reads (caused by adapter trimming)
        if name:
            yield name, r1_seq, r1_qual
        else:
        # if not r1_seq or not r2_seq:
            break

    read_stream.close()

#open configuration file
#put both barcode lists, param file, and fastq files in the same folder
dir_path = os.path.dirname(os.path.realpath(__file__)) #get the path of the current directory
config = configparser.ConfigParser()
config.read(dir_path+'/param') #read config file named param

#input to terminal will be the name of the directory
dir_name = sys.argv[1]

with open(dir_path+'/'+dir_name+'/barcode_read_counts.pickle','rb') as pkl_file:
    barcode_read_counter = pickle.load(pkl_file)
    num_barcodes = hist(barcode_read_counter, dir_name, 10)
    print(num_barcodes)

    while True:
        threshold = int(input("Looking at the histogram, what sequencing depth (reads/barcode) do you want to set the minimal threshold? "))
        #For even, set to 1000; num_barcodes from 30278 to 1294
        #For log1, set to 4000; num_barcodes from 59014 to 1178
        #For log2, set to 4500; num_barcodes from 69669 to 1485
        #For ecbs, set to 5000; num_barcodes from 109459 to 2464
        barcode_above_thres = filter(threshold)
        print('This will reduce the number of barcodes from %i to %i.\n' %(num_barcodes,len(barcode_above_thres)))
        confirm = input('Is that okay? (Type "Yes" if so,) ').strip().lower()

        if confirm == 'yes':
            print("You said 'Yes'. Filtering fastq file...")
            break
        elif confirm == 'no':
            print("You said 'No'.")
        else:
            print("Invalid input. Please enter 'Yes' or 'No'.")

##read in filtered_fastq file
#if the fastq name has the correct barcode name, keep it, if not, put it into another file to be looked at for analysis
with open(dir_path+'/'+dir_name+'/above_threshold.fastq','w') as output_fastq, open(dir_path+'/'+dir_name+'/below_threshold.fastq','w') as trash_fastq:
    for r1_name, r1_seq, r1_qual in read_fastq('/'+dir_name+'/filtered.fastq'):
        barcode_index = r1_name.find(":")
        barcode = r1_name[1:barcode_index]
        if barcode in barcode_above_thres:
            output_fastq.write(r1_name+"\n"+r1_seq+"\n+\n"+r1_qual+"\n")
        else:
            trash_fastq.write(r1_name+"\n"+r1_seq+"\n+\n"+r1_qual+"\n")

#push through vsearch
#./vsearch --fastx_filter even/above_threshold.fastq --fastaout even/even.good.fasta --fastqout_discarded even/even.bad.fastq
#./vsearch --fastx_filter log1/above_threshold.fastq --fastaout log1/log1.good.fasta --fastqout_discarded log1/log1.bad.fastq
#./vsearch --fastx_filter log2/above_threshold.fastq --fastaout log2/log2.good.fasta --fastqout_discarded log2/log2.bad.fastq
#./vsearch --fastx_filter ecbs/above_threshold.fastq --fastaout ecbs/ecbs.good.fasta --fastqout_discarded ecbs/ecbs.bad.fastq

##For below threshold(trash) sequences
#./vsearch --fastx_filter even/below_threshold.fastq --fastaout even/even.trash.fasta --fastqout_discarded even/even.trash.bad.fastq
#./vsearch --fastx_filter log1/below_threshold.fastq --fastaout log1/log1.trash.fasta --fastqout_discarded log1/log1.trash.bad.fastq
#./vsearch --fastx_filter log2/below_threshold.fastq --fastaout log2/log2.trash.fasta --fastqout_discarded log2/log2.trash.bad.fastq
#./vsearch --fastx_filter ecbs/below_threshold.fastq --fastaout ecbs/ecbs.trash.fasta --fastqout_discarded ecbs/ecbs.trash.bad.fastq

##Run on multiple threads
##salloc --account=ninalin1 --cpus-per-task=4 --mem-per-cpu=4gb --time=4:00:00 --partition=standard
#./vsearch --usearch_global even/even.good.fasta --db mock_community_v4_ref.fasta  --mothur_shared_out even/even.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 1807622 of 1870277 (96.65%)
#./vsearch --usearch_global log1/log1.good.fasta --db mock_community_v4_ref.fasta  --mothur_shared_out log1/log1.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 11118631 of 11450686 (97.10%)
#./vsearch --usearch_global log2/log2.good.fasta --db mock_community_v4_ref.fasta  --mothur_shared_out log2/log2.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 8314130 of 8915003 (93.26%)
#./vsearch --usearch_global ecbs/ecbs.good.fasta --db ecbs_v4_ref.fasta  --mothur_shared_out ecbs/ecbs.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 22683747 of 24479988 (92.66%)

#./vsearch --usearch_global even/even.trash.fasta --db mock_community_v4_ref.fasta  --mothur_shared_out even/even.trash.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 1772244 of 1839383 (96.35%)
#./vsearch --usearch_global log1/log1.trash.fasta --db mock_community_v4_ref.fasta  --mothur_shared_out log1/log1.trash.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 3540024 of 3876692 (91.32%)
#./vsearch --usearch_global log2/log2.trash.fasta --db mock_community_v4_ref.fasta  --mothur_shared_out log2/log2.trash.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 8030263 of 8593076 (93.45%)
#./vsearch --usearch_global ecbs/ecbs.trash.fasta --db mock_community_v4_ref.fasta  --mothur_shared_out ecbs/ecbs.trash.shared  --maxseqlength 260 --minseqlength 100 --id 0.95 --threads 8
# Matching unique query sequences: 8177972 of 8816057 (92.76%)