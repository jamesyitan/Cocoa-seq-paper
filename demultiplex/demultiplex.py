import os
import sys
from collections import defaultdict
from itertools import product, combinations
import configparser
import yaml
import pickle

#functions

def string_hamming_distance(str1, str2):
    """
    Fast hamming distance over 2 strings known to be of same length.
    In information theory, the Hamming distance between two strings of equal 
    length is the number of positions at which the corresponding symbols 
    are different.

    eg "karolin" and "kathrin" is 3.
    """
    dist = 0
    for i, j in zip(str1, str2):
        if i != j:
            dist += 1
    return dist

def from_fastq(handle):
    while True:
        name = next(handle).rstrip()[1:] #Read name
        seq = next(handle).rstrip() #Read seq
        next(handle) #+ line
        qual = next(handle).rstrip() #Read qual
        if not name or not seq or not qual:
            break
        yield name, seq, qual

# reverse complement function
def rev_comp(seq):
    tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(tbl[s] for s in seq[::-1])

def to_fastq_lines(bc, umi, seq, qual, read_id, v4_pos):
    """
    Return string that can be written to fastQ file
    """
    #Output V4 sequence after v4 primer
    seq_cut = seq[v4_pos:]
    qual_cut = qual[v4_pos:]
    return '@'+bc+':'+umi+':'+str(read_id)+'\n'+seq_cut+'\n+\n'+qual_cut+'\n'

# build barcode neighborhoods due to sequencing error
def build_barcode_neighborhoods(barcode_file):
    """
    Given a set of barcodes, produce sequences which can unambiguously be
    mapped to these barcodes, within 2 substitutions. If a sequence maps to 
    multiple barcodes, get rid of it. However, if a sequences maps to a bc1 with 
    1change and another with 2changes, keep the 1change mapping.
    Dict structure is with key being the barcode possibilities (including mutations)
    with values being the original barcode value
    """

    def seq_neighborhood(seq, n_subs=1):
        """
        Given a sequence, yield all sequences within n_subs substitutions of 
        that sequence by looping through each combination of base pairs within
        each combination of positions.
        """
        for positions in combinations(range(len(seq)), n_subs):
        # yields all unique combinations of indices for n_subs mutations
            for subs in product(*("ATGCN",)*n_subs):
            # yields all combinations of possible nucleotides for strings of length
            # n_subs
                seq_copy = list(seq)
                for p, s in zip(positions, subs):
                    seq_copy[p] = s
                yield ''.join(seq_copy)
    
    # contains all mutants that map uniquely to a barcode
    clean_mapping = dict()

    # contain single or double mutants 
    mapping1 = defaultdict(set)
    mapping2 = defaultdict(set)
    
    #Build the full neighborhood and iterate through barcodes
    with open(barcode_file, 'r', newline = '') as f:
        # iterate through each barcode (rstrip cleans string of whitespace)
        for line in f:
            barcode = rev_comp(line.rstrip()) #rc because the oligos provided in the list are the oligos used in extension so the actual sequences are the RC of those

            # each barcode obviously maps to itself uniquely
            clean_mapping[barcode] = barcode

            # for each possible mutated form of a given barcode, either add
            # the origin barcode into the set corresponding to that mutant or 
            # create a new entry for a mutant not already in mapping1
            # eg: barcodes CATG and CCTG would be in the set for mutant CTTG
            # but only barcode CATG could generate mutant CANG
            for n in seq_neighborhood(barcode, 1):
                mapping1[n].add(barcode)
            # same as above but with double mutants
            for n in seq_neighborhood(barcode, 2):
                mapping2[n].add(barcode) 
    
    # take all single-mutants and find those that could only have come from one
    # specific barcode
    for k, v in mapping1.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    
    for k, v in mapping2.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]

    del mapping1
    del mapping2
    return clean_mapping

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

##for individual reads
def process_reads(read, valid_bc1s={}, valid_bc2s={}):
    """
    Returns either:
        True, (barcode, umi)
            (if read passes filter)
        False, name of filter that failed
            (for stats collection)
    
    R1 anatomy: BBBBBBBB[BBB]WWWWWWWWWWWWWWWWWWWWWWCCCCCCCCUUUUUUSSSSSSSSSSSSSSSSSSSSPPPPPPPPPPPPPPPPPPP<read1>
        B = Barcode1, can be 8, 9, 10 or 11 bases long.
        W = 'W1' sequence, specified below
        C = Barcode2, always 8 bases
        U = UMI, always 6 bases
        S = CustomSeq adaptor, specified below
        P = 16S V4 forward primer, specified below
    """

    hamming_threshold_for_W1_matching = 3

    w1 = "GAGTGATTGCTTGTGACGCCTT"

    # 44-47 is the length of BC1+W1+BC2+UMI
    #BC1: 8-11 bases
    #W1 : 22 bases
    #BC2: 8 bases
    #UMI: 6 bases
    
    #Check for W1 adapter
    #Allow for up to hamming_threshold errors
    if w1 in read:
        w1_pos = read.find(w1)
        if not 7 < w1_pos < 12:
            # print_to_log(name)
            return False, 'No_W1'
    else:
        #Try to find W1 adapter at start positions 8-11
        #by checking hamming distance to W1.
        for w1_pos in range(8, 12):
            if string_hamming_distance(w1, read[w1_pos:w1_pos+22]) <= hamming_threshold_for_W1_matching:
                break
            else:
                return False, 'No_W1'
            
    bc2_pos=w1_pos+22
    umi_pos=bc2_pos+8
    adaptor_pos = umi_pos + 6
    primer16s_pos = adaptor_pos + 20
    v4start = primer16s_pos + 19

    # want to cut from start of bc1 to v4_pos ie seq[v4_pos:]
    bc1 = str(read[:w1_pos])
    bc2 = str(read[bc2_pos:umi_pos])
    umi = str(read[umi_pos:adaptor_pos])

    #Validate barcode (and try to correct when there is no ambiguity)
    if valid_bc1s and valid_bc2s:
        # Check if BC1 and BC2 can be mapped to expected barcodes
        if bc1 in valid_bc1s:
            # BC1 might be a neighboring BC, rather than a valid BC itself. 
            bc1 = valid_bc1s[bc1]
        else:
            return False, 'BC1'
    
        if bc2 in valid_bc2s:
            bc2 = valid_bc2s[bc2]
        else:
            return False, 'BC2'
    
    bc = '%s_%s'%(bc1, bc2)

    return True, (bc, umi, v4start)


#open configuration file
#put both barcode lists, param file, and fastq files in the same folder
dir_path = os.path.dirname(os.path.realpath(__file__)) #get the path of the current directory
config = configparser.ConfigParser()
config.read(dir_path+'/param') #read config file named param
anchor1 = config.get('default','anchor1')
umi_length = config.get('default','umi_length')
primer_seq = config.get('default', 'primer_sequence')
#sequencing_depth_threshold = config.get('default', 'sequencing_depth_threshold')

#input to terminal will be the name of the directory
dir_name = sys.argv[1]

#demultiplex based on barcodes (put information in sequence identifier information)

#start barcode metric counts
barcode_read_counter = defaultdict(int)

#counter for total reads
i = 0
#counter for kept reads
kept_reads = 0
#counter for reads thrown out due to no valid W2
no_W2_reads = 0
#counter for reads thrown out due to no valid barcode
no_valid_bc = 0  

#Get barcode neighborhoods
bc1s = build_barcode_neighborhoods(dir_path+'/barcode_lists/gel_barcode1_list')
bc2s = build_barcode_neighborhoods(dir_path+'/barcode_lists/gel_barcode2_list')

##open output file for writing
with open(dir_path+'/'+dir_name+'/filtered.fastq', 'w') as output_fastq:
    #counter for reads
    read_id = 0

    #for each line in the fastq file
    for r1_name, r1_seq, r1_qual in read_fastq('/'+dir_name+'/fwd.fastq'):
        
        #process each line to determine if it is a standard or V4 amplicon and also provide sequence information (bc, qual, sec)
        keep, result = process_reads(r1_seq, bc1s, bc2s)
        read_id+=1

        #keep track and output to terminal
        i += 1
        if i%1000000 == 0:
            print('%d reads parsed, kept %d reads' % (i, kept_reads))

        #report V4 amplicons to a fastq file
        if keep:
            bc, umi, v4_pos = result
            kept_reads += 1
            barcode_read_counter[bc] += 1
            output_fastq.write(to_fastq_lines(bc, umi, r1_seq, r1_qual, read_id, v4_pos))
        if not keep:
            if result == 'No_W1':
                no_W2_reads += 1
            elif (result == 'BC1') or (result == 'BC2'):
                no_valid_bc += 1

##output statistics
#output barcode data as a pickle
barcode_data = open(dir_path+'/'+dir_name+'/barcode_read_counts.pickle', 'wb') #needs to be binary write
pickle.dump(dict(barcode_read_counter), barcode_data)
barcode_data.close()

num_barcode_samps = len(barcode_read_counter.keys())

#output filtering statistics in output file
filtering_statistics = {
    'Total Reads' : i,
    'Valid Reads' : kept_reads,
    'Rejected Reads' : i - kept_reads,
    'Valid Fraction of Reads' : float(kept_reads)/i,
    'Total Droplets Analyzed' : num_barcode_samps,
    'Rejected Reads Due to Invalid W2' : no_W2_reads,
    'Rejected Reads Due to Invalid Barcode' : no_valid_bc

}

with open(dir_path+'/'+dir_name+'/filtering_statistics.yaml', 'w') as output_txt:
    yaml.dump(filtering_statistics, output_txt, default_flow_style=False)