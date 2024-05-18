#!/usr/bin/env python

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=parse_nsIME_RNAseq
#SBATCH --output=/nobackup1/ejkk/Cluster/IME/nsIME/RNAseq/parse_nsIME_RNAseq.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=256000
#SBATCH --mail-type=END
#SBATCH --mail-user=ejkk@mit.edu


###################################
### IMPORTS:
###################################


import sys
import time
import pickle
import pandas as pd
import numpy as np

import sys
sys.path.append("/home/ejkk/Cluster/IME/scripts")
from functions import *
from constants import *


###################################
### PATHS/DEFINITIONS/CONSTANTS:
###################################


# path to dump outputs:
outpath = '/nobackup1/ejkk/Cluster/IME/nsIME/RNAseq/'

# working data directory:
wdd = '/home/ejkk/Cluster/IME/nsIME/'

# paths to real reads:
fpath, rpath = outpath + '220309Bur_D22-3612_1_sequence.fastq', outpath + '220309Bur_D22-3612_2_sequence.fastq'

# paths to test reads:
#fpath, rpath = outpath + 'test_1m_fwd.fastq', outpath + 'test_1m_fwd.fastq'


###################################
### OTHER STUFF:
###################################


# import "final pool" of introns ordered in 300bp library and assoc. barcodes
FP = pickle.load(open(wdd+'nsIME_final_oligo_pool_10-15-20.pkl', 'rb'))

# import long intron pool
LOP = pickle.load(open(wdd+'nsIME_final_long_oligo_pool_10-21-20.pkl', 'rb'))

# merge into one table and get list of barcodes
IP = pd.concat([FP, LOP])
IP0 = IP[IP.PAS==0]
barcodes = list(IP0.bc)

# expand list of barcodes to all one-offs without Ns and store as dict mapping back to original bc
barcode_oos = {}
for b in barcodes:
    barcode_oos[b] = b
    for o in one_offs(b):
        barcode_oos[o] = b


# make dict of sample ids and one-offs mapped to corresponding samples
sample_ids = {}

def add_sample_oos(sample_dict, sample_id, sample_name):
    for x in [sample_id]+[i for i in one_offs(sample_id)]:
        sample_dict[x] = sample_name 
    
add_sample_oos(sample_ids, 'ATCACG', 'RNA1')
add_sample_oos(sample_ids, 'CGATGT', 'RNA2')
add_sample_oos(sample_ids, 'TTAGGC', 'RNA3')
add_sample_oos(sample_ids, 'TGACCA', 'DNA1')
add_sample_oos(sample_ids, 'ACAGTG', 'DNA2')
add_sample_oos(sample_ids, 'GCCAAT', 'DNA3')


###################################
### LOAD READS:
###################################


# function to parse reads

def parse_nsIME_RNAseq(fpath, rpath):
    start = time.time()
    umis = []
    bcs = []
    ids = []
    reporters = []
    freads = []
    rreads = []
    uhoh1 = 0
    uhoh2 = 0
    uhoh3 = 0
    uhoh4 = 0
    uhoh5 = 0
    
    with open(fpath, 'r') as fwd:
        with open(rpath, 'r') as rev:
            lines = []
            for fline,rline in zip(fwd,rev):
                lines.append(fline.rstrip())
                lines.append(rline.rstrip())
                if len(lines) == 8:
                    fwd = lines[2]    # grab fwd read
                    rev = lines[3]    # grab associated rev read

                    umi = fwd[:9]    # grab umi from fwd read
                    
                    if 'TGACGT' not in fwd[:30]:    # check fwd read for 6mer that should precede barcode
                        uhoh1+=1
                    else:
                        bc = rc(fwd.split('TGACGT')[1][:12])    # if present, grab barcode
                        if bc not in barcode_oos:
                            uhoh2+=1
                        else:
                            reporter = rc(fwd.split('TGACGT')[1][12:17])    # get 5mer following barcode to check GFP/dTom
                            if reporter not in ['ATCGG','GTGAC']:
                                uhoh3+=1
                            elif 'CTGTGA' not in rc(fwd):    # check for 6mer following sample index
                                uhoh4+=1
                            else:
                                sample = rc(fwd).split('CTGTGA')[0][-6:]    # if present grab sample index
                                if sample not in sample_ids:
                                    uhoh5+=1
                                else:
                                    bcs.append(bc)
                                    reporters.append(reporter)
                                    umis.append(umi)
                                    ids.append(sample)
                                    freads.append(fwd)
                                    rreads.append(rev)
                    lines = []
    
    reads_parsed = len(ids)+np.sum([uhoh1,uhoh2,uhoh3,uhoh4,uhoh5])
    printf('parsed %i reads' % (reads_parsed))
    printf('%i passed all filters' % (len(bcs)))
    printf('\nLOSSES:')
    printf('lost %i due to no match for 6mer preceding barcode' % (uhoh1))
    printf('lost %i due to invalid barcode' % (uhoh2))
    printf('lost %i due to invalid reporter id' % (uhoh3))
    printf('lost %i due to no match for 6mer preceding sample id' % (uhoh4))
    printf('lost %i due to invalid sample id' % (uhoh5))
    end = time.time()
    printf('this took %i seconds' % (end - start))
    
    return ids, bcs, umis, reporters, freads, rreads

###################################
### BUILD OUTPUT:
###################################

printf('\nloading reads...\n')
ids, bcs, umis, reporters, freads, rreads = parse_nsIME_RNAseq(fpath,rpath)
nsIME_data = zip(ids, bcs, umis, reporters, freads, rreads)


###################################
### PICKLES:
###################################

# store output
printf('storing all outputs as a pkl...\n')
pickle.dump(nsIME_data,open(outpath+'all_of_it.pkl','wb'))

#printf('let\'s do a csv too...')

#import csv

#with open(outpath+'all_of_it.csv', "w") as f:
#    writer = csv.writer(f)
#    for row in nsIME_data:
#        writer.writerow(row)

printf('done!\n')