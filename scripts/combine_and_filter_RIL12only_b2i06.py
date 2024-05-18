#!/usr/bin/env python

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=combine_and_filter_RIL12only_b2i06
#SBATCH --output=/home/ejkk/Cluster/IME/outfiles/combine_and_filter_RIL12only_b2i06.out
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
import argparse

import sys
sys.path.append("/home/ejkk/Cluster/IME/scripts")
from functions import *
from constants import *
printf('hello!')


###################################
### PATHS/DEFINITIONS/CONSTANTS:
###################################


# paths for i/o:
IME_path = '/pool001/ejkk/RIL2/'
reads_path = IME_path+'RNAseq_no_barcode_filters/no_spl_filter/'

test = False

###################################
### LOAD & FILTER READS:
###################################


printf('loading barcode info...')

# load in barcode-to-intron dictionary from round 2 (bottlenecked)
with open(IME_path+'b2i_06_with_spikeins_and_flanking_seq.pkl', 'rb') as f:
    b2i = pickle.load(f)

printf('barcode info loaded!')
printf('loading reads...')

# load in data from RIL1 and RIL2
# stored as a dictionary of barcodes->samples->[GFP, unspliced, dTomato] counts

with open(IME_path+'RNAseq_no_barcode_filters/2019_screen_b2r.pkl', 'rb') as f:
    RIL1 = pickle.load(f)
with open(reads_path+'RIL2_1_b2rc.pkl', 'rb') as f:
    RIL2_1 = pickle.load(f)
with open(reads_path+'RIL2_2_b2rc.pkl', 'rb') as f:
    RIL2_2 = pickle.load(f)

printf('reads loaded!')
printf('filtering reads...')

# let's merge the 3 dicts and do error-correcting for barcodes
from Levenshtein import distance as lev

RIL1_bulk_samples = list(range(0,5))
RIL2_bulk_samples = list(range(1,6))
FACS_samples = list(range(7,19))

RIL12_bulk = {b:{n:np.array([0,0,0]) for n in list(range(1,11))} for b in b2i}
RIL2_FACS = {b:{(n-6):np.array([0,0,0]) for n in FACS_samples} for b in b2i}

def print_status(list_of_lists, checkpoint):
    bcs_done = np.sum(list_of_lists)
    if (bcs_done+1) % checkpoint == 0:
        printf(f'done {bcs_done}...')
        return True
    else:
        return False


for i, RIL_dict in enumerate([RIL1, RIL2_1, RIL2_2]):

    printf(f'\nprocessing reads from dict {i+1}, containing {len(RIL_dict)} barcodes')
    matched, saved, multimap, nomap = 0,0,0,0

    for b in RIL_dict:
        if print_status([matched, saved, multimap, nomap], 25000) & test:
            break

        if b in b2i:
            if i==0:    # if this is RIL1
                for n in RIL1_bulk_samples:
                    RIL12_bulk[b][n+1] += RIL_dict[b][n]
            if i>1:    # if this is RIL2
                for n in RIL2_bulk_samples:
                    RIL12_bulk[b][n+5] += RIL_dict[b][n]
                for n in FACS_samples:
                    RIL2_FACS[b][n-6] += RIL_dict[b][n]
            matched+=1

        if b not in b2i:
            dists = [(lev(b, bc), bc) for bc in b2i if lev(b, bc) <= 2]
            if len(dists) == 0:
                nomap+=1
            elif len(dists)>1:
                multimap+=1
            elif len(dists) == 1:
                new_b = dists[0][1]

                if i==0:    # if this is RIL1
                    for n in RIL1_bulk_samples:
                        RIL12_bulk[new_b][n+1] += RIL_dict[b][n]
                if i>1:    # if this is RIL2
                    for n in RIL2_bulk_samples:
                        RIL12_bulk[new_b][n+5] += RIL_dict[b][n]
                    for n in FACS_samples:
                        RIL2_FACS[new_b][n-6] += RIL_dict[b][n]
                saved+=1
            
    print(f'{matched} barcodes recognized from dict')
    print(f'{saved} barcodes saved from one- or two-off')
    print(f'{nomap} barcodes could not be mapped')
    print(f'{multimap} barcodes could not be mapped uniquely')
      
for b in b2i:
    if np.sum([x for x in RIL12_bulk[b].values()])==0:
        del RIL12_bulk[b]
    if np.sum([x for x in RIL2_FACS[b].values()])==0:
        del RIL2_FACS[b]

printf('reads filtered!')


###################################
### PICKLES:
###################################


printf('storing filtered b2r dicts...')

with open(reads_path+'RIL12_bulk_b2i06_b2r.pkl','wb') as f:
    pickle.dump(RIL12_bulk, f)
with open(reads_path+'RIL2_FACS_b2i06_b2r.pkl','wb') as f:
    pickle.dump(RIL2_FACS, f)

printf('done! :)')
