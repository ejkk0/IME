#!/usr/bin/env python

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=preprocess_RNAseq_nbc
#SBATCH --output=/home/ejkk/Cluster/IME/outfiles/preprocess_RNAseq_nbc.out
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
### ARGS:
###################################


# Set up the ArgumentParser
parser = argparse.ArgumentParser(description='Process some variables.')

# Add arguments
parser.add_argument('round', type=str, nargs='?', default='RIL2', help='Which round of RNA seq, test, RIL2 or TRIS? (default: "test")')
parser.add_argument('run', type=int, nargs='?', default=1, help='Which run? (default: 1)')

# Parse the arguments
args = parser.parse_args()

rnd = args.round
run = args.run

###################################
### PATHS/DEFINITIONS/CONSTANTS:
###################################


# paths for i/o:
inpath = '/pool001/ejkk/RIL2/'
outpath = '/pool001/ejkk/RIL2/RNAseq_no_barcode_filters/'

# path to real reads or test reads:
if rnd=='test':
    printf('this is a test run.')
    fpath, rpath = inpath+'element/test_sets/elt_'+str(run)+'_1_100k.fastq', inpath+'element/test_sets/elt_'+str(run)+'_2_100k.fastq',
else:
    printf('this is NOT a test run!')
    printf(f'processing reads from {rnd} round {run}')
    if rnd == 'RIL2':
        fpath, rpath = inpath+'element/231107Bur_D23-16724-'+str(run)+'_1_sequence.fastq', inpath+'element/231107Bur_D23-16724-'+str(run)+'_2_sequence.fastq'
    if rnd == 'TRIS':
        fpath, rpath = inpath+'TRIS/231208Bur_D23-19001-'+str(run)+'_1_sequence.fastq', inpath+'TRIS/231208Bur_D23-19001-'+str(run)+'_2_sequence.fastq'

# sample indexes
indexes = ['GATC','AGCT','GTAT','CCAA','TTGG','GGCA',
           'AACC','CTTG','TATA','ATGC','CATT','TCCC',
           'GAGA','CGAT','ACTC','TGGA','AGTG','CTGA']
index_dict = {indexes[n]:(n+1) for n in range(18)}
#indexes.remove('GGCA')


###################################
### LOAD & FILTER READS:
###################################


# function to parse these fastqs, extract relevant features, filter out unwanted reads
# and produce dictionary of barcodes -> read counts for later processing

def parse_element_RNAseq(fpath, rpath):
    start = time.time()
    bcs, ids, reps, spliced, cryptic, freads, rreads = [], [], [], [], [], [], []
    reporter_fail = 0
    sample_fail = 0
    total = 0
    
    # first take in paired ends of each read at the same time and grab seqs
    
    with open(fpath, 'r') as fwd:
        with open(rpath, 'r') as rev:
            lines = []
            for fline,rline in zip(fwd,rev):              
                lines.append(fline.rstrip())
                lines.append(rline.rstrip())
                if len(lines) == 8:
                    total += 1
                    fwd = lines[2]    # grab fwd read
                    rev = lines[3]    # grab associated rev read
                    
                    # then apply filters to see if read contains expected sequences...
                    
                    if fwd[2:6] not in indexes:   # 1. toss read unless sample id matches valid sample id
                        sample_fail+=1
                    else:
                        current_read_sample_id = index_dict[fwd[2:6]]
                        
                        if rev[:6] not in reporter_id:   # 2. toss read unless reporter id matches GFP or dTom expected seq
                            reporter_fail+=1
                        else:
                            current_read_reporter = reporter_id[rev[:6]]
                            if current_read_reporter == 'GFP':
                                bc = rc(fwd[32:50])   # grab barcode of GFP read
                                splicecheck = rev[19:26]
                            elif current_read_reporter == 'dTom':
                                bc = rc(fwd[37:55])   # grab barcode of dTom read
                                splicecheck = 'TCGTGAA'

                            spl = splicecheck in int_bc_spacer   # all reads that do NOT match exon 2 are called "unspliced" for now
                            if spl:
                                cryp = False
                            else:
                                cryp = splicecheck not in fiveSS
                            ids.append(current_read_sample_id)
                            bcs.append(bc)
                            reps.append(current_read_reporter)
                            spliced.append(spl)
                            cryptic.append(cryp)
                            freads.append(fwd)
                            rreads.append(rev)
                    lines = []
                    
    reads_parsed = total
    printf('parsed %i reads' % (reads_parsed))
    printf('%i passed all filters (%.2f%%)' % (len(bcs),len(bcs)*100/total))
    printf('\nLOSSES:')
    printf('lost %i from no sample match' % (sample_fail))
    printf('lost %i from no reporter match' % (reporter_fail))
    printf('\nSTATS:')
    printf('%.2f%% of passing reads are from GFP' % (sum([r=='GFP' for r in reps])*100/len(bcs)))
    printf('%.2f%% of those are unspliced' % (sum([s==False for s in spliced])*100/sum([r=='GFP' for r in reps])))
    printf('%.2f%% of unspliced reads look like they could be cryptic splicing' % (sum([c==True for c in cryptic])*100/sum([s==False for s in spliced])))
    printf('%i unique barcodes in union of all samples' % (len(set(bcs))))
    GFP_bcs = [b for b,r in zip(bcs,reps) if r=='GFP']
    dTom_bcs = [b for b,r in zip(bcs,reps) if r=='dTom']
    printf('%i barcodes have at least 1 GFP and 1 dTom read' % (len(set(GFP_bcs).intersection(dTom_bcs))))
    end = time.time()
    printf('this took %i seconds' % (end - start))
    
    return bcs, ids, reps, spliced, cryptic, freads, rreads


###################################
### DO THE THING!:
###################################

printf('\nloading reads...')
bcs,ids,reps,spl,cryp,freads,rreads = parse_element_RNAseq(fpath,rpath)


###################################
### BUILD DICTS:
###################################


#printf('\nbuilding dictionary of barcodes:reads for each sample...')

RIL2_b2rc = {b:{n:[0,0,0] for n in np.arange(1,19)} for b in set(bcs)}
RIL2_b2unspl_canonical = {b:[] for b in set(bcs)}
RIL2_b2unspl_cryptic = {b:[] for b in set(bcs)}

for b,i,r,s,c,fr,rr in zip(bcs,ids,reps,spl,cryp,freads,rreads):
    if r == 'dTom':
        RIL2_b2rc[b][i][2] += 1
    if r == 'GFP':
        if s:
            RIL2_b2rc[b][i][0] += 1
        elif not c:
            RIL2_b2rc[b][i][1] += 1
            RIL2_b2unspl_canonical[b].append((fr,rr))
        elif c:
            RIL2_b2unspl_cryptic[b].append((fr,rr))

for b in set(bcs):
    if len(RIL2_b2unspl_canonical[b])==0:
        del RIL2_b2unspl_canonical[b]
    if len(RIL2_b2unspl_cryptic[b])==0:
        del RIL2_b2unspl_cryptic[b]

###################################
### PICKLES:
###################################


# store barcodes -> read counts dict
printf('storing dictionary...')
with open(outpath+rnd+'_' +str(run)+'_b2rc.pkl','wb') as f:
    pickle.dump(RIL2_b2rc, f)

# store barcodes -> unspliced reads dict
printf('storing unspliced...')
with open(outpath+rnd+'_'+str(run)+'_canonical_unspliced_dict_both_reads.pkl','wb') as f:
    pickle.dump(RIL2_b2unspl_canonical,f)
with open(outpath+rnd+'_'+str(run)+'_cryptic_unspliced_dict_both_reads.pkl','wb') as f:
    pickle.dump(RIL2_b2unspl_cryptic,f)

printf('\ndone!\n')
