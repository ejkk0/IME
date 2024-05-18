#!/usr/bin/env python3

###################################
#                                 #
#            FUNCTIONS            #
#                                 #
###################################

import sys
import os
import pandas as pd
import numpy as np
import numpy as np
import scipy as sp
from scipy.stats import mode
from scipy.spatial.distance import pdist, squareform
import time
import csv

###################################

# shortcuts

def vcs(x):
    return pd.Series(x).value_counts()

def dicthead(d, n=10):
    return list(d.items())[:n]
    
###################################
# print and flush to outfile

def printf(string):
    print(string, flush=True)

###################################

# pass sequences from string to int arrays

TO_INT = dict(A=0,C=1,T=2,G=3,N=4)
TO_INT['-'] = 5
TO_STR = {v:k for k,v in TO_INT.items()}
def toint(s):
    return [TO_INT[x] for x in s]
def tostr(i):
    return ''.join([TO_STR[x] for x in i])

def hamming(seq1, seq2):
    return int(pdist((toint(seq1), toint(seq2)), metric='hamming')*len(seq1))

###################################

# get consensus of sequences in int array format

def cons(seqs):
    return sp.stats.mode(seqs,axis=0)[0][0]

###################################

# flatten list

def flatten(l):
    return [item for sublist in l for item in sublist]

def list_to_csv(listname,filename):
  with open(filename, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    for val in listname:
      writer.writerow([val])

###################################

# parse fastq files into lists of read + index

def parse_fastq_to_lists(filename):
    start = time.time()
    introns = []
    bcs = []

    with open(filename, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 4:
                bcs.append((lines[0].split('#')[1].split('/'))[0])
                introns.append(lines[1])
                lines = []
    printf('parsed %i reads' % (len(bcs)))
    printf('%i unique barcodes' % (len(set(bcs))))
    printf('%i unique introns' % (len(set(introns))))
    end = time.time()
    printf('this took %i seconds' % (end - start))
    
    return bcs, introns

def parse_fastq_to_list(filename):
    start = time.time()
    introns = []

    with open(filename, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 4:
                introns.append(lines[1])
                lines = []
    printf('parsed %i reads' % (len(introns)))
    printf('%i unique introns' % (len(set(introns))))
    end = time.time()
    printf('this took %i seconds' % (end - start))
    
    return introns

def parse_BC_HiSeq(filename):
    start = time.time()
    bcs = []
    ids = []
    
    with open(filename, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 4:
                read_id = (lines[1][2:6])
                ids.append(read_id)
                bc = rc(lines[1][32:])
                bcs.append(bc)
                lines = []
                
    printf('parsed %i reads' % (len(bcs)))
    printf('%i unique barcodes' % (len(set(bcs))))
    end = time.time()
    printf('this took %i seconds' % (end - start))
    
    return bcs, ids


###################################

# parse alignment sam file

def parse_sam(fp):
    df = pd.read_csv(fp,
                     sep='\t',
                     header=None,
                     skiprows=4,
                     index_col=False,
                     usecols=[0,3,
                            #9,11,14
                            ])
    df.columns = ['intron_id','mapped_to',
                #'seq','counts','dists'
                ]
    #df.counts = [int(y) for y  in [df.counts[x].split(':')[2] for x in df.index]]
    #df.dists = [int(y) for y  in [df.dists[x].split(':')[2] for x in df.index]]

    #only add this part if you're sure of direct mapping!
    #df.mapped_to = [int(y) for y  in [x/170 for x in df.mapped_to]]

    return df

###################################

# get reverse complement of a sequence

tab = str.maketrans("ACTG", "TGAC")

def rc(seq):
    return seq.translate(tab)[::-1]

def gc(seq):
    return (seq.upper().count('G')+seq.upper().count('C'))/len(seq)

def percent_N(N,seq):
    return (seq.upper().count(N))/len(seq)

def count_N(N, seq):
    count = 0
    for x in range(len(seq) - len(N) + 1):
        if seq[x:x+len(N)] == N:
            count += 1
    return count

###################################

# get consensus of a list of sequences

def consensus(reads,length):  
    cons = ''
    nt=0
    while nt < length:
        ls = [read[nt] for read in reads]
        cons=cons+max(ls,key=ls.count)
        nt+=1
    return cons

###################################

# get all one-off seqs from a given seq (includes all bases plus N)

bases = ['A','C','G','T','N']

def one_offs(seq):
    one_off_list = []
    for i in range(len(seq)):
        for base in bases:
            if base != seq[i]:
                new_seq = seq[:i]+base+seq[i+1:]
                one_off_list.append(new_seq)
    return one_off_list

def two_offs(seq):
    return flatten([one_offs(x) for x in one_offs(seq)])

bc_bases = ['A','C','T']

def bc_one_offs(seq):
    one_off_list = []
    for i in range(len(seq)):
        for base in bc_bases:
            if base != seq[i]:
                new_seq = seq[:i]+base+seq[i+1:]
                one_off_list.append(new_seq)
    return one_off_list

####################################

# extract counts of all kmers in a sequence

def extract_kmers(intron,k):
    hex_dict = {}
    for j in range(0,len(intron)-(k-1)):
        kmer = intron[j:(j+k)]
        if 'N' not in kmer:
            if kmer not in hex_dict:
                hex_dict.update({kmer:1})
            else:
                hex_dict[kmer]+=1
    return hex_dict

#####################################

# misc functions

def de_nan(df):
    return df[np.isfinite(df).all(1)]

def longest_hairpin(stringy):
    # stringy has to be the structure output of RNA.fold
    s2 = [s for s in stringy.split('(') if len(s)>0]
    s3 = [s for s in flatten([x.split('.') for x in s2]) if len(s)>0]
    ls = [len(t) for t in s3]
    return max(ls)

#####################################



