#!/usr/bin/env python

#SBATCH --partition sched_mit_cburge
#SBATCH --job-name=get_nanopore_reads_bcs_ints
#SBATCH --output=/home/ejkk/Cluster/IME/outfiles/get_nanopore_reads_bcs_ints.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16000
#SBATCH --mail-type=END
#SBATCH --mail-user=ejkk@mit.edu

# this is a modified version of my original nanopore analysis script to extract
# barcodes and introns from the minimap2-aligned reads based on their location

# I created this version on 11/1/23 to reduce the number of things this script
# attempts to accomplish - it now returns a simplified table as follows:

# for each read that maps contiguously from dTomato barcode through GFP barcode,
# get the read id, the read length, the dTom barcode sequence, the GFP barcode 
# sequence, and the intron sequence

# I will output these as a .csv

# 11/10/23: edited script to make following changes:
# - measure length of query seq rather than length of alignment for read_len
# - filter on levenshtein distance <=3 rather than identical barcodes
# - output fasta of filtered_introns so i can feed it directly into usearch


###################################
### IMPORTS:
###################################


import sys
import csv
import time
import pickle
import pysam
import pandas as pd
import numpy as np
from Levenshtein import distance as lev
sys.path.append("/home/ejkk/Cluster/IME/scripts")
from functions import *
from constants import *
printf('hello hello')


###################################
### PATHS & CONSTANTS:
###################################


# path to dump outputs:
outpath = '/pool001/ejkk/RIL2/'

# path to bam file:
bampath = outpath+'RIL_ORI_alignment_2.sorted.bam'

# test bam
#bampath = '/orcd/nese/cburge/001/ejkk/RIL2/test_10K_alignment.sorted.bam'

# coordinates of barcodes and intron in plasmid map
dt_start, dt_end = 1918, 1935
g_start, g_end = 3837, 3854
i_start, i_end = 3627, 3786

colnames = ['read_id','read_len','dt_bc','g_bc', 'bc_dists', 'intron']
max_length = 8900

###################################
### FUNCTIONS:
###################################

# the main event
def get_raw_barcode_intron_table(bam_file):
    read_ids, read_lens, dt_bcs, g_bcs, bc_dists, introns, filtered_introns = [], [], [], [], [], [], []
    total_reads = 0
    
    # iterate through reads, taking only reads that overlap random regions
    with pysam.AlignmentFile(bampath, "rb", check_sq=False) as bam_file:
        for alignment in bam_file.fetch('RIL_ORI',dt_start,g_end):
            total_reads+=1

            if alignment.reference_start<=dt_start and alignment.reference_end>=g_end:
                read_id = alignment.query_name
                read_seq = alignment.query_sequence

                # for each region of interest get coords in read that map to desired window in reference
                if read_seq != None:
                    read_len = len(alignment.query_sequence)
                    aln_pairs = alignment.get_aligned_pairs(with_seq=True)
                    dt_bases = [x for (x,y,z) in aln_pairs if y != None and y >= dt_start and y <= dt_end]
                    dt = ''.join(['-' if x == None else read_seq[int(x)] for x in dt_bases])
                    g_bases = [x for (x,y,z) in aln_pairs if y != None and y >= g_start and y <= g_end]
                    g = ''.join(['-' if x == None else read_seq[int(x)] for x in g_bases])
                    i_bases = [x for (x,y,z) in aln_pairs if y != None and y >= i_start and y <= i_end]
                    intron = ''.join(['-' if x == None else read_seq[int(x)] for x in i_bases])
                    bc_dist = lev(g,dt)

                    # return results
                    read_ids.append(read_id)
                    read_lens.append(read_len)
                    dt_bcs.append(dt)
                    g_bcs.append(g)
                    bc_dists.append(bc_dist)
                    introns.append(intron)

                    # check if barcodes match AND read is not concatemer
                    if read_len < max_length and bc_dist<=3 and '-' not in g and '-----' not in intron:
                        filtered_introns.append((g,intron))

    printf('total alignments considered: %i' % total_reads)
    printf('alignments mapping to barcodes + intron: %i' % len(read_ids))
    printf('alignments passing filters: %i' % len(filtered_introns))
    return read_ids, read_lens, dt_bcs, g_bcs, bc_dists, introns, filtered_introns

###################################
### DO THE THING:
###################################

printf('\nparsing reads to get barcode and intron sequences...')
read_ids, read_lens, dt_bcs, g_bcs, bc_dists, introns, filtered_introns = get_raw_barcode_intron_table(bampath)
printf('done!\n')

# export .csv of all reads
printf('storing csv of all barcodes and introns from reads...')
csv_file = outpath+'raw_barcode_intron_table.csv'
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(colnames)
    for row in zip(read_ids, read_lens, dt_bcs, g_bcs, bc_dists, introns):
        writer.writerow(row)

# export fasta of filtered_introns only
printf('storing fasta of barcodes and introns from trusted reads only...')
fasta_file = outpath + 'filtered_introns.fa'

with open(fasta_file, 'w') as file:
    for read_id, read_sequence in filtered_introns:
        file.write(f'>{read_id}\n{read_sequence}\n')

printf('done!\n')
