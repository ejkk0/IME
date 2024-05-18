#!/usr/bin/env python

#SBATCH --partition defq
#SBATCH --job-name=generate_barcodes
#SBATCH --output=/scratch/users/ejkk/Cluster/IME/nsIME/generate_barcodes.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=128000
#SBATCH --mail-type=END
#SBATCH --mail-user=ejkk@mit.edu


###################################
### IMPORTS:
###################################


import numpy
import time
import random
import itertools
import Levenshtein as lev

import sys
sys.path.append("/home/ejkk/Cluster/IME/scripts")
from functions import *
from constants import *


###################################
### PATHS:
###################################


# path for file i/o:
outpath = '/scratch/users/ejkk/Cluster/IME/nsIME/'


###################################
### GENERATE BARCODES:
###################################

start = time.time()
printf('\nhello\n')

# params

maxlen = 5
mindist = 3
bc_len = 12
lib_size = 100000

# generate barcodes
printf('generating barcodes...\n')
barcodes = []

while len(barcodes) < lib_size:
	random_sequence=''
	for j in range(0,bc_len):
		random_sequence+=random.choice(dna)
	if ('A'*maxlen in random_sequence)|('T'*maxlen in random_sequence)|('C'*maxlen in random_sequence)|('G'*maxlen in random_sequence)|('ATG' in random_sequence):
		continue
	for b in barcodes:
		dist = lev.hamming(random_sequence,b)
		if dist<mindist:
			continue
	barcodes.append(random_sequence)
	if len(barcodes)%10000==0:
		printf(len(barcodes))

# save barcodes

list_to_csv(barcodes,outpath+'100k_12mer_barcodes.csv')

end = time.time()
printf('\nthis took %i seconds' % (end - start))

printf('\ndone!\n')
