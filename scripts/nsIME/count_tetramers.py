#!/usr/bin/env python

#SBATCH --partition sched_mem1TB
#SBATCH --job-name=count_4mers
#SBATCH --output=/scratch/users/ejkk/Cluster/IME/nsIME_screen/count_4mers.out
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
import pickle
import itertools
import argparse

import sys
sys.path.append("/home/ejkk/Cluster/IME/scripts")
from functions import *
printf('\nhello')


###################################
### PATHS:
###################################


# path for file i/o:
outpath = '/scratch/users/ejkk/Cluster/IME/nsIME_screen/'


###################################
### BUILD PENTAMER DICT:
###################################


start = time.time()

# load intron seqs
printf('\nloading intron sequences...')
with open(outpath+'intseqs.pkl','rb') as f:
	intseqs = pickle.load(f)

# enumerate pentamers
k=4
dna = ["A","G","C","T"]
pentamers = sorted([''.join(p) for p in itertools.product(dna, repeat=k)])

# generate intron pentamer tables
printf('generating tetramer tables...')
intron_pentamers = pd.DataFrame(index=list(intseqs.keys()), columns=pentamers)
ints = list(intron_pentamers.index)

for p in pentamers:
    p_counts = []
    for i in ints:
        p_counts.append(intseqs[i].count(p))
    intron_pentamers[p] = p_counts
    
end = time.time()
printf('\nthis took %i seconds\n' % (end - start))

with open(outpath+'ns_tetramer_counts.pkl','wb') as f:
    pickle.dump(intron_pentamers, f)

printf('done!\n')