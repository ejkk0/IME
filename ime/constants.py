#!/usr/bin/env python3

###################################
#                                 #
#            CONSTANTS            #
#                                 #
###################################

import sys
import time
import pickle
import pandas as pd
import numpy as np
import scipy as sp
import itertools
#import Levenshtein as lev
#import matplotlib.pyplot as plt
#from matplotlib.ticker import PercentFormatter
#from scipy.spatial.distance import pdist, squareform
sys.path.append("/home/ejkk/Cluster/IME/scripts")
from functions import *

###################################

# DNA

bases = ['A','C','G','T','N']
dna = ["A","G","C","T"]
dna_H = ["A","C","T"]

sevenmers = sorted([''.join(p) for p in itertools.product(dna, repeat=7)])
hexamers = sorted([''.join(p) for p in itertools.product(dna, repeat=6)])
pentamers = sorted([''.join(p) for p in itertools.product(dna, repeat=5)])
tetramers = sorted([''.join(p) for p in itertools.product(dna, repeat=4)])

# reporter id (first 6nt of reverse read)

reporter_id = {'GCTAGC':'dTom',
              'TGATCG':'GFP'}
reporter_id.update({OO:'dTom' for OO in one_offs('GCTAGC')})
reporter_id.update({OO:'GFP' for OO in one_offs('TGATCG')})

# seqs of constant region at beginning of intron and at beginning of 3' exon (for splicing analysis)

int_bc_spacer = one_offs('TCGTGAA')
int_bc_spacer.append('TCGTGAA')

fiveSS = one_offs('AGTAGCG')
fiveSS.append('AGTAGCG')

# barcodes of intronless and int-containing spike in controls

controls = {'pEK20':'TATATTATATCTTCTTTA',   # intronless 1
            'pEK21':'CCCCCTCTCATCCCACAA',   # intronless 2
            'pEK22':'CTCATAATCACTTTTCCT',   # AdmL-IgG intron
            'pEK23':'ATACCATAACTACCCACT'}   # UbC intron

random_ints = {'ACACTACTAATATACTCA':('CTTTATCTTGACACCCTGGGGTAGTGTGTTGCTAGCCTTTTGTAGTTTATTGTATTTCACGCATAGTACCACACTGACCGTTAAGTGATGGGACCGTTGATATCCTCGAAGGAAAATGGCAATGATAATGCAGTCCGTACGCTGCGTGCACTTGATGGGG',1),
           'ATCACTTTCTATCCCCTT':('TTCTGACACTGAAATTGAAGAGACGTCCCTGATAAGGCCCTGCTACATCTCATGGCACTTGACCCGGCGTTAGGCTTCATTTGAGTAACAAGATTGTTCCCATGATTAGGTCGTTCCTGCTCCTGGGGTGCGGGCGTCATACGAGTTACTAGCGTTGTGC',1),
           'CCTCAAACTTTAACCACC':('CTGTAGGTTGTAGGGGTCCCGCGTCTGTCGGCTATACGGTGGATCGAACATTTGGTTTCGTTCATGTCAAAAATACGCGTTAAGGGCTTGGTTAAAGGTGGCGCGTCTTTAATGTCGTAATAGCTTGGATGTATAGGGGGGACTGGGAGGCGACACCGTA',1),
           'CCATCATAAATTCTCAAA':('AACAGGGCATCTTCGCTCAAAGAATACTTCAAGACGACGTTATCAAGGCGCGACCTGTCCCGATTTACTTGAAGGACGTGAGTTGCCTTTGCCTATTCACGACCAAACATCCAGTCAAGGGCGTTCGCCTTTGGTTGCGGTTGCTCCTTACGTGCCTCAT',1),
           'ATACTTAACCATACCCTA':('TCCCGCCAAATGGTAGGAACCCCTATCATTCGCGATAACTCTGGCTTATATGACGAGAACTGTTTAATGGGAATAGAGACGTCACTGCGAAATGGTAATGCCACATATGTTCACTTATAGCGGTAGTGCGATGAGGAAGTTGGAAGCACGCAAATCC',1),
           'CCACATTAACACCAATCA':('GCCAATGAGTAGTTTAGTGGTTCATAGAGCTTAAGGAATCGGACTGTTTCAAAGAAAAAATGTTCACTTGCATATAAGTGCTGTGAAATAATGGTCTTAACAATCAGAACCGATACCCGGATCCACGAGAAGCTCGCGATTTAAGTGTTCGAGTTTTGGG',1),
           'CCCACTCCTCTCTTTCAC':('TTACAGGACATAATTTGGTAATTTGTTGCAAAGAGCGGTGCGGGTTCAGGACATCATCGGGAAGCGGGGTAGGGTCAAACAAGAACTCGTTATGTGTCAACATTAGGGGTGGTGTTCGAGGTGGACCGGCCTCGTGGAGGGTCACGGACTCGTGGGCCGT',1),
           'TCTAATACATACCATCCA':('ACAGTGCTCAACGATGGGGGTCGCTAGGCGTGCTTCCGTTATAGGCATGAGTTTTACGCGTGGCGGGGATTCACCGTGCTTACTGTTTAATTTTTGAGTAGTCACCCGGAGTTAGATAGACTTATTGTTTCTTTAGGGAGGCACACCCGCGCTTTCTTGG',1),
           'CCCTTAACCTTTATATAA':('CTTAGGCGAATCACCCACATAGTCTACTATGCGTTTTAAGGGTAGAAATCTGCTTACATGACGTTTATATGGCTATTGGGACCTAGCTCGAGCTTTGCGCGTCATCTATCTGGTGAAGCTTCAGAATGTACAGTGTGCCCGAACTACGCAAGGAGGTGTG',1)}

intname = {'pEK20':'intronless 1',
           'pEK21':'intronless 2',
           'pEK22':'AdmL-IgG',
           'pEK23':'UbC',
           'RIL':'Random Library',
          'RILnoLV':'Random Library',
          'neg':'Negative'}

UbC = 'ctgctgggctggccggggctttcgtggccgccgggccgctcggtgggacggaggcgtgtggagagaccgccaagggctgtagtctgggtccgcgagcaaggttgccctgaactgggggttggggggagcgcagcaaaatggcggctgttcccgagtcttgaatggaagacgcttgtgaggcgggctgtgaggtcgttgaaacaaggtggggggcatggtgggcggcaagaacccaaggtcttgaggccttcgctaatgcgggaaagctcttattcgggtgagatgggctggggcaccatctggggaccctgacgtgaagtttgtcactgactggagaactcggtttgtcgtctgttgcgggggcggcagttatggcggtgccgttgggcagtgcacccgtacctttgggagcgcgcgccctcgtcgtgtcgtgacgtcacccgttctgttggcttataatgcagggtggggccacctgccggtaggtgtgcggtaggcttttctccgtcgcaggacgcagggttcgggcctagggtaggctctcctgaatcgacaggcgccggacctctggtgaggggagggataagtgaggcgtcagtttctctggtcggttttatgtacctatcttcttaagtagctgaagctccggttttgaactatgcgctcggggttggcgagtgtgttttgtgaagttttttaggcaccttttgaaatgtaatcatttgggtc'.upper()
chimera = 'tctcaaaagcgggcatgacttctgcgctaagattgtcagtttccaaaaacgaggaggatttgatattcacctggcccgcggtgatgcctttgagggtggccgcgtccatctggtcagaaaagacaatctttttgttgtcaagcttgaggtgtggcaggct'.upper()


###################################

# constant regions for secondary structure calculation

flank5utr = 'TAGTTCCGTCGCAGCCGGGATTTGGGTCGCAGTTCTTGTTTGTGGATCGCTGTGATCGTCGTCACcaGgtgagtagcggg'.upper()
flank3toBC = 'gcctcagacagtggttcaaagtttttttcttccatttcagGTGTCGTGAA'.upper()
BCflank3 = 'ACgtCaCCTAaggcGCGCCAAA'.upper()
splicedPreBC = 'TAGTTCCGTCGCAGCCGGGATTTGGGTCGCAGTTCTTGTTTGTGGATCGCTGTGATCGTCGTCACcaGGTGTCGTGAA'.upper()
dTpreBC = 'TAGTTCCGTCGCAGCCGGGATTTGGGTCGCAGTTCTTGTTTGTGGATCGCTGgctagctcgtGAttcagGTGTCGTGAA'.upper()
dTpostBC = 'ACgtCaCATAaacaGCGCCAAA'.upper()
gORF60 = 'atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggac'.upper()
dORF60 = 'atggtgagcaagggcgaggaggtcatcaaagagttcatgcgcttcaaggtgcgcatggag'.upper()


###################################

# sample indexes in primers (as they would be read in the fwd HiSeq read)

sample_id = {'GATC':0,
             'AGCT':1,
             'GTAT':2}
sample_id.update({OO:0 for OO in one_offs('GATC')})
sample_id.update({OO:1 for OO in one_offs('AGCT')})
sample_id.update({OO:2 for OO in one_offs('GTAT')})

pool_id = sample_id

# sample indexes of sort samples only

sort_id = { 'CCAA':0,
            'TTGG':1,
            'GGCA':2,
            'AACC':3,
            'CTTG':4,
            'TATA':5,
            'ATGC':0,
            'CATT':1,
            'TCCC':2,
            'GAGA':3,
            'CGAT':4,
            'ACTC':5}
sort_id1 = {'CCAA':0,
            'TTGG':1,
            'GGCA':2,
            'AACC':3,
            'CTTG':4,
            'TATA':5}
sort_id2 = {'ATGC':0,
            'CATT':1,
            'TCCC':2,
            'GAGA':3,
            'CGAT':4,
            'ACTC':5}
sort_sample_id = {'GATC':'RILwc1',
             'AGCT':'RILwc2',
             'GTAT':'RILwc3',
             'CCAA':'dTom++1',
             'TTGG':'dTom+1',
             'GGCA':'dTom1',
             'AACC':'GFP1',
             'CTTG':'GFP+1',
             'TATA':'GFP++1',
             'ATGC':'dTom++2',
             'CATT':'dTom+2',
             'TCCC':'dTom2',
             'GAGA':'GFP2',
             'CGAT':'GFP+2',
             'ACTC':'GFP++2'}

###################################

# sample indexes in all primers

primer_id = {'EK287':'GATC',
             'EK288':'AGCT',
             'EK289':'GTAT',
             'EK302':'GATC',
             'EK303':'AGCT',
             'EK304':'GTAT'}