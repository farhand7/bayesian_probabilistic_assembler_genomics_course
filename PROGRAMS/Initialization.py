'''
Dan Adler, Farhan Damani
Computational Genomics Final

Initialization.py

Methods:

1) _process()
2) _read_fa_file()

The idea is that we are going to take our reads, make a reads dict and then randomly
map them to contigs with a specified length and number.  The original number of contigs are
predetermined by the parameters below.
'''

# IMPORTS
from random import randint
import math

# GLOBAL VARIABLES
CONTIG_LENGTH = 1000
NUM_CONTIGS = 3
DEFAULT_PROB = 1*10**(-100)
K_MER = 100


'''
reads dictionary maps read -> 
[ [s_1,o_1,p_1,offset],[s_2,o_1,p_2,offest],...[s_n,o_n,p_n,offset] ]
where n is the number of mappings of different reads 
from different viruses but are the same sequence

@param f: file of reads
@return: reads_dict
'''
def _process(f):
    reads = _read_fa_file(f)
    # map read -> [s,o]
    reads_dict = {}
    for r in reads:
        if r[0] not in reads_dict:
            reads_dict[r[0]] = []
        l = []
        # randomly assign read to one of the 7 contigs
        s = randint(0,NUM_CONTIGS - 1)
        o = randint(0,CONTIG_LENGTH-len(r)+1)
        l.append(s),l.append(o),l.append(-math.log(DEFAULT_PROB)),l.append(r[1])
        reads_dict[r[0]].append(l)
    return reads_dict

'''
Read fa file and return list of reads

@param f is the filename
@return the reads in a list
'''
def _read_fa_file(f):
    
    reads = []

    for i,l in enumerate(f):
        if l[0] == '>':
            curr = l[1:-1]
        w = l.strip().split()
        if len(w) > 0:
            if i%2 == 1: reads.append([''.join(w),curr])
    return reads
