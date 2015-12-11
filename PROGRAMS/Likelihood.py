'''
Dan Adler, Farhan Damani
Computational Genomics Final

Likelihood
'''

import Read_Mapping as rm
import math

def _likelihood(reads_dict, contig_list):
    '''
        Compute total data log likelihood

        likelihood = log p(x,y | s,o,b) + log p(b) + log p(o|s,b) + log p(s)
    '''
    sum = 0
    for r in reads_dict:
        for i in reads_dict[r]:
            p = float(i[2])
            sum = sum + p
    return sum