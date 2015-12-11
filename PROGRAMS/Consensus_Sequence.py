'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for consensus sequencing

Goes through all reads and creates a counter for reads in a certain offset with most frequent position
Returns new contig list based upon that.

Methods:
1. consensus_sequence()
2. compute_new_contigs()
3. run_consensus()

USES the following dictionary data structure (dictionary within dictionary)
Could imagine it's like a tree
{ CONTIG : {OFFSET : Counter()} }
'''

# IMPORTS
from collections import Counter

'''
Creates data structure object and iterates through each read to 
get a consensus sequence information.

@param reads_dict is a current dictionary with {read : [[s1, o1],...]}
@return the frequency information for each contig
'''
def consensus_sequence(reads_dict):
    # Object
    frequency_info = {}

    # For each read
    for r in reads_dict:
        for l in reads_dict[r]:
            s = l[0]
            o = l[1]
            # Check if contig number is in structure
            if s not in frequency_info:
                frequency_info[s] = {}
            for nt in r:
                # Check if offset in the contig number dictionary
                if o not in frequency_info[s]:
                    frequency_info[s][o] = Counter()
                # Add counter to that nucleotide in the contig at an offset
                frequency_info[s][o][nt] += 1
                o += 1

    # Send dictionary to compute the new contigs
    return frequency_info

'''
Compute new contigs.

@param the frequency_info structure with structure { CONTIG : {OFFSET : Counter()} }
@return the new contig sequence
'''
def compute_new_contigs(frequency_info):

    # New list
    new_contigs = []

    # Get contigs
    contigs = frequency_info.keys()
    contigs.sort()

    # For each contig
    for s in contigs:
        # Start new contig
        curr = ''
        offsets = frequency_info[s].keys()
        offsets.sort()
        for o in offsets:
            # Get most common nucleotide at certain (s,o)
            nt = frequency_info[s][o].most_common(1)
            curr += nt[0][0]

        new_contigs.append(curr)

    return new_contigs


'''
Main function to compute new consensus sequence

@param reads_dict is the dictionary of reads with mappings to [s,o] placements
@return the new contig list
'''
def run_consensus(reads):
    frequency_info = consensus_sequence(reads)
    return compute_new_contigs(frequency_info)
