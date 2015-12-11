'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for contig merging.

Methods:
1) merge_check_global()
2) suffixPrefixMatch()
3) suffix_filter()
4) extract_bbl()
5) trace_contigs()
6) merge_contigs()
7) reverse_reads_dict()
8) change_reads_on_merge()
9) run_merge() which runs everything above

Essentially the basic idea is that we are trying to merge contigs using an overlap graph
and then map reads to new locations.

'''

# IMPORTS
from collections import Counter

'''
Checks all nchoosek(contig_length,2) contigs using merge_check_local() to see which contigs
should be merged using merge_contigs()

@param contigs_list is the list of all current contigs
@param OVERLAP_LENGTH (defaults to 15)

@return the overlaps or a blank list
'''
def merge_check_global(contig_list, OVERLAP_LENGTH):
    # Merge will be 0, and a 1 is added each time merge_local() comes back as a needed merge
    merge = 0
    # Dictionary to hold permutations
    contig_dict = {}
    # Go through each contig
    for suffix in contig_list:
        # Check if a merge is needed
        contig_dict[suffix] = Counter()
        for prefix in contig_list:
            # Go through all prefixes for each suffix
            if prefix != suffix:
                contig_dict[suffix][prefix] = suffixPrefixMatch(suffix, prefix, OVERLAP_LENGTH)
                # If every overlap comes out to 0, this parameter remains at 0 (i.e. no merges need to occur)
                merge += contig_dict[suffix][prefix]

    # Check if merge_contigs needs to occur
    if merge > 0:
        return contig_dict
    else:
        return []

'''
Finds the longest suffix that is a prefix of
another string.
Code was taken from Ben's github (also used on HW4)
http://nbviewer.ipython.org/github/BenLangmead/comp-genomics-class/blob/master/projects/UnpairedAsmChallenge.ipynb

@param str1 is the string whose suffix will be checked
@param str2 is the string whose prefix will be checked
@param min_overlap is the minimum overlap length
@param OVERLAP_LENGTH (defaults to 15)
@returns the length of the match
'''
def suffixPrefixMatch(str1, str2, min_overlap):
    ''' Returns length of longest suffix of str1 that is prefix of
        str2, as long as that suffix is at least as long as min_overlap. '''
    #print str2
    if len(str2) < min_overlap: return 0
    str2_prefix = str2[:min_overlap]
    str1_pos = -1
    while True:
        str1_pos = str1.find(str2_prefix, str1_pos + 1)
        if str1_pos == -1: return 0
        str1_suffix = str1[str1_pos:]
        if str2.startswith(str1_suffix): return len(str1_suffix)

'''
Suffix filter (i.e. choose greatest prefix match for a suffix)

@param contigs_dict is a dictionary for each contig as a suffix that has a counter object
@return each contig in a dictionary with their respective right best buddy
'''
def suffix_filter(contig_dictionary):
    bbr = {}

    # Go through each suffix
    for suffix in contig_dictionary:
        # Check if there exists a match
        if len(contig_dictionary[suffix]) > 0:
            # Get most common match
            c = contig_dictionary[suffix].most_common(1)
            if c[0][1] > 0:
                bbr[suffix] = [c[0][0], c[0][1]]
    return bbr

'''
Go through best buddy list.
If entry exists as bbr twice, we take the one with the higher score.

Delete the other one.
@param bb is the best buddy list.
@return a bb list with only relevant terms
'''
def extract_bbl(bbr):
    bbl_set = set()
    bbl = {}
    # Find all possible left best buddies
    for p in bbr:
        bbl_set.add(bbr[p][0])
    # Go through each one and take its right buddy
    for l in bbl_set:
        c = Counter()
        for p in bbr:
            if bbr[p][0] == l:
                c[p] = bbr[p][1]
        # Get most common and see if the bbl is unique, keep if it is.
        m = c.most_common(2)
        if len(m) > 1 and m[0][1] == m[1][1]:
            pass
        else:
            bbl[l] = [m[0][0],m[0][1]]
    return bbl

'''
Trace the new contigs.

@param bbl is bbl dictonary
@return contigs, the contig traces
'''
def trace_contigs(bbl):
    right = set()
    left = set(bbl.keys())
    # Figure out ending points
    for p in bbl:
        right.add(bbl[p][0])
    end = left.difference(right)
    # Traceback
    contigs = [None]*len(end)
    i = 0
    for e in end:
        curr = e
        contigs[i] = []
        while curr in bbl:
            contigs[i].insert(0,[curr, bbl[curr][1]])
            curr = bbl[curr][0]
        contigs[i].insert(0,[curr, 0])
        i += 1

    return list(contigs)

'''
Perform merges

@param the contigs_list with overlap
@param the original contigs
@return the new contigs
'''
def merge_contigs(contigs_list, contigs):
    new_contigs = []
    for c in contigs_list:
        new = ''
        # Each list element is [contig, prefix-overlap with previous]
        for j in c:
            if contigs is not None and j[0] in contigs:
                contigs.remove(j[0])
            new += j[0][j[1]:]
        new_contigs.append(new)
    # Take unmerged/touched contigs and append them onto the list
    if contigs is not None:
        for c in contigs:
            new_contigs.append(c)

    return new_contigs

'''
Reverse reads dictionary

@param reads_dict is the reads dictionary
@return the reversed dictionary with {contig: [[r_1,o_1],...]}
'''
def reverse_reads_dict(reads_dict):
    reverse_dict = {}
    # Go through each read
    for r in reads_dict:
        # Take contig and make it key, then read gets added as a value
        for pair in reads_dict[r]:
            if pair[0] not in reverse_dict:
                reverse_dict[pair[0]] = []
            reverse_dict[pair[0]].append([r,pair[1],pair[2],pair[3]])

    return reverse_dict

'''
Change position of reads during a merge

@param reverse_dict is the dictionary of reads (assuming {read: [[s, o]]})
@param reads_dict is the original dictionary
@param contigs_list are the list of contigs to be joined and their overlap
@param contigs are the original contigs list
@return new_reads_dict once updated
'''
def change_reads_on_merge(reverse_dict, reads_dict, contigs_list, contigs, new_contigs):
    new_reads_dict = {}
    read_length = len(contigs[0])

    # Go through the contig traces
    for s in range(len(contigs_list)):
        o = 0
        # Find the trace values
        for i in range(len(contigs_list[s])):
            # Find reads that match up with that value and add them to the new dictionary
            for reads in reverse_dict[contigs.index(contigs_list[s][i][0])]:
                if reads[0] not in new_reads_dict:
                    new_reads_dict[reads[0]] = []
                if [s,o + reads[1],reads[2],reads[3]] not in new_reads_dict[reads[0]]:
                    new_reads_dict[reads[0]].append([s,o + reads[1],reads[2],reads[3]])
            if i < len(contigs_list[s]) - 1:
                o += (read_length - contigs_list[s][i+1][1])
    # Get reads that have not yet been added
    for i in reads_dict:
        if i not in new_reads_dict:
            new_reads_dict[i] = []
            for j in reads_dict[i]:
                new_reads_dict[i].append([new_contigs.index(contigs[j[0]]),j[1],j[2],j[3]])
            
    return new_reads_dict

'''
Run a merge contigs with default overlap_length of 15.

@param contigs is the contigs list
@param reads is the reads list
@param OVERLAP_LENGTH is the minimum length of an overlap we will accept
@return the new contigs list and reads list
'''
def run_merge(contigs, reads, OVERLAP_LENGTH = 15):
    contig_dict = merge_check_global(contigs, OVERLAP_LENGTH)
    if contig_dict == []:
        return contigs, reads
    bbr = suffix_filter(contig_dict)
    bbl = extract_bbl(bbr)
    contig_trace = trace_contigs(bbl)
    new_contigs = merge_contigs(contig_trace, contigs[:])
    reverse_reads = reverse_reads_dict(reads)
    new_reads = change_reads_on_merge(reverse_reads, reads, contig_trace, contigs, new_contigs)
    return new_contigs, new_reads



