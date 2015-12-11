'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for read mapping

Maps each read to a contig, offset in contig, and the optimal alignment of read to the contig.
'''

# IMPORTS
import math
import random
import Initialization as init

# GLOBAL VARIABLES 
geometric_prob = 1.0/3.0

'''
Compute p(s,o,y | . )
We simplify this computation to
p(s,o|y*,.) = p(s_i = s | s_-i) * p(o_i = o | s_i = s) * p(x_i, y*_so | s_i = s, o_i = o, b)

@param x_i,s_i,o_i,y_i are the random variables
@return the approximate joint probability (or really -log(p))
'''
def _compute_p_s_o_y(x_i,s_i,o_i,y_i):
    a =  _compute_p_s(s_i)+_compute_p_o(o_i)+_compute_p_x_y(x_i,y_i)
    return a

'''
Compute p(s_i = s | s_-i)
P(s) is modeled as a geometric distribution
P(s) = (1-p)^(k-1) * p -> this is the probability that the first occurrence of success
requires k number of independent trials, each with success probability p.

@param s_i : contig for read i
@return: probability of read i belonging to contig s_i (or -log(p))
'''
def _compute_p_s(s_i):
    # geometric distribution
    return -math.log((1-geometric_prob)**(s_i)*geometric_prob)


'''
Compute p(o_i = o | s_i = s)

@param o_i is the offset
@return the probability of that offset (or -log(p))
'''
def _compute_p_o(o_i):
    # uniform distribution
    return -math.log(1.0/(init.CONTIG_LENGTH-init.K_MER)) 

'''
Compute p(x,y*_so | s_i = s, o_i = 0, b)
x_i, y_i ~ A(s_i, o_i, b, p_mis)
SCORE(read_i) = log(0.5) + n_hit * log(1-p_mis) + n_mis * log(p_mis / (|B|-1))

@param x: read i
@param y: optimal alignment of x to contigs
@return: the probability of that alignment (or -log(p))
'''
def _compute_p_x_y(x,y):
    n_hit = 0
    p_mis = 0.1
    for i in range(len(x)):
        if x[i] == y[i]: n_hit = n_hit + 1
    # Binomial like distribution
    return -(math.log(0.5)+n_hit*math.log(1-p_mis) + (len(x)-n_hit)*math.log(p_mis / (3)))

'''
Compute y* = arg min (hamming distance) over all possible (s,o)

@param x: read i
@param contigs: list of contigs
@return y* sequence
'''
def _compute_y_star(x, contigs):
    min_hamming = len(x) + 1
    y_star = ''
    s_star = -1
    o_star = -1
    for s,c in enumerate(contigs):
        # from 0 to n-k+1
        for i in xrange(0,len(c)-len(x)+1):
            # x and substring starting at index i up to length of x
            hamming_dist = _compute_hamming(x,c[i:i+len(x)])
            if hamming_dist < min_hamming:
                min_hamming = hamming_dist
                y_star = c[i:i+len(x)]
                s_star = s
                o_star = i
    return y_star, s_star, o_star

'''
Compute hamming distance of x and y

@param x,y : two sequences
@return hamming score
'''
def _compute_hamming(x,y):
    score = 0
    for i in range(len(x)):
        if x[i] != y[i]: score = score + 1
    return score


'''
Main method for this Library

@param reads_dict is a dictionary of reads
@param contigs is the current contigs list based upon these reads
@return the new reads dictionary
'''
def run(reads_dict, contigs):
    # Number of reads we will sampled
    percent_sampling = .6
    key_list = reads_dict.keys()
    order = []
    # Pick a random 60% (without replacement)
    for i in range(len(key_list)):
        if random.random() > percent_sampling:
            order.append(i)

    for i in order:
        candidate_probabilities = {}
        # take a random read, x
        x = key_list[i]
        # take random mapping from that read (if it maps to > 1 spot)
        if len(reads_dict.get(x)) > 1:
            read_index = random.randint(0,len(reads_dict.get(x))-1)
        else:
            read_index = 0

        '''GET Y_STAR'''
        # Get the list value
        read_mapping_list = reads_dict.get(x)[read_index]
        # Get current s,o,p_prev
        s,o,p_prev = read_mapping_list[0],read_mapping_list[1],read_mapping_list[2]
        # Compute new s,o at optimal alignment
        y_star,s_star,o_star = _compute_y_star(x,contigs)
        max_prob = _compute_p_s_o_y(x,s_star,o_star,y_star)

        '''APPROXIMATE INFERENCE'''
        u = random.random()
        max_prob = float(math.exp(-max_prob))
        p_prev = float(math.exp(-p_prev))
        alpha = min(1,float(max_prob/p_prev))
        # accept
        if u < alpha:
            reads_dict[x][read_index][0],reads_dict[x][read_index][1],reads_dict[x][read_index][2] = s_star, o_star, -math.log(max_prob)
    
    return reads_dict