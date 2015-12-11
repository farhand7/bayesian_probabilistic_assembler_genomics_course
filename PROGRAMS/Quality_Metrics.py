__author__ = 'farhan_damani'

import matplotlib.pyplot as plt
import numpy as np


def _create_likelihood_figure(ll_list):
    iters = []
    for i in range(len(ll_list)): iters.append(i)
    plt.plot(np.array(iters),np.array(ll_list), label ='Likelhood')
    plt.xlabel('Iteration')
    plt.ylabel('-log(likelihood)')
    plt.show()


def _process_likelihood_output_file(f):
    ll = []
    for i,l in enumerate(f):
        if i == 0 or i == 1: continue # skip ll = 0
        w = l.split('\t')
        if len(w) > 0:
            ll.append(float(w[1]))
    return ll



def _percent_change_of_contigs():
    '''
        Evaluate percent chance of each contig over time steps.
        Use edit distance to show this.

    '''

    f = open('../SRC_OUTPUT/trial_three/contig.txt', 'r')
    c_prev = ''
    contig_one_list = []
    contig_two_list = []
    contig_three_list = []
    for l in f:
        w = l.split('\t')
        if len(w) > 0:
           if str(w[1].strip())=='merge': continue
           if str(w[1].strip()) =='end': continue
           contig_one_list.append(str(w[2].strip())), contig_two_list.append(str(w[3].strip())), contig_three_list.append(str(w[4].strip()))

    for i in xrange(0,len(contig_one_list)):
        if i == 0: continue
        print edDistDp(contig_one_list[i], contig_one_list[i-1]), edDistDp(contig_two_list[i], contig_two_list[i-1]), \
            edDistDp(contig_three_list[i], contig_three_list[i-1])

def edDistDp(x, y):
    """ Calculate edit distance between sequences x and y using
        matrix dynamic programming.  Return distance.

    citation: http://nbviewer.ipython.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_DP_EditDist.ipynb

    """
    D = np.zeros((len(x)+1, len(y)+1), dtype=int)
    D[0, 1:] = range(1, len(y)+1)
    D[1:, 0] = range(1, len(x)+1)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            delt = 1 if x[i-1] != y[j-1] else 0
            D[i, j] = min(D[i-1, j-1]+delt, D[i-1, j]+1, D[i, j-1]+1)
    return D[len(x), len(y)]

def _create_change_in_contigs_figure(f):
    c_1,c_2,c_3,iters= [],[],[],[]
    for i,l in enumerate(f):
        w = l.split(' ')
        if len(w) > 0:
            c_1.append(int(w[0])), c_2.append(int(w[1])), c_3.append(int(w[2]))
            iters.append(i)
    plt.plot(np.array(iters),np.array(c_1), label ='contig one')
    plt.plot(np.array(iters),np.array(c_2), label ='contig two')
    plt.plot(np.array(iters),np.array(c_3), label ='contig three')
    plt.xlabel('Iteration')
    plt.ylabel('Edit Distance')
    plt.legend(loc='lower left')
    plt.show()


def create_scatterplot():
    '''
        Create a boxplot of k-mer positions at contig offsets
    '''
    f = open('../SRC_OUTPUT/Dan_trial_two/Dan_trial_two_reads_0.txt', 'r') # random data initialized
    f2 = open('../SRC_OUTPUT/Dan_trial_two/Dan_trial_two_reads_20.txt', 'r') # after 20 iters

    labels_0,labels_20, data_0,data_20 = [],[],[],[]
    for i,l in enumerate(f):
        if i == 0: continue
        w = l.split(' ')
        if len(w) > 0:
            label = w[0]+','+w[1].split('-')[0]

            #label = w[0]+','+w[1][:-1]
            data = []
            for n in w[2:]: data.append(int(n))
            labels_0.append(label), data_0.append(data)
    for i,l in enumerate(f2):
        if i == 0: continue
        w = l.split(' ')
        if len(w) > 0:
            #label = w[0]+','+w[1][:-1]
            label = w[0]+','+w[1].split('-')[0]

            data = []
            if len(w) > 3:
                for n in w[2:]: data.append(int(n))
            else:
                data.append([])
            labels_20.append(label), data_20.append(data)



    plt.figure()
    plt.boxplot(data_0)
    plt.xlabel('Contig, Offset Bin')
    plt.ylabel('K-mer Offset')
    iters = []
    for i in range(len(data_0)): iters.append(i)


    plt.xticks(iters,labels_0, rotation='45')

    plt.show()


    plt.figure()
    plt.boxplot(data_20)
    plt.xlabel('Contig, Offset Bin')
    plt.ylabel('K-mer Offset')
    plt.xticks(iters,labels_20, rotation='45')

    plt.show()


def _driver():
    #f = open('../SRC_OUTPUT/trial_one/likelihood.txt', 'r')
    #_create_likelihood_figure(_process_likelihood_output_file(f))
    #_percent_change_of_contigs()
    #f = open('../figures/change_in_contigs_data.txt', 'r')
    #_create_change_in_contigs_figure(f)
    create_scatterplot()
