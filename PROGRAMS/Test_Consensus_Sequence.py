'''
Dan Adler, Farhan Damani
Computational Genomics Final

Tests for Consensus Sequencing.
'''

import unittest
import Consensus_Sequence as cs 

class Test_Consensus_Seuqence(unittest.TestCase):

    # Test with one contig
    def test_single_contig_1(self):
        contigs_1 = ['ATCGCTGATT']
        reads_1 = {\
        'ATCG':[[0,0,0.5]],'TCGC':[[0,1,0.5]],'CCCT':[[0,2,0.5]],'GCTG':[[0,3,0.5]],\
        'TTCA':[[0,4,0.5]],'TGAT':[[0,5,0.5]],'GATT':[[0,6,0.5]]}
        freq = cs.consensus_sequence(reads_1)
        self.assertEqual(contigs_1,cs.compute_new_contigs(freq))

    # Test with single contig and equal probabilities of certain nts
    def test_single_contig_2(self):
        contigs_2 = ['ATCGCTGATT','ATCCCTCATT','ATCGCTCATT','ATCCCTGATT']
        reads_2 = {\
        'ATCC':[[0,0,0.5]],'TCGC':[[0,1,0.5]],'CCCT':[[0,2,0.5]],'GCTC':[[0,3,0.5]],\
        'TTCA':[[0,4,0.5]],'TGAT':[[0,5,0.5]],'GATT':[[0,6,0.5]]}
        freq = cs.consensus_sequence(reads_2)
        self.assertTrue(cs.compute_new_contigs(freq)[0] in contigs_2)

    # Test with multiple contigs
    def test_multiple_contigs_1(self):
        contigs_3 = ['ATCGCTGATT','TTTACGATGC']
        reads_3 = {\
        'ATCG':[[0,0,0.5]],'TCGC':[[0,1,0.5]],'CCCT':[[0,2,0.5]],'GCTG':[[0,3,0.5]],\
        'TTCA':[[0,4,0.5]],'TGAT':[[0,5,0.5]],'GATT':[[0,6,0.5]],\
        'TTAA':[[1,0,0.5]],'TTAC':[[1,1,0.5]],'TGCG':[[1,2,0.5]],'ACGA':[[1,3,0.5]],\
        'AGAT':[[1,4,0.5]],'GTTG':[[1,5,0.5]],'AAGC':[[1,6,0.5]]}
        freq = cs.consensus_sequence(reads_3)
        self.assertEqual(contigs_3,cs.compute_new_contigs(freq))

    # Test with multiple contigs and equal probabilities of certain nts
    def test_multiple_contigs_2(self):
        contigs_4 = ['ATCGCTGATT','ATCCCTCATT','ATCGCTCATT','ATCCCTGATT',\
            'TTTACGATGC','TTTAAGATGC']
        reads_4 = {\
        'ATCG':[[0,0,0.5]],'TCGC':[[0,1,0.5]],'CCCT':[[0,2,0.5]],'GCTG':[[0,3,0.5]],\
        'TTCA':[[0,4,0.5]],'TGAT':[[0,5,0.5]],'GATT':[[0,6,0.5]],\
        'TTAA':[[1,0,0.5],[1,1,0.5]],'TGCG':[[1,2,0.5]],'CCGA':[[1,3,0.5]],\
        'AGAT':[[1,4,0.5]],'GTTG':[[1,5,0.5]],'AAGC':[[1,6,0.5]]}
        freq = cs.consensus_sequence(reads_4)
        self.assertTrue(cs.compute_new_contigs(freq)[0] in contigs_4\
            and cs.compute_new_contigs(freq)[1] in contigs_4)

if __name__ == '__main__':
    unittest.main()


