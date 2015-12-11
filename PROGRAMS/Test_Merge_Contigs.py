'''
Dan Adler, Farhan Damani
Computational Genomics Final

Tests for the Merge_Contigs.py file.
'''

import unittest
import Merge_Contigs as mc
from collections import Counter

class Test_Merge_Contigs(unittest.TestCase):

    # Test whether two contigs are merged
    def test_merge_check_1(self):
        c_1 = ['GATCTTTT', 'GATCGATC']
        test_1 = mc.merge_check_global(c_1, 3)
        truth_1 = {'GATCTTTT': Counter({'GATCGATC':0}), 'GATCGATC': Counter({'GATCTTTT':4})}
        self.assertEqual(test_1, truth_1)

    # Test whether no contigs are merged (overlap not big enough)
    def test_merge_check_2(self):
        c_1 = ['GATCTTTT', 'GATCGATC']
        test_1 = mc.merge_check_global(c_1, 5)
        truth_1 = []
        self.assertEqual(test_1, truth_1)

    # Test three contigs where two have overlap, but one has more overlap than the other
    def test_three_contig_1(self):
        c_1 = ['GATCTTTT', 'GATCGATC','ATCAAAAA']
        contig_dict_1 = mc.merge_check_global(c_1, 3)
        bbr_test_1 = mc.suffix_filter(contig_dict_1)
        bbr_truth_1 = {'GATCGATC':['GATCTTTT', 4]}
        self.assertEqual(bbr_test_1,bbr_truth_1)
        c_2 = ['GATCAAAA', 'CCCGGAAA','AAAAGATC']
        contig_dict_2 = mc.merge_check_global(c_2, 3)
        bbr_test_2 = mc.suffix_filter(contig_dict_2)
        bbr_truth_2 = {'CCCGGAAA':['AAAAGATC', 3], 'AAAAGATC':['GATCAAAA', 4], 'GATCAAAA':['AAAAGATC', 4]}
        self.assertEqual(bbr_test_2,bbr_truth_2)
        

    # Test three contigs where one best buddy should be eliminated when looking left
    def test_three_contig_2(self):
        c_1 = ['GATCAAAA', 'CCCGGAAA','AAAAGATC']
        contig_dict_1 = mc.merge_check_global(c_1, 3)
        bbr_test_1 = mc.suffix_filter(contig_dict_1)
        bbl_test_1 = mc.extract_bbl(bbr_test_1)
        bbl_truth_1 = {'AAAAGATC':['GATCAAAA', 4], 'GATCAAAA':['AAAAGATC',4]}
        self.assertEqual(bbl_test_1, bbl_truth_1)

    # Test three contigs traced
    def test_three_contig_3(self):
        c_1 = ['GATCTAAA', 'CCCGAAAA','AAAAGATC']
        contig_dict_1 = mc.merge_check_global(c_1, 3)
        bbr_test_1 = mc.suffix_filter(contig_dict_1)
        bbl_test_1 = mc.extract_bbl(bbr_test_1)
        contig_trace_test_1 = mc.trace_contigs(bbl_test_1)
        contig_trace_truth_1 = [[['CCCGAAAA',0],['AAAAGATC',4],['GATCTAAA',4]]]
        self.assertEqual(contig_trace_test_1, contig_trace_truth_1)

    # Test three contigs that should be merged
    def test_three_contig_4(self):
        c_1 = ['GATCTAAA', 'CCCGAAAA','AAAAGATC']
        contig_dict_1 = mc.merge_check_global(c_1, 3)
        bbr_test_1 = mc.suffix_filter(contig_dict_1)
        bbl_test_1 = mc.extract_bbl(bbr_test_1)
        contig_trace_test_1 = mc.trace_contigs(bbl_test_1)
        new_contigs_test_1 = mc.merge_contigs(contig_trace_test_1, c_1)
        new_contigs_truth_1 = ['CCCGAAAAGATCTAAA']
        self.assertEqual(new_contigs_test_1, new_contigs_truth_1)

    # Test reverse reads
    def test_read_reverse_1(self):
        reads_dict_1 = {'AAAAT':[[0,0,0.5,1]],'AAATC':[[0,1,0.5,1]],'AATCA':[[0,2,0.5,1],[1,0,0.5,1]],'ATCAG':[[1,1,0.5,1]],'TTTTT':[[2,0,0.5,1]]}
        reverse_dict_test_1 = mc.reverse_reads_dict(reads_dict_1)
        reverse_dict_truth_1 = {0:[['AAAAT',0,0.5,1],['AAATC',1,0.5,1],['AATCA',2,0.5,1]],1:[['AATCA',0,0.5,1],['ATCAG',1,0.5,1]],2:[['TTTTT',0,0.5,1]]}
        for i in reverse_dict_test_1:
            for j in reverse_dict_test_1[i]:
                self.assertTrue(j in reverse_dict_truth_1[i])

    # Test reverse reads
    def test_read_updates_1(self):
        c_1 = ['AAAATCA', 'AATCAGG', 'TTTTTTT']
        contig_dict_1 = mc.merge_check_global(c_1,3)
        bbr_test_1 = mc.suffix_filter(contig_dict_1)
        bbl_test_1 = mc.extract_bbl(bbr_test_1)
        contig_trace_test_1 = mc.trace_contigs(bbl_test_1)
        new_contigs_test_1 = mc.merge_contigs(contig_trace_test_1, c_1[:])
        reads_dict_1 = {'AAAAT':[[0,0,0.5,1]],'AAATC':[[0,1,0.5,1]],'AATCA':[[0,2,0.5,1],[1,0,0.5,1]],'ATCAG':[[1,1,0.5,1]],'TTTTT':[[2,0,0.5,1]]}
        reverse_dict_test_1 = mc.reverse_reads_dict(reads_dict_1)
        new_reads_dict_test_1 = mc.change_reads_on_merge(reverse_dict_test_1, reads_dict_1, contig_trace_test_1, c_1, new_contigs_test_1)
        new_reads_dict_truth_1 = {'ATCAG': [[0,3,0.5,1]], 'AATCA': [[0,2,0.5,1]], 'AAAAT': [[0,0,0.5,1]], 'TTTTT': [[1,0,0.5,1]], 'AAATC': [[0,1,0.5,1]]}
        self.assertEqual(new_reads_dict_truth_1, new_reads_dict_test_1)

    # Test Run through
    def test_run(self):
        c_1 = ['AAAATCA', 'AATCAGG', 'TTTTTTT']
        reads_dict_1 = {'AAAAT':[[0,0,0.5,1]],'AAATC':[[0,1,0.5,1]],'AATCA':[[0,2,0.5,1],[1,0,0.5,1]],'ATCAG':[[1,1,0.5,1]],'TTTTT':[[2,0,0.5,1]]}
        new_contigs_test_1, new_reads_test_1 = mc.run_merge(c_1,reads_dict_1,3)
        new_contigs_truth_1 = ['AAAATCAGG','TTTTTTT']
        new_reads_dict_truth_1 = {'ATCAG': [[0,3,0.5,1]], 'AATCA': [[0,2,0.5,1]], 'AAAAT': [[0,0,0.5,1]], 'TTTTT': [[1,0,0.5,1]], 'AAATC': [[0,1,0.5,1]]}
        self.assertEqual(new_reads_dict_truth_1, new_reads_test_1)
        self.assertEqual(new_contigs_truth_1,new_contigs_test_1)

if __name__ == '__main__':
    unittest.main()