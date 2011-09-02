import unittest
import collections

import sys
import itertools

sys.path.append("../scripts")
from flatfeature import Bed
#from pseudo import 
from mask_non_cds import mask_non_cds
from merge import merge_same_hits

class TestPseudo(unittest.TestCase):
    def setUp(self):
        self.sbed = Bed("/Users/gturco/rice_v6.bed", "/Users/gturco/rice_v6.fasta") ;self.sbed.fill_dict()

    def test_mask_non_cds(self):
        """test for test_mask_non_cds TODO mask inbetween genes currently does
        not mask inbetween genes"""
        handle = open("/Users/gturco/1.fasta")
        fh = handle.read()
        # double checked these are starts and ends acorrding to coge and match
        # with other fasta
        self.assertEqual(fh[71774],'N')
        self.assertEqual(fh[78937],'N')
        #first start
        self.assertEqual(fh[71902:71961],'ATGGCGGCTCCCAAGCCCATCTGGGTGCGCCAGGCCGAGGAGGCCAAGCTCAAGTCCGA')
        #first end
        self.assertEqual(fh[72868:72935],'CAGGCACTGCAGGCACATGCTGCGCAGATGCAGGCTGATTCCAAGGCTGCAGGAGGGGAGGCTTCAG')
        #second one start
        self.assertEqual(fh[73467:73509],'GTTCAACAGATAAAACTGATAAAGGTGATGTTCTGAAGAAAA')
        #third one end
        self.assertEqual(fh[75967:76008],'TTGATGGTGCATTTGTGGGCACTGAAGAGTCTGCTATATGA')

    def test_merge_same_hits(self):
        """confrim that genes are being grouped correctly"""
        new_bed = Bed('/Users/gturco/rice_v6.all.bed')
        self.assertEqual(new_bed.accn('rice_v6_3_13648324_13648522')['locs'][0][0],13648324)
        self.assertEqual(len(new_bed.accn('rice_v6_3_13648324_13648522')['locs']),1)
        self.assertEqual(len(new_bed.accn('rice_v6_1_19308831_19309020')['locs']),3)

if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
