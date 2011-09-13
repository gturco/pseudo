import unittest
import collections

import sys
import itertools

sys.path.append("../scripts")
from flatfeature import Bed
#from pseudo import 
from mask_non_cds import mask_non_cds
from merge import merge_same_hits, merge
from pseudo import group_cds, best_hit

class TestPseudo(unittest.TestCase):
    def setUp(self):
        self.qallbed = Bed('data/rice_v6_setaria64/rice_v6.all.bed')
        self.sallbed = Bed('data/rice_v6_setaria64/setaria64.all.bed')
        self.sbed = Bed("data/rice_v6_setaria64/rice_v6.bed", "data/rice_v6_setaria64/rice_v6.fasta") ;self.sbed.fill_dict()
        self.missed_bed = Bed("data/rice_v6_setaria64/missed_rice_v6_from_setaria64.bed") 
        self.matches = "data/rice_v6_setaria64/missed_rice_v6_from_setaria64.matches.txt"
    #def test_mask_non_cds(self):
    #    """test for test_mask_non_cds TODO mask inbetween genes currently does
    #    not mask inbetween genes"""
    #    for seqid, seq in mask_non_cds(self.sbed,'N', None):
    #        f = ("data/rice_v6_setaria64/{0}.fasta".format(seqid))
    #        fh = open(f, "wb")
    #        print >> fh,seq
    #        fh.close()
    #    handle = open("data/rice_v6_setaria64/1.fasta")
    #    fh = handle.read()
    #    # double checked these are starts and ends acorrding to coge and match
    #    # with other fasta
    #    self.assertEqual(fh[71774],'N')
    #    self.assertEqual(fh[78937],'N')
    #    #first start
    #    self.assertEqual(fh[71902:71961],'ATGGCGGCTCCCAAGCCCATCTGGGTGCGCCAGGCCGAGGAGGCCAAGCTCAAGTCCGA')
    #    #first end
    #    self.assertEqual(fh[72868:72935],'CAGGCACTGCAGGCACATGCTGCGCAGATGCAGGCTGATTCCAAGGCTGCAGGAGGGGAGGCTTCAG')
    #    #second one start
    #    self.assertEqual(fh[73467:73509],'GTTCAACAGATAAAACTGATAAAGGTGATGTTCTGAAGAAAA')
    #    #third one end
    #    self.assertEqual(fh[75967:76008],'TTGATGGTGCATTTGTGGGCACTGAAGAGTCTGCTATATGA')

    #def test_merge_same_hits(self):
    #    """confrim that genes are being grouped correctly"""
    #    merge_same_hits(self.missed_bed, self.matches, self.sbed)
    #    merge(self.sbed, Bed("data/rice_v6_setaria64/missed_from_rice.bed"),"data/rice_v6_setaria64/rice_v6.all.bed")
    #    new_bed = Bed('data/rice_v6_setaria64/rice_v6.all.bed')
    #    self.assertEqual(new_bed.accn('rice_v6_3_13648324_13648522')['locs'][0][0],13648324)
    #    self.assertEqual(len(new_bed.accn('rice_v6_3_13648324_13648522')['locs']),1)
    #    self.assertEqual(len(new_bed.accn('rice_v6_1_19308831_19309020')['locs']),3)
	## groups to find new genes
	#self.assertEqual(len(new_bed.accn('Os09g38268')['locs']),6)
	## groups exsting genes
	#self.assertEqual(new_bed.accn('Os03g61610')['start'],34920211)
	## confrim not renaming var and giving wrong start
    ## add intersection test...
    def test_group_cds(self):
        qallbed = Bed('data/rice_v6_setaria64/rice_v6.all.bed')
        qaccn = self.qallbed.accn('Os01g01890')
        d,psudo = group_cds('blast_res', qaccn)
        print d.keys()
        self.assertEqual(len(d.keys()),3)
        total_values = []
        for key in d.keys():
            values = len(d[key])
            total_values.append(values)
            print d[key]
        self.assertEqual(sum(total_values),37)
    def test_best_hit(self):
        qaccn = self.qallbed.accn('Os01g01890')
        saccn = self.sallbed.accn('Si000834m')
        d,psudo = group_cds('blast_res', qaccn)
        best_hit(d,qaccn,saccn,psudo)

if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
