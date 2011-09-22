import unittest
import collections

import sys
import itertools
import logging

logging.basicConfig(level=logging.INFO)
sys.path.append("../scripts")
from flatfeature import Bed
from pseudo import bites, group_cds,main, remove_crossing_hits, append_to_included_groups
from orf import find_orf

class TestPseudo(unittest.TestCase):
    def setUp(self):
        self.qallbed = Bed('data/rice_v6_setaria64/rice_v6.all.bed','data/rice_v6_setaria64/rice_v6.fasta') ;self.qallbed.fill_dict()
        self.sallbed = Bed('data/rice_v6_setaria64/setaria64.all.bed','data/rice_v6_setaria64/setaria64.fasta') ;self.sallbed.fill_dict()
        self.saccn = self.sallbed.accn('Si000834m')
        blastfh = open('blast_res')
        self.blast = blastfh.read()
        self.d, self.pseudo = group_cds(self.blast,self.saccn)

    def test_group_cds_1(self):
        self.assertEqual(len(self.d.keys()),4)
        total_values = []
        for key in self.d.keys():
            values = len(self.d[key])
            total_values.append(values)
        self.assertEqual(sum(total_values),38)

    def test_group_cds_2(self):
        blast_2fh = open('blast_2')
        blast_2 = blast_2fh.read()

        d,pseudo = group_cds(blast_2,self.sallbed.accn('Si002524m'))

        self.assertEqual(len(d.keys()),5)
        for key in d.keys():
            #logging.info('key: {0}'.format(key))

            self.assertEqual (1, len(d[key]))

    def test_append_to_included_groups(self):
        locs = [1,2,3,4]
        group_dict = {(2,5):[], (3,6):[], (9,8):[]}

        result_dict = append_to_included_groups(locs, group_dict)
        expected = {(2,5):[(1, 2, 3, 4)], (3,6):[(1, 2, 3, 4)], (9,8):[]}

        self.assertEquals(expected, result_dict)

    def test_remove_crossing_hit(self):
        qaccn = self.qallbed.accn('Os01g01890')
        for group_key in self.d.keys():
            exon_hits = self.d[group_key]
            non_crossing =remove_crossing_hits(exon_hits,qaccn,self.saccn)
            if len(non_crossing) > 1:
                mid,start,stop =bites(non_crossing)

    def test_find_orf(self):
        qaccn = self.qallbed.accn('Os01g01295')
        orf = find_orf(self.qallbed, qaccn)
        self.assertEqual(orf+1, 141084)

    def test_find_orf_neg(self):
        saccn = self.sallbed.accn('Si001539m')
        orf = find_orf(self.sallbed,saccn)
        self.assertEqual(orf,7662777)
if __name__ == '__main__':
    unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSequenceFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)
