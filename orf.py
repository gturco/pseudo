import os
import os.path as op
import commands
import operator
from flatfeature import Bed
from pyfasta import Fasta
import numpy as np
import sys
sys.path.append("/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/")
from find_cns import get_pair, remove_crossing_cnss
import re
from Bio.Seq import Seq
pool = None

def remove_fake_orfs(orfs,orf_string,hit_start):
    """weeds out orfs with stop codons and then picks closest """
    real_orfs = []
    for orf_start in orfs:
        first_exon_string = orf_string[orf_start:]
        real_stop_codons = [m.start() for m in re.finditer('TGA',first_exon_string) if (m.start() + orf_start)%3 == 0]
        if len(real_stop_codons) > 0: continue
        #print >> sys.stderr, orf_start
        # which orf_start is closest to the hit start
        real_orfs.append((abs(orf_start - hit_start),orf_start))
    if len(real_orfs) > 0:
        return min(real_orfs)[1]
    else:
        return 'stop codon in orf'

def find_orf(bed, qaccn):
    """ takes query 'gene' and searches 2000bp up first locs for start codon \
        choses based on closest start codon and least stop codons"""
        ### once a reading frame is found there can not be a stop codon within
        ##the frame
    f = Fasta(bed.fasta.fasta_name)
    chromsome = f[qaccn['seqid']]
    if qaccn['strand'] == '+':
        start_pos, search_end = min(qaccn['locs'])
        search_start = start_pos - 2000
        assert search_start < search_end
        hit_start = start_pos
        orf_string = chromsome[search_start:search_end]
        orfs = [m.start() + search_start for m in re.finditer('ATG', orf_string)]
    elif qaccn['strand'] =='-':
        search_start,end_pos = max(qaccn['locs'])
        search_end = end_pos + 2000
        hit_start = end_pos
        my_seq = chromsome[search_start:search_end]
        orf_string = str(Seq(my_seq).reverse_complement())
        orfs = [abs(m.start() - search_end) for m in re.finditer('ATG', orf_string)]
    if len(orfs) > 1:
        best_orf = remove_fake_orfs(orfs,orf_string,hit_start)
        return best_orf
    else: return None

