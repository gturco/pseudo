#Python pseudo.py \
#        -n 8 \
#        --qfasta \
#        --qbed \
#        --sfasta \
#        --sbed \
#        --miss \ > test.txt




"""tblastx: faster,what coge uses, done have exact postion where gene protein
starts.. makes easier for that"""
import os
import os.path as op
import commands
import operator
from flatfeature import Bed
import numpy as np
from processing import Pool
import sys
sys.path.append("/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/")
from find_cns import get_pair, remove_crossing_cnss
from mask_non_cds import mask_non_cds
from pyfasta import Fasta
import re
import logging
from shapely.geometry import LineString
pool = None
logging.basicConfig(level=logging.INFO)

def group_cds(blast_str, ref_gene):
    """groups blast hits to the ref gene
    creates two dict the first is a list of locational hits for each exon of the ref
    genome the other contatins addtional information for theses hits"""
    cds = ref_gene['locs']
    group_list = [(start,end) for (start,end) in cds]
    d = {key:[] for key in group_list}
    psudo = {}
    #blast_fh = open(blast_str)
    #blast_str = blast_fh.read()
    for line in blast_str.split("\n"):
        if not line: continue
        if "WARNING" in line: continue
        if "ERROR" in line: continue
        line = line.split("\t")
        locs  = map(int, line[6:10])
        locs.extend(map(float,line[10:]))
        print >>sys.stderr,'llllllllllllllll{0}'.format(line)
        percent,length =map(float,line[2:4])
        psudo[str(tuple(locs[:-1]))] = (length,percent)
        #logging.info("locs:{0}".format(locs))
        #logging.info("d:{0}".format(d))
        append_to_included_groups(locs, d)
    return d, psudo

def append_to_included_groups(locs, d):
    for group_key in d.keys():
        ref_gene_poly = LineString([(0.0, float(group_key[0])),(0.0, float(group_key[1]))])
        locs_poly = LineString([(0, locs[2]), (0, locs[3])])
        if ref_gene_poly.intersects(locs_poly):
            d[group_key].append(tuple(locs))
    return d

def remove_crossing_hits(exon_hits,qfeat,sfeat):
        """uses find cns remove overlaping hits and corssing on each exon"""
        qgene =[qfeat['start'], qfeat['end']]
        sgene =[sfeat['start'], sfeat['end']]
        orient = qfeat['strand'] == sfeat['strand'] and 1 or -1
        if orient == -1:
            exon_hits = list(exon_hits)
            for i, hit in enumerate(exon_hits):
                hit = list(hit)
                hit[2] *= -1
                hit[3] *= -1
                hit[i] = tuple(hit)
            sgene[0] *= -1
            sgene[1] *= -1
        non_crossing_hits = [(c[0], c[1], c[2], c[3], c[-2]) for c in remove_crossing_cnss(exon_hits,qgene,sgene)]
        if orient == -1:
            non_crossing_hits == [(c[0],c[1],-c[2], -c[3], c[-1]) for c in remove_crossing_cnss(exon_hits,qgene,sgene)]
        #non_crossing_dict = {str(locs):psudo[str(locs)] for locs in non_crossing_hits}
        return non_crossing_hits
        
def bites(hit_dict):
    """for each exon compares blast info to find bites
    returns start and stop of loc"""
    x = [(qstart,qend) for qstart,qend,sstart,send,evalue in hit_dict]
    #print >>sys.stderr,x
    x.sort()
    mid = [abs(hstop-x[n+1][0]) for n,(hstart,hstop) in enumerate(x[:-1])]
    start = min(x)[0]
    stop= max(x)[1]
    return mid,start,stop

#        #best_hit = max(non_crossing_dict.iteritems(), key=operator.itemgetter(1))[0]
#        hit_lengths = [float(length) for (length,percent) in non_crossing_dict.values()]
#        best_hits.append(sum(hit_lengths))
#    total_cds = sum([abs(qstart-qend) for (qstart,qend) in d.keys()])
#    total_hits = sum(best_hits)*3
#    x = (total_hits/float(total_cds))
#    print >>sys.stderr, x
##    if avg > 70:
##        print

def get_mask_non_cds(bed):
    f = bed.fasta.fasta_name
    fname = op.splitext(op.basename(f))[0]
    d = op.dirname(f) + "/{0}_split".format(fname)
    try: os.mkdir(d)
    except OSError: pass

    fastas = {}
    for seqid, seq in mask_non_cds(bed, 'N', None):
        f = d + "/{0}.fasta".format(seqid)
        fastas[seqid] = f
        if op.exists(f): continue
        fh = open(f,"wb")
        print >>fh, seq
        fh.close()
    return fastas


def split_fastas(bed):
   """mkdir fasta_split and puts fastas split by seqid/chr in it"""
   # do we want non features masked?
   f = bed.fasta.fasta_name
   fname = op.splitext(op.basename(f))[0]
   d = op.dirname(f) + "/{0}_split".format(fname)
   try: os.mkdir(d)
   except OSError:pass

   fastas = {}
   for seqid in bed.fasta.keys():
       f = d + "/{0}.fasta".format(seqid)
       fastas[seqid] = f
       if op.exists(f):continue
       fh = open(f,"wb")
       fa = np.array(bed.fasta[seqid], dtype='c')
       seq = fa.tostring()
       print >> fh, seq
       fh.close()
   return fastas
   """pyfasta split --header "%(seqid)s.fasta" original.fasta"""

#sbed = Bed('/Users/gturco/setaria64.bed','/Users/gturco/setaria64.fasta'); sbed.fill_dict()
#split_fastas(sbed)
def get_pairs_file(missed_pair_file):
    """takes the missed_gene.match file
    and creates a file that can be used in get pair by creating a one to one
    pair ratio"""
    handle = open(missed_pair_file)
    fh = handle.read()
    d = missed_pair_file.split("/")
    old_name = d[-1].split('.')[0]
    d= '/'.join(d[:-1])
    f = "{0}/{1}.pairs.txt".format(d, old_name)
    w = open(f,"wb")
    for line in fh.split('\n')[:-1]:
        query, subject = line.split('\t')
        for gene in subject.split(','):
            new_pair = "{0}\t{1}".format(query,gene)
            print >> w, new_pair
    w.close()
    return f

#get_pairs_file("/Users/gturco/missed_exons/missed_rice_v6_from_sorghum_v1.matches.txt")
def main(qbed,sbed,missed_pairs, ncpu):
    """run tblastx on missed pairs..."""
    pool = Pool(ncpu)

    pairs_file = get_pairs_file(missed_pairs)

    #tblastx = "/home/gturco/src/blast-2.2.25/bin/bl2seq -p tblastx -G 11 -E 1 -W 3 -e 1e-10 -D 1 -i {0} -j {1} -I {2},{3} -J {4},{5}"
    blastn = "/Users/gturco/blast-2.2.25/bin/bl2seq -p blastn -G 5 -E 2 -W 7 -q -2 -e 0.001 -D 1 -i {0} -j {1} -I {2},{3} -J {4},{5} | grep -v '#' | grep -v 'WARNING' | grep -v 'ERROR' "
    qfastas = split_fastas(qbed)#remians unmasked
    sfastas = split_fastas(sbed) #mask noncoding

    pairs = [True]
    _get_pair_gen = get_pair(pairs_file,"pair", qbed,sbed)

    def get_pair_gen():
        try: return _get_pair_gen.next()
        except StopIteration: return None
        
    while any(pairs):
        pairs = [get_pair_gen() for i in range(ncpu)]
        
        def get_blastn_cmd(pair):
            """creates the dictionary values used to fill in blast cmd"""
            if pair is None: return None
            hit, gene = pair
            hstart, hstop = abs(3000 - hit['start']), (3000 + hit['end'])
            # double check fasta to make sure i dont need to add or remove one
            gstart,gstop = gene['start'],gene['end']
            # checks the entire gene...
            # do we want to add padding?
            query_file = qfastas[hit['seqid']]
            subject_file = sfastas[gene['seqid']]

            cmd = blastn.format(query_file, subject_file, hstart, hstop, gstart, gstop)
            #print >> sys.stderr,'{0},{1},{2}'.format(hit['accn'],gene['accn'],cmd)
            
            return cmd, hit, gene
        cmds = [c for c in map(get_blastn_cmd, [l for l in pairs if l]) if c]
        results = (r for r in pool.map(commands.getoutput,[c[0]for c in cmds]))
        for res, (cmd, hit, gene) in zip(results,cmds):
            d,psdu = group_cds(res, gene)
            gap_list =[]
            intron_list = []
            for group_key in d.keys():
                exon_hits = d[group_key]
                non_crossing = remove_crossing_hits(exon_hits,hit,gene)
                if len(non_crossing) > 1:
                    gaps,hstart,hend =bites(non_crossing)
                    gap_list.append(map(str,gaps))
                    intron_list.append(group_key[0])
                if len(non_crossing) == 1:
                    intron_list.append(group_key[0])
            introns = "{0}/{1}".format(len(intron_list),len(gene['locs']))
            gap_list = map(str,gap_list)
            w ="{0},{1},{2},{3}".format(hit['accn'],gene['accn'],"/t".join(gap_list),introns)
            print >> sys.stdout, w
#            if float(values[2]) >= 70.0:
#                blastn_hits.append((percent_int,evalue))
#            if len(tblastx_hits) == 0: continue
#            w = "{0}\t{1}\t{2}\t".format(hit['accn'], gene['accn'], ",".join(tblastx_hits))
#            print >> sys.stdout, w
#


if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser("usage: %prog [options] ")
    parser.add_option("-n", dest="ncpu", help="parallelize to this many cores", type='int', default=2)
    parser.add_option("--qbed", dest="qbed", help="query bed file/bed file with  hits/new genes")
    parser.add_option("--sbed", dest="sbed", help="subject bed file/ exsiting genes" )
    parser.add_option("--qfasta", dest="qfasta", help="query fasta file/hits/new genes")
    parser.add_option("--sfasta", dest="sfasta", help="subject fasta file /exsiting genes")
    parser.add_option("--miss", dest="missed_pairs", help="missed pairs file query_missed_from_subject.matches") 
    (options, _) = parser.parse_args()

    sbed = Bed(options.sbed, options.sfasta); sbed.fill_dict()
    qbed = Bed(options.qbed, options.qfasta); qbed.fill_dict()
    main(qbed,sbed,options.missed_pairs,options.ncpu)


