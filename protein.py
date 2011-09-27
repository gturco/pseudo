import string
from Bio.Seq import Seq
import commands
import sys
import re
from shapely.geometry import LineString
sys.path.append("/Users/gturco/code/freeling_lab/find_cns_gturco/pipeline/scripts/")
from find_cns import remove_crossing_cnss
from flatfeature import Bed
_complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
complement  = lambda s: s.translate(_complement)

def remove_intersecting_hits(hit_list):
    """takes hits and remove if intersect keeping the best bit score"""
    for hit in hit_list:
        hit_list_missing = list(hit_list)
        hit_list_missing.remove(hit)
        qref_hit = LineString([(0.0,hit[0]),(0.0,hit[1]-30)])
        sref_hit = LineString([(0,hit[2]),(0,hit[3]-30)])
        for addt_hit in hit_list_missing:
            qaddt_hit = LineString([(0,addt_hit[0]),(0,addt_hit[1])])
            saddt_hit = LineString([(0,addt_hit[2]),(0,addt_hit[3])])
            if qaddt_hit.intersects(qaddt_hit) or sref_hit.intersects(saddt_hit):
                try:
                    if qref_hit >= qaddt_hit:
                        remove_hit = qaddt_hit
                    else:
                        remove_hit = qref_hit
                        hit_list.remove(remove_hit)
                        hit_list_missing.remove(remove_hit)
                except ValueError: continue
    return hit_list

def protein_fasta(bed,gene,translate,f):
    fh = open(f,"wb")
    #name = ">{0}".format(gene['accn'])
    if gene['strand'] == "+":
        seq = bed.row_cds_sequence(gene['accn'])
        if translate is True:
            codding_dna = Seq(seq)
            protein =  codding_dna.translate()
            seq = protein
        print >>fh, seq
    elif gene['strand'] == "-":
        seq = bed.row_cds_sequence(gene['accn'])
        if translate is True:
            rv_seq = complement(bed.row_cds_sequence(gene['accn']))[::-1]
            codding_dna = Seq(rv_seq)
            protein = codding_dna.translate()
            seq = protein
        print >> fh, seq

def protein_parse(hit,gene,gene_bed, hit_bed):
    "creates a protein fasta and non translated exon fasta \
            blastx them and parse the results"
    hit_fasta = "{0}q.fasta".format('/Users/gturco/code/freeling_lab/pseudo/data/rice_v6_setaria64/')
    gene_fasta = "{0}s.fasta".format('/Users/gturco/code/freeling_lab/pseudo/data/rice_v6_setaria64/')
    if len(re.findall('X',gene_bed.row_cds_sequence(gene['accn']))) >0:
        return "masked", "masked","masked"
    else:
        protein_fasta(hit_bed,hit,False,hit_fasta)
        protein_fasta(gene_bed,gene,True,gene_fasta)
        #cmd = "/Users/gturco/blast-2.2.25/bin/bl2seq -p blastx -G 11 -E 1 -W 3 -e 0.001 -D 1 -i {0} -j {1}  | grep -v '#' | grep -v 'WARNING' | grep -v 'ERROR'".format(hit_fasta,gene_fasta)
        cmd = "/Users/gturco/ncbi-blast-2.2.25+/bin/blastx -gapopen 11 -gapextend 1  -word_size 3 -evalue 0.001 -outfmt '7 gaps qframe std' -query {0} -subject {1} | grep -v '#' | grep -v 'WARNING' | grep -v 'ERROR'".format(hit_fasta,gene_fasta)
        #print >>sys.stderr, "{1} {2} cmd : {0} ".format(cmd,gene,hit)
        res = commands.getoutput(cmd)
        print >>sys.stderr, res
        frame_dict = {'1':{"alignment":[],"qstart":[],"gaps":[]},
                '2':{"alignment":[],"qstart":[],"gaps":[]},
                '3':{"alignment":[],"qstart":[],"gaps":[]},'-1':{"alignment":[],"qstart":[],"gaps":[]},
                '-2':{"alignment":[],"qstart":[],"gaps":[]},
               '-3':{"alignment":[],"qstart":[],"gaps":[]}}
        qhit =[hit['start'], hit['end']]
        sgene =[gene['start'], gene['end']]
        locs_list = []
        for line in res.split("\n"):
            if not line: continue
            if "WARNING:" in line: continue
            if "ERROR" in line: continue
            line = line.split("\t")
            locs = map(int, line[8:12])
            locs.extend(map(float, line[12:]))
            frame = line[1]
            gaps = line[0]
            qstart = min(locs[0],locs[1])
            length = line[5]
            frame_dict[frame]["alignment"].append(length)
            frame_dict[frame]["qstart"].append(qstart)
            frame_dict[frame]["gaps"].append(gaps)
       #frame_lengths = [(sum(frame_dict[key]),key) for key in frame_dict.keys()]
       #frame_lengths.sort()
       #largest_frame = frame_lengths[-1][1]
       #if largest_frame < len....:
       #    frame_shift
      #find stop codon from largest frame + start site....

            locs = tuple(locs)
            locs_list.update((locs,))
        #print >>sys.stderr, "locs_list: {0}".format(locs_list)
        #non_crossing = [(c[0], c[1], c[2], c[3]) for c in remove_intersecting_hits(list(locs_list))]
        non_crossing = [(c[0], c[1], c[2], c[3], c[4]) for c in remove_crossing_cnss(list(locs_list),qhit,sgene)]
        frame_shift = False
        if len(non_crossing) > 1:
            frame_shift = False
        total_hit_len = sum([abs(q_start-q_end) for q_start,q_end,s_start,s_end, evalu in non_crossing])
        total_gene_len = sum(abs(s_start-s_end) for q_start,q_end, s_start,s_end, evalu in non_crossing)
        print >>sys.stderr,non_crossing
        ref_hit_len = len(hit_bed.row_cds_sequence(hit['accn']))
        ref_gene_len = len(gene_bed.row_cds_sequence(gene['accn']))
        print >>sys.stderr,"hit_total {0} \n gene_len {1}".format(total_hit_len,ref_hit_len)
        #print >>sys.stderr, total_hit_len
        hit_len = total_hit_len/float(ref_hit_len)
        gene_len = total_gene_len/float(ref_gene_len/3)
        return hit_len, gene_len, frame_shift
