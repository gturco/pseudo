import numpy as np
from flatfeature import Bed

def mask_non_cds(bed, mask_with, out):
    mask_with = np.array(mask_with, dtype='c')
    if out is not None and isinstance(out, basestring):
        out = open(out, 'wb')
    for seqid in sorted(bed.fasta.keys()):
        fa = np.array(bed.fasta[seqid], dtype='S1')
        s = fa.shape[0]
        for row in bed[bed['seqid'] == seqid]:
            cds = row['locs']
            non_cds = [(stop,cds[n+1][0]) for n,(start,stop) in enumerate(cds[:-1])]
            non_cds.append((row['start']-1,cds[0][0]))
            non_cds.append((cds[-1][1], row['end']+1))
            for start ,end in non_cds:
                #assert start <= end, (start,end,row)
                #write these to a file
                fa[start: end-1] = mask_with
        assert s == fa.shape[0]
        if out is None:
            yield seqid, fa.tostring()
        else:
            print >> out, ">%s\n%s" % (seqid, fa.tostring())


#bed = Bed('/Users/gturco/setaria64.bed','/Users/gturco/setaria64.fasta'); bed.fill_dict()
##for seqid, seq in mask_non_cds(qbed,'N',None):
##    f = "/Users/gturco/{0}.fasta".format(seqid)
##    fh = open(f,'wb')
##    print >> fh, seq
##    fh.close
#
#for seqid in bed.fasta.keys():
#    f = "/Users/gturco/setaria64_split/{0}.fasta".format(seqid)
#    fh = open(f,'wb')
#    fa = np.array(bed.fasta[seqid], dtype='c')
#    seq = fa.tostring()
#    print >> fh, seq
#    fh.close()


