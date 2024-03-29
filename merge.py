from flatfeature import Bed
from intersection import Intersecter, Feature
import sys

def get_intervening_genes(start,end,seqid, main,saccn):
    inbetween = main.get_features_in_region(seqid, start, end)
    accns = [i['accn']for i in inbetween ]
    intervening_genes = True
    if len(accns)==0 or saccn in accns and len(accns)==1:
        intervening_genes = False
    if abs(end - start) > 7500:
        intervening_genes = True
    return intervening_genes

def merge_same_hits(missed, fh_match, main):
    """ groups genes that hit more then once """
    d = {}
    handle = open(fh_match)
    matches = handle.read()
    main_path = main.path
    path = main_path.split('/')
    dirc = '/'.join(path[:-1])
    org = path[-1]
    fh = open('{0}/missed_from_{1}'.format(dirc,org), "wb")
    for match in matches.split('\n')[:-1]:
        qaccn,saccn = match.split('\t')
        #create dictionary
        try:
            seqid = missed.accn(qaccn)['seqid']
        except KeyError: continue
        if (seqid,saccn) not in d.keys():
            #append whole dict to keys
            d[(seqid,saccn)]= missed.accn(qaccn)
        else:
            #else add locs to exsting one
            gene_start = min(d[(seqid,saccn)]['locs'])[0]
            gene_end = max(d[(seqid,saccn)]['locs'])[1]
            missed_end = missed.accn(qaccn)['locs'][0][1]
            missed_start = missed.accn(qaccn)['locs'][0][0]
            if missed_end < gene_start:
                # if no intervening genes and they are close together...
                intervening_genes = get_intervening_genes(missed_end,gene_start,seqid, main, d[(seqid,saccn)]['accn'])
                if intervening_genes is False:
                    d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
                    d[(seqid,saccn)]['start'] = missed_start
                    if 'Os' in qaccn:
		    	        d[seqid,saccn]['accn'] = qaccn
                else:
                    d[(seqid,qaccn)] = missed.accn(qaccn)
            elif gene_end < missed_start:
                intervening_genes = get_intervening_genes(gene_end,missed_start,seqid, main,d[(seqid,saccn)]["accn"])
                if intervening_genes is False:
                    d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
                    d[(seqid,saccn)]['end'] = missed_end
                    if 'Os' in qaccn:
                        d[seqid,saccn]['accn'] = qaccn
                else:
                    d[(seqid,qaccn)]= missed.accn(qaccn)
            else:
                d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
        
    for key in d.keys():
        new_row = d[key]['locs'].sort()
        row = d[key]
        print >>fh, Bed.row_string(row)



def merge(main, missed, merge_file):
    """creates blast.all file and updates everything"""
    merge_fh = open(merge_file, "w")
    #cds_missed = missed[missed['ftype'] == 'CDS']
    #count = main.shape[0] + missed[missed['ftype'] !='CDS'].shape[0]
    new_rows = []
    seen_accns = {}
    # CDS added to existing gene.
    for row_missed in missed:
        if row_missed['accn'] in seen_accns: continue
        try:
            main_row = main.accn(row_missed['accn'])
             # it's a CDS
        except KeyError:
            #its a new gene
            new_rows.append(row_missed)
            seen_accns[row_missed['accn']] = True
            continue
        locs_interval = Intersecter()
        [locs_interval.add_interval(Feature(start,stop)) for start,stop in main_row['locs']]
        for missed_start,missed_end in row_missed['locs']:
            if len(locs_interval.find(missed_start,missed_end)) > 0:
#                print >>sys.stderr, main_row['accn']
                locs_intersects = {(l.start,l.stop) for l in locs_interval.find(missed_start,missed_end)}
                [main_row['locs'].remove(locs_intersect) for locs_intersect in locs_intersects]
                locs_intersects.add((missed_start,missed_end))
                locs_start = min([start for start,end in locs_intersects])
                locs_end = max([end for start,end in locs_intersects])
                main_row['locs'] = main_row['locs'] + [(locs_start,locs_end)]
                row_missed['locs'].remove((missed_start,missed_end))

        main_row['locs'] = main_row['locs'] + row_missed['locs']
        #print >>sys.stderr, "{0},{1}".format(row_missed['accn'], locs)
        main_row['locs'].sort()
        main_row['start'] = min(min([start for start,end in main_row['locs']]), main_row['start'])
        main_row['end'] = max(max([end for start,end in main_row['locs']]), main_row['end'])
        new_rows.append(main_row)
        seen_accns[main_row['accn']] =True

    for main_rw in main:
        if main_rw['accn'] not in seen_accns:
            new_rows.append(main_rw)
            seen_accns[main_rw['accn']] =True

    def row_cmp(a,b):
        return cmp(a['seqid'], b['seqid']) or cmp(a['start'], b['start'])


    new_rows.sort(cmp=row_cmp)
    #print >>merge_fh, "\t".join(Bed.names)
    for i, row in enumerate(new_rows):
        print >>merge_fh, Bed.row_string(row)

#merge(Bed('data/peach_coco/peach.bed'),Bed('data/peach_coco/missed_from_peach.bed'),'data/peach_coco/peach.all.bed')
#merge_same_hits(Bed('data/peach_coco/missed_peach_from_coco.bed'),'data/peach_coco/missed_peach_from_coco.matches.txt',Bed('data/peach_coco/peach.bed'))

