from flatfeature import Bed
import sys

def merge_same_hits(missed, fh_match, main):
    """ groups genes that hit more then once """
    d = {}
    handle = open(fh_match)
    matches = handle.read()
    fh = open('/Users/gturco/missed_from_rice.bed', "wb")
    for match in matches.split('\n')[:-1]:
        qaccn,saccn = match.split('\t')
        #create dictionary
        try:
            seqid = missed.accn(qaccn)['seqid']
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
                    inbetween = main.get_features_in_region(seqid, missed_end, gene_start)
                    accns = [i['accn']for i in inbetween ]
                    # if accn in ones inbtween or there are non inbetween
                    # update otherwise give new own name
                    if len(accns)==0 or d[(seqid,saccn)]['accn'] in accns and len(accns)==1:
                        d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
                    else:
                        d[(seqid,qaccn)]= missed.accn(qaccn)

                elif gene_end < missed_start:
                    inbetween = main.get_features_in_region(seqid, gene_end, missed_start)
                    accns = [i['accn']for i in inbetween]
                    if len(accns)==0 or d[(seqid,saccn)]['accn'] in accns and len(accns)==1:
                        d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
                    else:
                        d[(seqid,qaccn)]= missed.accn(qaccn)

                else:
                    d[(seqid,saccn)]['locs'] =  d[(seqid,saccn)]['locs'] + missed.accn(qaccn)['locs']
            if 'Os' in qaccn:
                d[seqid,saccn]['accn'] = qaccn
        except KeyError:continue
    for key in d.keys():
        new_row = d[key]['locs'].sort()
        d[key]['start'] = min(d[key]['locs'])[0]
        d[key]['end'] = max(d[key]['locs'])[1]
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
        main_row['locs'] = main_row['locs'] + row_missed['locs']
        #print >>sys.stderr, "{0},{1}".format(row_missed['accn'], locs)
        main_row['locs'].sort()
        main_row['start'] = min(min(main_row['locs'])[0], main_row['start'])
        main_row['end'] = max(max(main_row['locs'])[1], main_row['end'])
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

merge(Bed('/Users/gturco/rice_v6.bed'),Bed('/Users/gturco/missed_from_rice.bed'),'/Users/gturco/rice_v6.all.bed')
#merge_same_hits(Bed('/Users/gturco/missed_rice_v6_from_setaria64.bed'),'/Users/gturco/missed_rice_v6_from_setaria64.matches.txt',Bed('/Users/gturco/rice_v6.bed'))


