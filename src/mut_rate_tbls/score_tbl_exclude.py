#!/usr/bin/env python2.7

import os
import sys
import gzip
from collections import defaultdict

CHRS=set(['1','2','3','4','5','6','7','8','9',
          '10','11','12','13','14','15','16','17','18','19',
          '20','21','22','X'])

def main():

    global CHRS

    try:
        ARGS = sys.argv[1:]
        score_tbl_tsv = ARGS[0]
        exclude_varlist_txt = ARGS[1]

    except:
        print(
              "score_tbl_exclude.py " + \
              "<scores.synmislof.tsv> " + \
              "<exclude.varlist>")
        sys.exit(1)

    # init scores dict
    exclude_varlist = set([])

    # read exclude varlist to dict
    varlist_fh = open(exclude_varlist_txt,"r")
    for line in varlist_fh:
        exclude_varlist.add(line.rstrip())
    varlist_fh.close()


    # go through score tbl, store max(score) into data struct for missense
    # go through score tbl, recode as '4' if var is in exclude list
    score_tbl_fh = open(score_tbl_tsv, "r")
    header_str = score_tbl_fh.readline()
    header_str = header_str.rstrip()
    print(header_str)
    header = header_str.split()
    nts = header[4:8]
    for line in score_tbl_fh:
        data = line.rstrip().split()
        [gene, chrom, pos, ref] = data[:4]
        chrom_pos_ref = "-".join([chrom, pos, ref])
        for i in range(len(nts)):
            nt = nts[i]
            if nt == ref: continue
            varid = chrom_pos_ref + "-" + nt
            score_new = float(data[i + 4])
            if varid in exclude_varlist:
                score_new = 4
            data[i + 4] = str(score_new)

        line_out = "\t".join(data)
        print(line_out)
    
    score_tbl_fh.close()
    
    return

def format_intervals(intervals_str):
    intervals_str=intervals_str.replace("(","")
    intervals_str=intervals_str.replace(")","")
    intervals = intervals_str.split(",")
    for i in range(len(intervals)):
        intervals_j = intervals[i].split("..")
        intervals_j = [int(x) for x in intervals_j]
        intervals[i] = intervals_j
    return intervals

if __name__ == "__main__":
    main()
