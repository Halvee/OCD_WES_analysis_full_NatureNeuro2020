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
        annot_dbnsfp_vcf = ARGS[0]
        score_tbl_synlof_tsv = ARGS[1]

    except:
        print(
              "score_tbl_missense_dbnsfp.py " + \
              "<annot.dbnsfp.vcf> " + \
              "<score.synlof.tsv>")
        sys.exit(1)

    # init scores dict
    scores_dict = dict()

    # read annot_dbnsfp_vcf, 
    # store genename -> chrom-pos-ref-alt -> max(score)
    if annot_dbnsfp_vcf == "stdin":
        annot_fh = sys.stdin
    elif annot_dbnsfp_vcf.find(".gz") != -1:
        annot_fh = gzip.open(annot_dbnsfp_vcf, "rb")
    else:
        annot_fh = open(annot_dbnsfp_vcf, "r")
    x=0
    for line in annot_fh:

        if line[0]=="#": continue
        #x+=1
        #if x==100000: break
        
        data = line.rstrip().split()
        [chrom, pos, varid, ref, alt, qual, filter, info] = data[:8]
        
        # get INFO keysvals
        info_dict={}
        info_keysvals = info.split(";")
        for info_keyval_str in info_keysvals:
            info_keyval = info_keyval_str.split("=")
            info_dict[ info_keyval[0] ] = info_keyval[1]

        # skip row if dbNSFP_Polyphen2_HDIV_score not defined
        if "dbNSFP_Polyphen2_HDIV_score" not in info_dict: continue
        
        # get dbNSFP_Polyphen2_HDIV_score value
        ppn2_scores = info_dict["dbNSFP_Polyphen2_HDIV_score"]

        # derive max score from ppn2 scoreset
        score_max = -0.001
        for score_i in ppn2_scores.split(','):
            if score_i == ".": continue
            score_i = float(score_i)
            if score_i > score_max:
                score_max = score_i
        
        # if score_max is defined, store to scores_dict
        if score_max > -0.001:
            if chrom not in scores_dict: scores_dict[chrom] = dict()
            scores_dict[chrom][varid] = score_max

    # go through score tbl, store max(score) into data struct for missense
    score_tbl_fh = open(score_tbl_synlof_tsv, "r")
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
            score = float(data[i + 4])
            if score == 2 or score == 3: continue
            score_new = -0.001
            if varid in scores_dict[chrom]:
                score_new = scores_dict[chrom][varid]
            data[i + 4] = str(score_new)

        line_out = "\t".join(data)
        print(line_out)
    
    score_tbl_fh.close()



    # reopen score_tbl file, drop missense scores into appropriate SNVs,
    # print to stdout

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
