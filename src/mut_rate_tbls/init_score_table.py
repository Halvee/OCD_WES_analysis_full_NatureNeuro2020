#!/usr/bin/env python2.7

import sys
from collections import defaultdict
import re

CHRS=set(['1','2','3','4','5','6','7','8','9',
          '10','11','12','13','14','15','16','17','18','19',
          '20','21','22','X'])

def main():

    global CHRS

    try:
        ARGS = sys.argv[1:]
        geneset_txt = ARGS[0]
        ensts_ccds_txt = ARGS[1]
        snpeff_vcf = ARGS[2]
        mu_tbl_file = ARGS[3]
        assert len(ARGS) == 4
    except:
        print("init_score_table.py <geneset.txt> <ccds.enst.txt> <snpeff.vcf> <mu.tbl.tsv>")
        sys.exit(1)

    # read input geneset file, store symbols to set
    scoredict = dict()
    geneset_fh = open(geneset_txt, "r")
    for line in geneset_fh: 
        gene= line.rstrip()
        if gene not in scoredict:
            scoredict[gene] = {}
    geneset_fh.close()

    # read input enst file, store symbols to set
    ensts_ccds = set([])
    enst_fh = open(ensts_ccds_txt, "r")
    for line in enst_fh:
        enst = line.rstrip()
        ensts_ccds.add(enst)
    enst_fh.close()

    # read snpeff_vcf, store values to data structs
    if snpeff_vcf == "stdin":
        snpeff_vcf_fh = sys.stdin
    else:
        snpeff_vcf_fh = open(snpeff_vcf, "r")
    prev_chrom = None
    for line in snpeff_vcf_fh:
        if line[0] == "#": continue
        data = line.rstrip().split()
        [chrom, pos, varid, ref, alt] = data[:5]
        ann = data[7]
        if chrom != prev_chrom:
            # print("loading annots from chromosome " + chrom + ".")
            prev_chrom = chrom
        # trim 'ANN=' from ann string
        anns = ann.split(";")
        for ann in anns:
            ann_keyval=ann.split("=")
            if ann_keyval[0] != "ANN": 
                continue
            tx_annots = ann_keyval[1].split(",")
            score_max = dict()
            for tx_annot_str in tx_annots:
                tx_annot = tx_annot_str.split("|")
                eff = tx_annot[1]
                gene_symbol = tx_annot[3]
                enst = re.sub("\.[0-9]*$","",tx_annot[6])
                if enst not in ensts_ccds: continue
                if gene_symbol not in score_max: score_max[gene_symbol] = -1
                if eff == "stop_gained":
                    score = 2
                elif eff.find("splice_acceptor_variant") != -1:
                    score = 2
                elif eff.find("splice_donor_variant") != -1:
                    score = 2
                elif eff.find("missense") != -1:
                    score = 1
                elif eff.find("synonymous") != -1:
                    score = 0
                else:
                    score = -1
                if score > score_max[gene_symbol]: 
                    score_max[gene_symbol] = score
            for gene_symbol in score_max:
                varid = "-".join([chrom, pos, ref, alt]) 
                if gene_symbol not in scoredict:
                    continue
                if score_max[gene_symbol] == 2:
                    scoredict[gene_symbol][varid] = '2'
                elif score_max[gene_symbol] == 0:
                    scoredict[gene_symbol][varid] = '3'
                elif score_max[gene_symbol] == -1:
                    scoredict[gene_symbol][varid] = '3'
                elif score == 1:
                    scoredict[gene_symbol][varid] = '-0.001'
                else:
                    continue
                #print(score_max[gene_symbol], scoredict[gene_symbol][varid])
    snpeff_vcf_fh.close()

    # open filehandle to mu table
    mu_tbl_fh = open(mu_tbl_file, "r")

    # get header
    mu_header_str = mu_tbl_fh.readline()
    mu_header_str = mu_header_str.rstrip()
    print(mu_header_str)
    mu_header = mu_header_str.split()
    # Gene    Chr     Pos     Ref     A       C       G       T
    nts = mu_header[4:8]
    for line in mu_tbl_fh:
        data = line.rstrip().split()
        [gene, chrom, pos, ref] = data[:4]
        if gene not in scoredict: continue
        scores = ["0", "0", "0", "0"]
        for i in range(len(nts)):
            if nts[i] == ref: continue
            varid = "-".join([chrom, pos, ref, nts[i]])
            if varid in scoredict[gene]:
                scores[i] = scoredict[gene][varid]
        out = [gene,chrom,pos,ref] + scores
        out_str = "\t".join(out)
        print(out_str)
    mu_tbl_fh.close()

    return

if __name__ == "__main__":
    main()
