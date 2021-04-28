#!/usr/bin/env python2.7

import sys
import pandas as pd

def main():
    ARGS = sys.argv[1:]
    try:
        satterstrom_2020_ASD_xlsx = ARGS[0]
        coe_2018_NDD_geneset = ARGS[1]
        loeuf_tsv = ARGS[2]
        out_geneset_csv = ARGS[3]
        out_geneset_txt = ARGS[4]
    except:
        print("01.form_NDD_geneset_file.py" + \
              "<satterstrom_2020_ASD.xlsx> <coe_2018_NDD.txt> " + \
              "<loeuf.tsv> <out.geneset.csv> <out.geneset.txt>")
        sys.exit(1)

    # read sheet from satterstrom 2020 xlsx with 102 ASD risk genes
    asd_df = pd.read_excel(satterstrom_2020_ASD_xlsx, 
                           header=0, sheetname="102_ASD")
    asd_genes=set(asd_df.loc[0:101,"gene"])

    # read geneset from Coe 2018 NDD (n=124 exomewide signif genes)
    ndd_genes = set()
    fh = open(coe_2018_NDD_geneset, "r")
    for line in fh: ndd_genes.add(line.rstrip())
    fh.close()
    
    # read tsv that should have LOEUF bins in it
    loeuf = pd.read_csv(loeuf_tsv, header=0, sep="\t")   

    # form ASD / DD geneset
    asd_ndd_genes = asd_genes.union(ndd_genes)
    asd_ndd_genes = list(asd_ndd_genes)
    asd_ndd_genes.sort()
    asd_ndd = set(asd_ndd_genes)
    
    # print number of genes per geneset to stdout
    print("N genes (ASD, Satterstrom et al. 2020) : " + str(len(asd_genes)))
    print("N genes (NDD, Coe et al. 2019) : " + \
          str(len(ndd_genes)))
    print("N genes (ASD only) : " + str(len(asd_genes.difference(ndd_genes))))
    print("N genes (NDD only) : " + str(len(ndd_genes.difference(asd_genes))))
    print("N genes (ASD or NDD) : " + str(len(asd_ndd_genes)))    
    asd_ndd_intersect=ndd_genes.intersection(asd_genes)
    print("N genes (ASD and NDD) : " + str(len(asd_ndd_intersect)))

    # healthy siblings only
    # df = df[df.inChild.isin(["sF","sM"])]

    # make stats df that lists for a given gene which genesets it's in
    # get asd genes intersect with ASD / DD joint geneset
    asd_genes_intersect=[]
    asd_genes_set=set(asd_genes)
    for gene in asd_ndd_genes:
        overlap = 0
        if gene in asd_genes_set: overlap = 1
        asd_genes_intersect.append(overlap)
    ndd_genes_intersect=[]
    ndd_genes_set=set(ndd_genes)
    for gene in asd_ndd_genes:
        overlap = 0
        if gene in ndd_genes_set: overlap = 1
        ndd_genes_intersect.append(overlap)

    # make stats df that lists for a given gene which genesets it's in
    g_df = pd.DataFrame({"gene":asd_ndd_genes,
                         "ASD_Satterstrom_2020":asd_genes_intersect,
                         "NDD_Coe_2019":ndd_genes_intersect})

    # merge with loeuf
    g_df = g_df.merge(loeuf, on="gene", how="left")
    g_df = g_df[["gene","gene_id","oe_lof_upper_bin",
                 "ASD_Satterstrom_2020","NDD_Coe_2019"]]
    g_df.rename(columns={"oe_lof_upper_bin":"gnomAD_oe_lof_upper_bin"})

    # write geneset csv to file
    g_df.to_csv(path_or_buf=out_geneset_csv, sep=",", index=False)

    # write geneset to file
    asd_ndd_genes = list(asd_ndd_genes)
    asd_ndd_genes.sort()
    out_fh = open(out_geneset_txt,"w")
    for gene in asd_ndd_genes:
        out_fh.write(gene+"\n")
    out_fh.close()    

    return

if __name__ == "__main__":
    main()
