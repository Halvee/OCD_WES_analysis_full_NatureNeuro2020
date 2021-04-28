#!/usr/bin/env python2.7

import sys
import pandas as pd

# get args
ARGS = sys.argv[1:]
try:
    in_xlsx = ARGS[0]
    genes_recurr_q_lt30_txt = ARGS[1]
    genes_hit_txt = ARGS[2]
except:
    print("get_Wang_2018_TS_genes.py " + \
         "<in.xlsx> " + \
         "<genes_recurr_q_lt30.txt> <genes_hit_txt>")
    sys.exit(1)

# read de novo mutation callset
df = pd.read_excel(in_xlsx, header=0, sheet_name="Non-multiplex families")

# get genes hit at least 2x with q < 0.3
df['dn.lofmis3'] = df['dn.lof'] + df['dn.mis3']
df_recurr = df[df['dn.lofmis3'] > 1]
df_recurr = df_recurr[df_recurr['qval'] < 0.3]
recurr_q_lt30_genes = list(df_recurr['gene.id'])
recurr_q_lt30_genes.sort()
recurr_q_lt30_genes_set = set(recurr_q_lt30_genes)

# write this geneset to file
out_fh = open(genes_recurr_q_lt30_txt, "w")
for gene in recurr_q_lt30_genes: out_fh.write(gene + "\n")
out_fh.close()

# take subset of gene-level results where n denovo mis3>1 OR lof>1
df = df[df['dn.lofmis3']>0]
df = df[df['gene.id'].isin(recurr_q_lt30_genes)==False]
hit_genes = df['gene.id']

# write this geneset to file
out_fh = open(genes_hit_txt, "w")
for gene in hit_genes: out_fh.write(gene+"\n")
out_fh.close()


