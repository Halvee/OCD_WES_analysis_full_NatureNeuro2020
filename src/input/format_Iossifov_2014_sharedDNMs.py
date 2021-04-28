#!/usr/bin/env python2.7

import sys
import pandas as pd

# get args
ARGS = sys.argv[1:]
try:
    in_xlsx = ARGS[0]
    out_tsv = ARGS[1]
except:
    print("format_Iossifov_2014_sharedDNMs.py <in.xlsx> <out.tsv>")
    sys.exit(1)

# read de novo mutation callset
df = pd.read_excel(in_xlsx, header=0, sheet_name="Table S2A")

# healthy siblings only
df = df[df.inChild.isin(["sFpF","sFpM", "pFsF", "pMsF",
                         "sMpF","sMpM", "pFsM", "pMsM"])]

# form output df of familyID, chrom, pos, vcfVariant, ref, alt
out_df = df[["familyId", "vcfVariant"]]
out_df = out_df.reset_index()
out_df["chrom"] = None
out_df["pos"] = None
out_df["ref"] = None
out_df["alt"] = None
for i in range(out_df.shape[0]):
    var = str(out_df.loc[i, "vcfVariant"])
    var_info = var.split(":")
    out_df.loc[i, "chrom"] = var_info[0]
    out_df.loc[i, "pos"] = var_info[1]
    out_df.loc[i, "ref"] = var_info[2]
    out_df.loc[i, "alt"] = var_info[3]
    out_df.loc[i, "vcfVariant"] = var.replace(":","-")
out_df=out_df[["familyId", "chrom","pos","vcfVariant","ref","alt"]]
out_df.to_csv(path_or_buf=out_tsv,
              sep="\t", header=False, index=False)
