#!/usr/bin/env python2.7

import sys
import pandas as pd

# read input
try:
    ARGS =sys.argv[1:]
    in_xlsx = ARGS[0]
    out_ca_tsv = ARGS[1]
except:
    print("format_Cappi_2017_DNMs.py <in.xlsx> <out.ca.tsv>")
    sys.exit(1)

# read table of DNMs from Cappi et al. 2017
df = pd.read_excel(in_xlsx,
                   header=0, sheet_name="ALL_DE_NOVO_VARIANTS")
# reindex
df = df.reset_index()
# add variant ID column
df.Variant_ID = "NA"
for i in range(df.shape[0]):
    varid_i = "-".join([str(df.loc[i,"chr"]), str(df.loc[i,"pos_vcf"]),
                        str(df.loc[i,"ref_vcf"]), str(df.loc[i,"alt_vcf"])])
    df.loc[i,"Variant_ID"] = varid_i
# get case subset - ocd samples only, no exclusions
df_ca = df[(df.batch == "ocd")]
df_ca = df_ca[df_ca.exclude == 0]
df_ca = df_ca.reset_index()
df_ca = df_ca.loc[:, ["Proband","chr","pos_vcf","Variant_ID","ref_vcf","alt_vcf"]]
# get control subset - ssc unaffected siblings, no exclusions
df_co = df[(df.batch == "sscsib")]
df_co = df_co[df_co.exclude == 0]
df_co = df_co.reset_index()
df_co = df_co.loc[:, ["Proband","chr","pos_vcf","Variant_ID","ref_vcf","alt_vcf"]]
# write output files
df_ca.to_csv(path_or_buf=out_ca_tsv,
             sep="\t", header=False, index=False)
