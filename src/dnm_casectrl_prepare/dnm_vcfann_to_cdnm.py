
import sys
import gzip
from collections import defaultdict
import re

func_scores = {"synonymous_variant":0,
               "splice_region_variant&synonymous_variant":0,
               "conservative_inframe_insertion":0,
               "conservative_inframe_deletion":0,
               "disruptive_inframe_insertion":0,
               "disruptive_inframe_deletion":0,
               "missense_variant":1,
               "missense_variant&splice_region_variant":1,
               "frameshift_variant":2,
               "stop_gained":2,
               "stop_lost":0,
               "start_lost":0,
               "splice_donor_variant&intron_variant":2,
               "splice_acceptor_variant&intron_variant":2}
func_trans = {"synonymous_variant":"syn",
              "splice_region_variant&synonymous_variant":"syn",
              "missense_variant":"mis",
              "missense_variant&splice_region_variant":"mis",
              "conservative_inframe_deletion":"del",
              "conservative_inframe_insertion":"ins",
              "disruptive_inframe_insertion":"ins",
              "disruptive_inframe_deletion":"del",
              "frameshift_variant":"frameshift",
              "stop_gained":"non",
              "stop_lost":"stoplost",
              "start_lost":"startlost",
              "splice_donor_variant&intron_variant":"splice",
              "splice_acceptor_variant&intron_variant":"splice"}
maf_names=set(["dbNSFP_ExAC_NFE_AF",
               "dbNSFP_ExAC_SAS_AF",
               "dbNSFP_ExAC_Adj_AF",
               "dbNSFP_ExAC_AFR_AF",
               "dbNSFP_ExAC_AF",
               "dbNSFP_ExAC_FIN_AF",
               "dbNSFP_ExAC_AMR_AF",
               "dbNSFP_ExAC_EAS_AF"])
maf_thresh=0.0005

ARGS = sys.argv[1:]
try:
    enst_file = ARGS[0]
    dnm_file = ARGS[1]
    vcfann_file = ARGS[2]
    cdnm_file = ARGS[3]
except:
    print("dnm_vcfann_to_cdnm.py <enst.file> <dnm_file> <vcfann_file> <cdnm_file>")
    sys.exit(1)

## initialize cdnm output file
out_fh = open(cdnm_file, "w")
out_fh.write("")
out_fh.close()

## read list of ensts to subset on
ensts=set([])
enst_fh = open(enst_file, "r")
for line in enst_fh: ensts.add(line.rstrip())
enst_fh.close()

## read dnm file calls into memory
dnm_fh = open(dnm_file, "r")
dnm_calls = defaultdict(set)
for line in dnm_fh:
    data = line.rstrip().split()
    sample_id = data[0]
    var_id = data[3]
    dnm_calls[ var_id ].add(sample_id)
dnm_fh.close()

## read annotated vcf file
if vcfann_file.find(".gz") != -1:
    vcfann_fh = gzip.open(vcfann_file, "rb")
else:
    vcfann_fh = open(vcfann_file, "r")
for line in vcfann_fh:
    if line[0] == "#": continue
    data = line.rstrip().split()
    [chrom, pos, varname, ref,
     alt, qual, filter, info] = data[:8]
    var_id = "-".join([chrom, pos, ref, alt])
    if var_id not in dnm_calls: continue

    ## get annotations from info
    ppn2_score = None
    maf_popmax = 0
    max_func = None
    max_func_trans = None
    max_func_gene = None
    info_list = info.split(";")
    score_max = -1
    func_max = -1
    gene_max = -1
    score_i = 0
    ppn2_hdiv_dmg = False
    for info_i in info_list:
        info_i_keyval = info_i.split("=")
        if len(info_i_keyval) != 2: continue
        i_key = info_i_keyval[0]
        i_val = info_i_keyval[1]
        if i_key in maf_names:
            i_val=max([float(x) for x in i_val.split(",")])
            if float(i_val) > maf_popmax: maf_popmax = float(i_val)
        elif i_key=="dbNSFP_Polyphen2_HDIV_pred":
            if i_val.find("D") != -1: ppn2_hdiv_dmg = True
        elif i_key=="ANN":
            for eff_tx_full in i_val.split(","):
                eff_tx_set = eff_tx_full.split("|")
                eff_tx = eff_tx_set[1]
                eff_tx_gene = eff_tx_set[3]
                transcript=re.sub("\.[0-9]*$","",eff_tx_set[6])
                if transcript not in ensts: continue
                if eff_tx not in func_scores: continue
                eff_tx_score = func_scores[ eff_tx ]
                if eff_tx_score > score_max:
                    score_max = eff_tx_score
                    func_max = eff_tx
                    gene_max = eff_tx_gene
                    max_func_trans = func_trans[ func_max ]
    
    # skip variant if popmax maf greater than threshold
    if maf_popmax > maf_thresh: continue

    if func_max == "missense_variant":
        if ppn2_hdiv_dmg == True:
            max_func_trans = "misD"
    
    if var_id not in dnm_calls: continue
    if max_func_trans == None: continue
    burdened_samples = list( dnm_calls[ var_id ] )
    burdened_samples.sort()
    for sample_id in burdened_samples:
        out_list = [sample_id, chrom, pos, var_id,
                    ref, alt, gene_max, max_func_trans]
        out_list = [str(i) for i in out_list]
        out_str = "\t".join( out_list )
        out_fh = open(cdnm_file, "a")
        out_fh.write(out_str + "\n")
        out_fh.close()


vcfann_fh.close()
