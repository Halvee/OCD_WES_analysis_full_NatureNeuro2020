
import os
import sys
import argparse
import gzip
from bx.intervals import IntervalTree

# global vars
DMG_THRESH = 0.957
NTS = ["A","C","G","T"]
FRAMESHIFT_NONSENSE_RATIO = 1.25 # taken from Samocha et al. 2015
ANNOTS=["length", "syn","misU","misB","misP","misD","misND","mis",
        "non","splice", "lof_snv","frameshift","lof","dmg"] 
def main():
    global DMG_THRESH
    global NTS
    global FRAMESHIFT_NONSENSE_RATIO
    global ANNOTS

    args = parse_args()

    if os.path.isfile(args.n_alleles):
        n_alleles = []
        alleles_tbl_fh = open(args.n_alleles, "r")
        for line in alleles_tbl_fh: 
            n_alleles.append( int(line.rstrip()) )
        alleles_tbl_fh.close()
    else:
        try:
            n_alleles = int( args.n_alleles )
        except:
            print("ERROR : if n_alleles is not a file, must be an int.")
            sys.exit(1)    

    # read parent geneset, init output df
    gs = set()
    gs_fh = open(args.parent_geneset, "r")
    for line in gs_fh: gs.add( line.rstrip() )
    gs_fh.close()

    # init data structures for storing mu rates
    exome_len = 0
    exome_lencds = 0
    gene_len = dict()
    gene_lencds = dict()
    mu_tot = dict()
    for gene in gs:
        gene_len[ gene ] = 0
        gene_lencds[ gene ] = 0
        mu_tot[ gene ] = dict()
        for annot in ANNOTS:
            mu_tot[ gene ][ annot ] = 0
    
    # init filehandles for mu and ppn2 tbls
    mu_fh = open(args.mu_tbl, "r")
    ppn2_fh = open(args.ppn2_tbl, "r")
    mu_header = mu_fh.readline()
    ppn2_header = ppn2_fh.readline()
    prev_gene = None
    h = -1
    while(1):
        h += 1
        mu_line = mu_fh.readline().rstrip()
        ppn2_line = ppn2_fh.readline().rstrip()
        if mu_line == "" and ppn2_line == "":
            break
        mu_data = mu_line.split("\t")
        ppn2_data = ppn2_line.split("\t")
        if len(mu_data) != 8 or len(ppn2_data) != 8:
            continue
        mu_gene = mu_data[0]
        ppn2_gene = ppn2_data[0]
        if mu_gene not in mu_tot: continue
        if mu_gene != prev_gene:
            prev_gene = mu_gene
        chrom = mu_data[1]
        pos = mu_data[2]
        exome_len += 1
        gene_len[ mu_gene ] += 1

        # if only looking at autosomal genes and chrom is not 1-22, skip
        if args.autosomes_only == True and (chrom == "X" or chrom == "Y"):
            continue

        if type(n_alleles) == int:
            n_alleles_h = n_alleles
        else:
            n_alleles_h = n_alleles[ h ]
        n_lof = 0
        splice_mu = 0
        ref_nt = mu_data[3]
        if ref_nt == "A":
            idx = (5,6,7)
        elif ref_nt == "C":
            idx = (4,6,7)
        elif ref_nt == "G":
            idx = (4,5,7)
        elif ref_nt == "T":
            idx = (4,5,6)
        else:
            continue
        for i in idx:
            alt_nt = NTS[i-4]
            varid = "-".join([chrom, pos, ref_nt, alt_nt])
            is_mis = True
            is_misd = False
            if mu_data[i] == "NA" or ppn2_data[i] == "NA":
                continue    
            mu_i = float( mu_data[i] ) * n_alleles_h
            ppn2_i = float( ppn2_data[i] )
            if ppn2_i == 4:
                continue
            if ppn2_i == 3:
                mu_tot[ mu_gene ][ "syn" ] += mu_i
                is_mis = False
            elif ppn2_i == 2:
                mu_tot[ mu_gene ][ "lof_snv" ] += mu_i
                is_mis = False
                n_lof += 1
                splice_mu += mu_i
            elif ppn2_i >= 0.957:
                mu_tot[ mu_gene ][ "misD" ] += mu_i
                mu_tot[ mu_gene ][ "mis" ] += mu_i
                is_misd = True
            elif ppn2_i >= 0.453:
                mu_tot[ mu_gene ][ "misP" ] += mu_i
                mu_tot[ mu_gene ][ "misND" ] += mu_i
                mu_tot[ mu_gene ][ "mis" ] += mu_i 
            elif ppn2_i >= 0:
                mu_tot[ mu_gene ][ "misB" ] += mu_i
                mu_tot[ mu_gene ][ "misND" ] += mu_i
                mu_tot[ mu_gene ][ "mis" ] += mu_i
            else:
                mu_tot[ mu_gene ][ "misU" ] += mu_i
                mu_tot[ mu_gene ][ "mis" ] += mu_i

        # capture info specific to splice donor/acceptor sites
        if n_lof == 3:
            mu_tot[ mu_gene ][ "splice" ] += splice_mu
        else:
            gene_lencds[ mu_gene ] += 1
            exome_lencds += 1

    # for each gene, derive mutation rate specific for nonsense.
    # if ratio defined, estimate frameshift rate based on nonsense
    # rate multiplied by ratio of frameshift DNMs to nonsense DNMs
    # in trios.
    rate_non_tot = 0
    rate_snv_tot = 0
    for gene in mu_tot:
        mu_gene = mu_tot[ gene ]
        mu_gene[ "non" ] = mu_gene[ "lof_snv" ] - mu_gene[ "splice" ]

        # iterate total snv rate
        rate_snv_tot += mu_gene[ "syn" ]
        rate_snv_tot += mu_gene[ "misB" ]
        rate_snv_tot += mu_gene[ "misP" ] 
        rate_snv_tot += mu_gene[ "misD" ]
        rate_snv_tot += mu_gene[ "misU" ]
        rate_snv_tot += mu_gene[ "lof_snv" ]

        # iterate total nonsense mu rate
        rate_non_tot += mu_gene[ "non" ]

    # set exomic exp. frameshift rate as exp_non_rate_total multiplied by
    # FRAMESHIFT_NONSENSE_RATIO. Set exp. rate per gene as total exp.
    # rate divided by number of coding bases in each gene
    rate_fs_tot = rate_non_tot * FRAMESHIFT_NONSENSE_RATIO
    for gene in mu_tot:
        mu_gene = mu_tot[ gene ]
        portion_exome = float(gene_lencds[ gene ]) / float(exome_lencds)
        mu_gene[ "frameshift" ] = rate_fs_tot * portion_exome
        mu_gene["lof"] = mu_gene["lof_snv"] +  mu_gene["frameshift"]
        mu_gene["dmg"] = mu_gene["lof"] + mu_gene["misD"]

        # add total basepairs for gene
        mu_gene["length"] = float(gene_len[ gene ])

    header = ["Gene"] + ANNOTS
    print("\t".join(header))
    for gene in gs:
        if gene not in mu_tot: continue
        out_list = [gene]
        for annot in ANNOTS:
            exp_rate = mu_tot[ gene ][ annot ]
            out_list.append( str( exp_rate ) )
        out_str = "\t".join(out_list)
        print(out_str)
    return

def parse_args():
    opts = argparse.ArgumentParser()
    opts.add_argument("--autosomes-only", action="store_true", default=False,
                      help="compile stats for autosomal genes only (chr1-22)")
    opts.add_argument("parent_geneset", help="file with parent geneset symbols")
    opts.add_argument("mu_tbl", help="file with mutation rates per possible variant.")
    opts.add_argument("ppn2_tbl", help="file with ppn2 / lof /syn scores per possible variant.")
    opts.add_argument("n_alleles", help="numbers of input trios OR " + \
                                        "file where each line has number of sufficiently covered " + \
                                        "alleles in input trio cohort, lines up with mu/ppn2 tbls.")
    args = opts.parse_args()
    return args

def open_file(filename):
    if filename.find(".gz") != -1:
        fh = gzip.open(filename, "rb")
    else:
        fh = open(filename,"r")
    return fh

if __name__ == "__main__":
    main()
