#!/usr/bin/env python2.7

import os
import sys
import gzip
import struct
try:
    from twobitreader import TwoBitFile
except:
    print("ERROR : module 'twobitreader' needed.")

CHRS=set(['1','2','3','4','5','6','7','8','9',
          '10','11','12','13','14','15','16','17','18','19',
          '20','21','22','X'])

def main():

    global CHRS

    try:
        ARGS = sys.argv[1:]
        mut_rate_matrix_file = ARGS[0]
        mer_nbase_flank = int(ARGS[1])
        in_2bit = ARGS[2]
        in_rgn = ARGS[3]
    except:
        print("build_mu_tbl.py <in.mut.rate.matrix>" + \
              "<mer.nbase.flank> <in.2bit> <in.rgn>")
        sys.exit(1)

    # load mutation rate mer data to dictionary
    mut_rates = dict()
    mut_rate_fh = open(mut_rate_matrix_file, "r")
    header = mut_rate_fh.readline()
    header = header.rstrip().split()
    nts = header[1:]
    for line in mut_rate_fh:
        data = line.rstrip().split()
        mer = data[0]
        rates = data[1:]
        mut_rates[mer]=dict()
        for i in range(len(rates)):
            nt = nts[i]
            mut_rate = rates[i]
            mut_rates[mer][nt] = mut_rate
    mut_rate_fh.close()

    # establish header for output table
    #Gene   Chr Pos Ref A   T   C   G
    out_header=["Gene","Chr","Pos","Ref"]
    out_header.extend(nts)
    out_header_str = "\t".join(out_header)

    # init TwoBitFile instance
    twobit = TwoBitFile(in_2bit)

    # print output table header
    print(out_header_str)

    # init rgn filehandle
    rgn_fh = open(in_rgn, "r")
    i = 0
    for line_i in rgn_fh:
        i += 1
        data = line_i.rstrip().split()
        try:
            [gene,chrom,intervals_str,nbase]=data[:4]
        except:
            print("ERROR : incorrect rgn formatting starting at line "+str(i))
            sys.exit(1)

        # skip if chrom isn't chrom 1-22 or X
        if chrom not in CHRS: continue

        intervals = format_intervals(intervals_str)
        for interval_j in intervals:
            seq = twobit[chrom].get_slice(interval_j[0]-1-mer_nbase_flank,
                                          interval_j[1]+mer_nbase_flank)
            zerosite=interval_j[0]
            for site in range(mer_nbase_flank, len(seq)-mer_nbase_flank):
                siteg = interval_j[0] + site - mer_nbase_flank
                ref = seq[site]
                mer = seq[(site-mer_nbase_flank):(site+mer_nbase_flank+1)]
                out = [gene,chrom,str(siteg),ref]
                for nt in nts:
                    out.append(mut_rates[mer][nt])
                out_str = "\t".join(out)
                print(out_str)
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
