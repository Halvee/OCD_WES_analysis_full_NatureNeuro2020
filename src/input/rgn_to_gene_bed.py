#!/usr/bin/env python2.7

import sys

def main():
    try:
        ARGS = sys.argv[1:]
        rgn_file = ARGS[0]
    except:
        print("rgn_to_gene_bed.py <in.rgn>")
        sys.exit(1)

    fh = open(rgn_file,"r")
    i = 0
    for line in fh:
        i += 1
        data = line.rstrip().split()
        try:
            [gene,chrom,intervals_str,nbase] = data[:4]
        except:
            print("ERROR : incorrect formatting of rgn file at line " + str(i))
            sys.exit(1)

        intervals_str = intervals_str.replace("(","")
        intervals_str = intervals_str.replace(")","")
        intervals = intervals_str.split(",")
        min_pos=float("inf")
        max_pos=-(float("inf"))
        for interval_str_i in intervals:
            interval_i = interval_str_i.split("..")
            interval_start_i=int(interval_i[0])
            interval_end_i=int(interval_i[1])
            min_pos_i = min([interval_start_i, interval_end_i])
            max_pos_i = max([interval_start_i, interval_end_i])
            if min_pos_i < min_pos: min_pos = min_pos_i
            if max_pos_i > max_pos: max_pos = max_pos_i
        
        out=[chrom, str(min_pos-1), str(max_pos), gene]
        out_str = "\t".join(out)
        print(out_str)

    fh.close()
    
    return


if __name__ == "__main__":
    main()
