import os
import sys
import argparse
import gzip
from bx.intervals import IntervalTree

def main():
    args = parse_args()

    """
    make sure right number of input files are provided
    """
    if len(args.cdnm_files_in_out) % 2 != 0:
        print("ERROR : input files need to be in format " + \
              "in_1.cdnm out_1.cdnm in_2.cdnm out_2.cdnm ..")
        sys.exit(1)

    """
    store every interval into an interval forest data structure,
    where each chromosome has its own interval tree
    """

    interval_forest = dict()
    fh = open_file(args.bed)
    print("reading bed file...")
    for line in fh:
        if line[0] == "#": continue
        data = line.rstrip().split()
        percentile=None
        [chrom, start, end] = data[:3]
        start = int(start)
        end = int(end)
        if args.is_ccr_bed == True:
            percentile=data[3]
            is_ccr_rgn = data[6]
            cov_score = float(data[9])
            if is_ccr_rgn == "VARTRUE": continue
            if args.norm_cvg_min != None:
                cov_score = float(data[9]) 
                norm_cvg = cov_score / (end-start)
                if norm_cvg < args.norm_cvg_min: continue
            if args.syn_dens_min != None:
                syn_dens = float(data[7])
                if syn_dens < args.syn_dens_min: continue
            percentile = float(percentile)
        if chrom not in interval_forest: interval_forest[chrom] = IntervalTree()
        interval_forest[chrom].add(start, end, percentile)
    fh.close()
    print("done.")

    """
    read cdnm file, if ANN annot falls in CCR region, remark annotation
    """
    i = 0
    while i < len(args.cdnm_files_in_out):
        cdnm_file_in = args.cdnm_files_in_out[i]
        cdnm_file_out = args.cdnm_files_in_out[i+1]

        print("Reading DNMs from file " + cdnm_file_in)
        print("Storing bed overlaps to " + cdnm_file_out)

        in_fh = open_file(cdnm_file_in)
        out_fh = open(cdnm_file_out, 'w')

        j = 0
        for line in in_fh:
            j += 1
            data = line.rstrip().split()
            try:
                [iid, chrom, pos, varid, ref, alt, gene, annot] = data[:8]
                pos = int(pos)
            except:
                print("ERROR : improper formatting of cdnm file " + \
                      args.cdnm_file + " at line " + str(i))
                sys.exit(1)
            if chrom not in interval_forest: 
                overlaps = []
            else:
                overlaps = interval_forest[chrom].find(pos-1, pos)
            if args.is_ccr_bed == True:
                if len(overlaps) == 0:
                    max_score = "NA"
                else:
                    max_score = -float('inf')
                    for score in overlaps:
                        if score > max_score: max_score = score
            else:
                if len(overlaps) == 0: 
                    max_score = "FALSE"
                else:
                    max_score = "TRUE"
            out = [iid, chrom, pos, varid, ref, alt, gene, annot, max_score]
            data.append(max_score)
            data = [str(x) for x in data]
            data = "\t".join(data)
            out_fh.write(data + "\n")

        in_fh.close()
        out_fh.close()

        i += 2

    return

def open_file(filename):
    if filename.find(".gz") != -1:
        fh = gzip.open(filename, "rb")
    else:
        fh = open(filename,"r")
    return fh

def parse_args():
    opts = argparse.ArgumentParser()
    opts.add_argument("--is-ccr-bed",action="store_true",default=False,
                      help="input BED file is for CCR data.")
    opts.add_argument("--norm-cvg-min",type=float,
                      default=None,
                      help="min normalized cvg for CCR rgn req for inclusion.")
    opts.add_argument("--syn-dens-min",type=float,
                      default=None,
                      help="min density of synonymous variation in rgn.")
    opts.add_argument("bed", help="BED file with input regions")
    opts.add_argument("cdnm_files_in_out",nargs="+",
            help=".cdnm file input to intersect with CCR, " + \
                 "and output file (<in_1.cdm> <out_1.cdnm> .. " + \
                 "<in_N.cdnm> <out_N.cdnm>")
    args = opts.parse_args()
    return args

if __name__ == "__main__":
    main()
