
import os
import sys
from bx.intervals.intersection import IntervalTree

def main():
    ARGS = sys.argv[1:]
    try:
        in_vcf = ARGS[0]
        in_rgn = ARGS[1]
    except:
        print("vcf_rgn_intersect.py <in.dnm> <in.rgn>")
        sys.exit(1)

    # read rgn data
    rgn_fh = open(in_rgn, "r")
    intervaltreeset = rgn_to_intervaltreeset(rgn_fh)

    if in_vcf != "stdin":
        vcf_fh = open(in_vcf, "r")
    else:
        vcf_fh = sys.stdin

    inputcoord_intersect(vcf_fh, intervaltreeset)
    rgn_fh.close()
    vcf_fh.close()
    return


def inputcoord_intersect(inputcoord_fh, intervaltreeset, 
                         input_type = "vcf",
                         inverse = False):

    i = 0
    for line in inputcoord_fh:
        i += 1
        if line[0] == "#": 
            print(line.rstrip())
            continue
        data = line.rstrip().split()
        [chrom,pos,varid,ref,alt] = data[:5]
        try:
            chrom = chrom.replace("chr","")
            start = int(pos)
            end = start + len(ref) - 1
        except:
            print("malformed inputcoord file at line "+str(i))
            sys.exit(1)

        overlaps = []
        if chrom not in intervaltreeset: continue
        overlaps = intervaltreeset[chrom].find(start, end)
        if (len(overlaps) > 0 and inverse == False): 
            print(line.rstrip())
        elif inverse == True and len(overlaps) == 0:
            print(line.rstrip())
    inputcoord_fh.close()
    return

def bed_to_intervaltreeset(bed_fh):
    intervaltreeset = dict()
    i = 0
    for line in bed_fh:
        i += 1
        data = line.rstrip().split()
        try:
            chrom = data[0].replace("chr","")
            start = int(data[1])
            end = int(data[2])
        except:
            print("error: malformed BED line "+str(i))
            sys.exit(1)
        if chrom not in intervaltreeset: 
            intervaltreeset[chrom] = IntervalTree()
        intervaltreeset[chrom].insert(start, end + 1, line)
    bed_fh.close()
    return intervaltreeset

def rgn_to_intervaltreeset(rgn_fh):
    intervaltreeset = dict()
    i = 0
    for line in rgn_fh:
        i += 1
        data = line.rstrip().split()
        try:
            rgnID = data[0]
            chrom = data[1].replace("chr","")
            intervals_str = data[2]
            rgn_length = data[3]
            intervals = rgn_intervals_str_to_list(intervals_str)
        except:
            print("error: malformed RGN line "+str(i))
            sys.exit(1)
        if chrom not in intervaltreeset:
            intervaltreeset[chrom] = IntervalTree()
        for interval in intervals:
            intervaltreeset[chrom].add(interval[0]-1, interval[1]+1,
                                       rgnID)
    rgn_fh.close()
    return intervaltreeset


def rgn_intervals_str_to_list(rgn_intervals_str):
    rgn_intervals_str = rgn_intervals_str.replace("(","")
    rgn_intervals_str = rgn_intervals_str.replace(")","")
    rgn_intervals = rgn_intervals_str.split(",")
    for i in range(len(rgn_intervals)):
        rgn_intervals[i] = rgn_intervals[i].split("..")
        rgn_intervals[i] = [int(rgn_intervals[i][0]),
                            int(rgn_intervals[i][1])]

    return rgn_intervals

if __name__ == "__main__":
    main()
