import os
import sys
import argparse
import gzip

CHRS_INCLUDE=set(['1','2','3','4','5','6','7','8','9','10',
                  '11','12','13','14','15',
                  '16','17','18','19','20',
                  '21','22','X'])

def main():
    args = parse_args()
    print("loading genders from sampleped file...")
    sampleped_fh = open(args.sampleped_file, "r")
    genders = dict()
    for line in sampleped_fh:
        data = line.rstrip().split()
        genders[ data[1] ] = int( data[4] )
    sampleped_fh.close()

    print("loading parent .lis file...")
    rgn_cov_tree = load_p_lis(args.parent_lis_filename)
    if args.samplecovpaths_file == "stdin":
        rgncov_filepaths_fh = sys.stdin
    else:
        rgncov_filepaths_fh = open(args.samplecovpaths_file, "r")
    print("loading rgn filepaths...")
    rgncov_filepaths = load_rgncov_filepaths( rgncov_filepaths_fh, genders )
   
    for gender in (1,2):
        for rgncov_filepath in rgncov_filepaths[ gender ]:
            print("loading "+rgncov_filepath+" to coverage tree..")
            rgn_cov_tree = load_rgncov( rgncov_filepath, gender, rgn_cov_tree )
    print("writing per-nucleotide coverage stats to output file")
    write_cov_samples_per_pos(rgn_cov_tree, 
                              args.parent_lis_filename,
                              args.out_filename,
                              full_loc_info=args.full_loc_info,
                              buff_size=args.buff_size)
    return

def parse_args():
    ''' load user-defined arguments and options '''
    desc='''for each position in parent region file, return the number of samples
            from list of paths to coverage files where coverage is sufficient 
            for de novo mutation analysis '''
    parser = argparse.ArgumentParser(prog="cov_samples_per_pos", 
                                     description=desc)
    parser.add_argument('--includeY', action="store_true", dest="includeY", 
                        default=False, help='include Y chromosome genes')
    parser.add_argument('--full-loc-info',action="store_true",default=False,
                        help="in output, write cols 1-3 for loc info.")
    parser.add_argument('--buff-size',action="store",type=int,default=50000,
                        help="size of buffer in writing to output file")
    parser.add_argument('parent_lis_filename', action="store", 
                        help='.lis file with all nucleotides to assess for coverage')
    parser.add_argument('samplecovpaths_file', action="store", 
                         help='table with IID as col1 and col2 as path to ' + \
                              'per-trio coverage files representing assessable sites')
    parser.add_argument('sampleped_file',action="store",
                        help='ped file with needed gender info')
    parser.add_argument('out_filename', action="store", help='name of analysis output file')
    args = parser.parse_args(sys.argv[1:])
    return args

def load_rgncov_filepaths(rgncov_filepaths_fh, genders_dict):
    ''' make sure rgncov_filepath points to an existing file and then add it to list of 
        files to assess, otherwise halt program '''
    rgncov_filepaths = {1:[], 2:[]}
    for line in rgncov_filepaths_fh:
        rgncov_filepath = line.rstrip()
        if os.path.isfile(rgncov_filepath) == False:
            print("rgncov file "+rgncov_filepath+" does not exist.  Exiting...")
            sys.exit(1)
        filename_base = os.path.basename(rgncov_filepath)
        fid_iid = filename_base.split(".")[0]
        [fid, iid] = fid_iid.split("_")
        gender = genders_dict[ iid ]
        rgncov_filepaths[ gender ].append( rgncov_filepath )
    return rgncov_filepaths


def load_p_rgn(p_rgn_fh):
    ''' load parent genomic intervals into tree of the following structure :
        /\
          gene
           \
            chrom
             \
              pos
    '''
    rgn_cov_tree = {}
    i = 0
    for line in p_rgn_fh:
        i += 1
        data = line.rstrip().split()
        try:
            gene,chrom,intervals_str = data[:3]
            intervals = load_intervals_str(intervals_str, i)
        except:
            print("malformed rgn file entry at line " + str(i) + ".  Exiting...")
            sys.exit(1)
        if chrom not in rgn_cov_tree: rgn_cov_tree[chrom] = {}
        if gene not in rgn_cov_tree[chrom]: rgn_cov_tree[chrom][gene] = {}
        for interval in intervals:
            for pos in xrange(interval[0], interval[1]+1):    
                rgn_cov_tree[chrom][gene][pos] = 0
        print(i)

    p_rgn_fh.close()
    return rgn_cov_tree

def load_p_lis(parent_lis_filename, skip_header=True):
    ''' load parent genomic sites from .lis file into tree of 
        the following structure :
        /\
          gene
           \
            chrom
             \
              pos
    '''

    p_lis_fh = open(parent_lis_filename, "r")
    rgn_cov_tree = {}
    i = 0
    for line in p_lis_fh:
        i += 1
        if i == 1 and skip_header == True: continue
        data = line.rstrip().split()
        try:
            gene,chrom,pos = data[:3]
            pos = int(pos)
        except:
            print("malformed .lis file entry at line " + str(i) + ".  Exiting...")
            sys.exit(1)
        if chrom not in rgn_cov_tree: rgn_cov_tree[chrom] = {}
        if gene not in rgn_cov_tree[chrom]: rgn_cov_tree[chrom][gene] = {}
        if pos not in rgn_cov_tree[chrom][gene]:
            rgn_cov_tree[chrom][gene][pos] = 0
    
    p_lis_fh.close()
    return rgn_cov_tree

def load_intervals_str(intervals_str, i = -1):
    ''' convert rgn formatted intervals to a list of interval start-stops '''
    intervals_str = intervals_str[1:-1]
    if len(intervals_str) == 0:
        return []
    intervals = intervals_str.split(",")
    for j in xrange(len(intervals)):
        intervals[j] = intervals[j].split("..")
        intervals[j] = [int(intervals[j][0]), int(intervals[j][1])]
        #try:
        #    intervals[j] = [int(intervals[j][0]), int(intervals[j][1])]
        #except:
        #    print("malformed intervals at line " + str(i) + ".  Exiting...")
        #    sys.exit(1)
    return intervals

def load_rgncov(rgncov_filename, rgncov_gender, 
                rgn_cov_tree, chrom_gene = True):
    ''' count sufficiently covered nucleotides from parent region file using trio-level
        rgn-cov file, store to rgn_cov_tree, return rgn_cov_tree '''
    global CHRS_INCLUDE
    if rgncov_filename.find(".gz") != -1:
        rgncov_fh = gzip.open(rgncov_filename, "rb")
    else:
        rgncov_fh = open(rgncov_filename, "r")
    i = 0

    # data struct to keep track of repeated chrom/pos passed in rgn file
    cov_chrpos=dict()

    for line in rgncov_fh:
        i += 1

        data = line.rstrip().split()
        if chrom_gene == True:
            chrom,gene,intervals_str = data[:3]
        else:
            gene,chrom,intervals_str = data[:3]
        intervals = load_intervals_str(intervals_str, i)
        if chrom not in CHRS_INCLUDE: continue
        if (chrom=="X" or chrom=="Y") and rgncov_gender==1:
            n_copies = 2
        else:
            n_copies = 2
        gene_node = rgn_cov_tree[chrom][gene]
        if chrom not in cov_chrpos: cov_chrpos[chrom]=set()
        for interval in intervals:
            for pos in range(interval[0], interval[1]+1):
                if pos not in cov_chrpos[chrom]:
                    gene_node[pos] += n_copies
                    cov_chrpos[chrom].add(pos)
    rgncov_fh.close()
    return rgn_cov_tree

def write_cov_samples_per_pos(rgn_cov_tree, 
                              parent_lis_filename,
                              out_filename, delim="\t",
                              buff_size=50000,
                              full_loc_info=False,
                              skip_header=True):
    ''' print number of samples sufficiently covered at each
        gene-chrom-pos in input parent rgn file '''
    out_fh = open(out_filename, "w")
    out_fh.write("")
    out_fh.close()
    
    p_lis_fh = open(parent_lis_filename, "r")
    b = 0

    i = 0
    b = 0
    out_fh = open(out_filename, "a")
    for line in p_lis_fh:
        i += 1
        if i == 1 and skip_header == True: continue
        data = line.rstrip().split()
        gene,chrom,pos = data[:3]
        pos = int(pos)
        n_samples_cov = rgn_cov_tree[chrom][gene][pos]
        if full_loc_info == False:
            out_row_str = str( n_samples_cov )
        else:
            out_row = [gene,chrom,pos,n_samples_cov]
            out_row_str = delim.join([str(item) for item in out_row])
        out_fh.write(out_row_str + "\n")
        b += 1
        if b == buff_size:
            out_fh.close()
            out_fh = open(out_filename, "a")
    out_fh.close()
    return    

if __name__ == '__main__':
    main()
