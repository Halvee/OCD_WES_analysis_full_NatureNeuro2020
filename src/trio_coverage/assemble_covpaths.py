

import sys
import os

def main():
    try:
        ARGS = sys.argv[1:]
        covdir = ARGS[0]
        sampleped = ARGS[1]
        cov = ARGS[2]
        covpaths_file = ARGS[3]
    except:
        print("assemble_covpaths.py <covdir> <trios.sampleped> <10x|20x> <covpaths.txt>")
        sys.exit(1)

    # load fid_iid of probands from sampleped
    fid_iid_keep = set([])
    ped_fh = open(sampleped,"r")
    for line in ped_fh:
        data = line.rstrip().split()
        [fid, iid, pid, mid, sex, phe] = data[:6]
        if phe == "2" and pid != "0" and mid != "0":
            fid_iid = fid + "_" + iid
            fid_iid_keep.add(fid_iid)
    ped_fh.close()

    # parse contents of covdir, write files to keep to output file
    covpaths_fh = open(covpaths_file, "w")
    for filename in os.listdir(covdir):
        if os.path.isfile(covdir + "/" + filename):
            filename_short = os.path.basename(filename)
            if filename_short.find(cov) == -1: continue
            file_fid_iid = filename_short.split(".")[0]
            if file_fid_iid in fid_iid_keep:
                filename_full = covdir + "/" + filename
                covpaths_fh.write(filename_full + "\n")
    covpaths_fh.close()
   
    return

if __name__ == "__main__":
    main()
