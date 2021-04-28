
library(data.table)

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T) 
  if (length(ARGS) < 4) {
    cat("get_scores_cdnm.R <scores.lis>",
        "<exclude.varlist>",
        "<in.1.cdnm> <out.1.cdnm>", "..",
        "<in.N.cdnm> <out.N.cdnm>\n")
    q()
  }

  scores <- fread(ARGS[1], header=T, data.table=F)
  exclude.varlist <- scan(ARGS[2], what=character())
  i <- 3
  while (i < length(ARGS)) {
    in.cdnm <- ARGS[i]
    out.cdnm <- ARGS[(i+1)]
    print(in.cdnm)
    cdnm.i <- read.table(in.cdnm, stringsAsFactors=F)

    # ensure that variant ids in exclude.varlist are filtered from cdnm
    # presumably varlist is gnomAD POPMAX >= 0.0005
    cdnm.i <- cdnm.i[(cdnm.i[,4] %in% exclude.varlist)==F, , drop=F] 

    # subset scores table on genes that are in cdnm file
    genes <- unique(cdnm.i[,7])
    scores.g <- scores[scores[,1] %in% genes, , drop=F]
    
    # iterate through every row in cdnm, if variant is missense look up
    # the annotation according to polyphen humdiv, misD is score >= 0.957
    for (j in 1:nrow(cdnm.i)) {
      if (grepl("mis", cdnm.i[j, 8]) == T) {
        print(cdnm.i[j,])
        chrom <- cdnm.i[j,2]
        pos <- cdnm.i[j,3]
        gene <- cdnm.i[j,7]
        alt <- cdnm.i[j,6]
        rowname.j <- paste(gene, chrom, pos, sep=":")
        scores.g.j <- scores.g[(scores.g[,1]==gene) &
                               (scores.g[,3]==pos), , drop=F]
        if (nrow(scores.g.j) != 1) { next }
        val <- scores.g.j[1, alt]
        if ((val >= 0.957)) {
          cdnm.i[j,8] <- "misD"
        } else if ((val < 0.957) & (val > 0.452)) {
          cdnm.i[j,8] <- "misP" 
        } else if ((val <= 0.452) & (val >= 0)) {
          cdnm.i[j,8] <- "misB"
        } else {
          cdnm.i[j,8] <- "misU"
        }
        write.table(cdnm.i, file=out.cdnm, row.names=F,
                    col.names=F, quote=F, sep="\t")
      }
    }
    i <- i + 2 
  }
}

if (interactive() == F) {
  Main()
}
