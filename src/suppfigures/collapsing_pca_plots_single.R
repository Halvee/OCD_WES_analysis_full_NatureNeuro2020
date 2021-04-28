#!/usr/bin/env Rscript

ETHNICITIES <- c("African", "Caucasian", "EastAsian",
                 "Hispanic","MiddleEastern","SouthAsian")
N_CA <- c(18,   925,  5,   40,  91,  9  )
N_CO <- c(2169, 4340, 165, 278, 697, 235)
NONCARRIER_COL <- "#999999"
CARRIER_COL <- "black"

Main <- function(){

  # get ARGS
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 3) {
    cat("collapsing_pca_plots_single.R <cohort.evec>",
        "<PC1,PC2,PC3,PC4..>",
        "<out.pdf>\n")
    q()
  }

  cohort.evec <- ARGS[1]
  pcs.str <- ARGS[2]
  out.pdf <- ARGS[3]

  # make sure out.pdf has .pdf at end, if not call error
  if (grepl(".pdf", out.pdf)==F) {
    cat("ERROR : file ", out.pdf, 
        " appears to be missing a required .pdf extension.")
    q()
  }

  # read evec file
  evec <- ReadEvec(cohort.evec)

  # get pcs
  pcs <- strsplit(pcs.str, ",")[[1]]
  n_pcs <- length(pcs)
  nrows <- length(pcs)/2

  # init output pdf
  pdf(out.pdf,width=7,height=7)
  par(mfrow=c(2,2))
  i <- 1
  while (i < n_pcs) {
    pc_i <- pcs[i]
    j <- i + 1
    pc_j <- pcs[j]
    xmin <- min(evec[[pc_i]])
    xmax <- max(evec[[pc_i]])
    ymin <- min(evec[[pc_j]])
    ymax <- max(evec[[pc_j]])
    evec.ca <- subset(evec, 
                      grepl("ocd", rownames(evec))
                     )
    evec.co <- subset(evec,
                      grepl("ocd", rownames(evec))==F
                     ) 
   
    # pc plot
    plot(evec.co[[pc_i]], evec.co[[pc_j]],
         xlim=c(xmin,xmax), ylim=c(ymin,ymax),
         xlab=pc_i, ylab=pc_j,
         main=paste0(pc_j," vs. ",pc_i),
         pch=20, col="#999999"
        )
    points(evec.ca[[pc_i]], evec.ca[[pc_j]], 
           pch=20, col='#000000'
        )

    # increment i
    i <- i + 2
  }

  # save pdf file
  dev.off()

}

ReadEvec <- function(evecfile) {
  evec <- read.table(evecfile, row.names=1, stringsAsFactors=F)
  evec[,ncol(evec)] <- NULL
  colnames(evec) <- paste0("PC",1:10)
  return(evec)
}

if (interactive() == F) {
  Main()
}
