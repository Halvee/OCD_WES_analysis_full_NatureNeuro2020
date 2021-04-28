

## DATA
GENES <- c("COL4A3","FN1","ZMYM2","CHD8","STARD13","INO80","NWD1","SCUBE1","ADIPOR2")

Main <- function() {

  ## get ARGS
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 3) {
    cat("plot_extTADA.R <genes.bed> <exttada.res.tsv> <out.pdf>\n")
    q()
  }
  genes_bed <- ARGS[1]
  exttada_res_tsv <- ARGS[2]
  out_pdf <- ARGS[3]

  df<-read.table(exttada_res_tsv, header=T, sep="\t", stringsAsFactors=F)
  genes <- read.table(genes_bed, header=F, sep="\t", stringsAsFactors=F)
  colnames(genes)<-c("chrom",   "start",   "end",     "Gene")
  genes <- genes[ duplicated(genes$Gene) == F, ]
  genes$midpt <- (genes$end + genes$start) / 2
  genes <- genes[, c("Gene","chrom","midpt")]
  genes$chrom <- gsub("chr","", genes$chrom)
  df <- merge(df, genes, by="Gene")
  df$chrom <- gsub("X",23,df$chrom)
  df$chrom <- as.numeric(df$chrom)
  df <- df[ order(df$chrom, df$midpt), ]
  df$midpt.rank <- rep(0, nrow(df))
  x<-1
  df$chrom <- as.numeric(df$chrom)
  chroms <- c()
  chroms.ticks <- c()
  for (chrom.i in seq(1,23,by=1)) {
    df.c <- subset(df, chrom==chrom.i)
    df.c <- df.c[ order(df.c$midpt), ]
    cids <- rownames(df.c)
    df[ cids, "midpt.rank" ] <- seq(x,
                                     (x+nrow(df.c)-1),
                                     by=1)
    df.c[ cids, "midpt.rank" ] <- seq(x,
                                      (x+nrow(df.c)-1),
                                      by=1)
    midpt.mid <- (min(df.c$midpt.rank) + max(df.c$midpt.rank)) / 2

    chroms <- c(chroms, x)
    chroms.ticks <- c(chroms.ticks, midpt.mid)
    x <- x + nrow(df.c) +1
  }
  chroms <- c(chroms, x)

  df$logq <- -log10( df$qvalue )


  pdf(out_pdf)

  df <- df[order(df$midpt.rank), ]

  add.to.plot <- F
  for (chrom in seq(1,23,by=1)) {
    if (chrom %% 2 == 1) {
      colchrom = rgb(0.0,0.3,0.2,1)
    } else {
      colchrom = rgb(0.0,0.0,0.8,1)
    }
   
    df.c <- df[ df$chrom == chrom, ]
    if (chrom == 1) {
      plot(df.c$midpt.rank, df.c$logq, type="p",
           xlim=c(1,max(df$midpt.rank)), xaxt='n',
           ylim=c(0,1.025),
           col=colchrom,
           xlab="", ylab="-log10(Q)", main="", pch=20)
      #plot(df.c$midpt.rank, df.c$logq, type="s",
      #     xlim=c(1,max(df$midpt.rank)), xaxt='n',
      #     ylim=c(0,1.2),
      #     col=colchrom,
      #     xlab="", ylab="-log10(Q)", main="")
    
    } else {
      points(df.c$midpt.rank, df.c$logq, col=colchrom, pch=20)
      #lines(df.c$midpt.rank, df.c$logq, type="s",
      #      col=colchrom)
    }
  }

  axis(1, at=chroms.ticks, labels=seq(1,23,by=1), cex.axis=0.6, las=3)

  abline(h= -log10( 0.3 ), col="black")
  abline(h= -log10( 0.1 ), col="firebrick4")
  text(500, 0.55, label="Q < 0.3", col='black')
  text(500, 1.03, label="Q < 0.1", col="firebrick4")

  df.p <- df[df$Gene %in% GENES, , drop=F]
  text(df.p$midpt.rank, df.p$logq+0.05, label=df.p$Gene)

  dev.off()

}

MapCoords <- function(df.x, coord.df) {
  coords.new <- c()
  for (i in 1:nrow(coord.df)) {
    midpt <- (coord.df[i, "end"] + coord.df[i, "start"]) / 2
    df.x.c <- df.x[ df.x$chrom == coord.df[i, "chrom"], ]
    x <- which.min(abs(df.x.c$midpt-midpt))
    midpt.rank.i <- df.x.c$midpt.rank[x]
    cat(midpt,x,midpt.rank.i,"\n")
    coords.new <- c(coords.new, midpt.rank.i)
  }
  return(coords.new)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

if (interactive() == F) {
  Main()
}
