

library(data.table)
COLS <- c("red","blue","green","purple","brown","goldenrod")

## PARAM
ABBREV.MODE <- F

ARGS <- commandArgs(trailingOnly=T)
if (length(ARGS) != 5) {
  cat("plot_cnv_gene.R <in.PCA_NORM_Z_SCORES.mat.txt>",
      "<capture.bed> <chr:bp1-bp2> <iid1,..,iidN> <out.pdf>\n")
  q()
}
in.pca_norm_z.mat.txt <- ARGS[1]
capture.bed <- ARGS[2]
interval.str<- ARGS[3]
iids.str <- ARGS[4]
out.pdf <- ARGS[5]

# get chrom, bp1, bp2
interval <- strsplit(interval.str, ":")[[1]]
chrom <- interval[1]
startend <- strsplit(interval[2], "-")[[1]]
bp1 <- as.numeric(startend[1])
bp2 <- as.numeric(startend[2])

# read cnv.gene bed file, get loci to keep
capture <- fread(capture.bed, sep="\t", 
                 data.table=F, stringsAsFactors=F)
colnames(capture)[1:3] <- c("CHR","BP1","BP2")
capture$LOCUS <- paste0(capture$CHR, 
                             ":", capture$BP1,
                             "-", capture$BP2)
capture <- subset(capture, (CHR==chrom) & 
                       (BP1 > bp1) &
                       (BP1 < bp2))
gene.loci <- capture$LOCUS
gene.loci <- gene.loci[duplicated(gene.loci)==F]

# read files
mat<-fread(in.pca_norm_z.mat.txt, header=T, 
           select=c("Matrix",gene.loci),
           sep="\t", data.table=F)
rownames(mat) <- mat[,1]
mat[,1] <- NULL

# rename matrix columns
mat.colnames<-colnames(mat)
positions <- c()
for (i in 1:ncol(mat)) {
  mat.col.i <- colnames(mat)[i]
  interval.str <- strsplit(mat.col.i,":")[[1]]
  interval <- strsplit(interval.str[2], "-")[[1]]
  positions <- c(positions, as.numeric(interval[1]))
}

# plot per sample
iids.ca <- strsplit(iids.str, ",")[[1]]
iids.co <- setdiff(rownames(mat), iids.ca)
pdf(out.pdf)
if (ABBREV.MODE==F) { 
  for (i in 1:length(iids.co)) {
    if (i==1) {
      plot(positions, mat[iids.co[i], ], 
           xlab='position', ylab='pca-norm rd z',
           ylim=c(min(mat), max(mat)),
           col='black')
    } else {
      lines(positions, mat[iids.co[i], ], col='black')
    }
  }
  for (i in 1:length(iids.ca)) {
    lines(positions, mat[iids.ca[i], ], col=COLS[i]) 
  }
} else {
   plot(positions, apply(mat[iids.co, ,drop=F], 2, median), 
        xlab='position', ylab='pca-norm rd z',
        ylim=c(min(mat), max(mat)),
        col='black', type='p', pch='.')
   lq <- apply(mat[iids.co, ,drop=F],2,
               function(x){return(quantile(x)[["25%"]])})
   uq <- apply(mat[iids.co, ,drop=F],2,
               function(x){return(quantile(x)[["75%"]])})
   for (i in 1:length(iids.ca)) { 
     points(positions, mat[iids.ca[i], ], col=COLS[i], pch='.') 
   }
   points(positions, as.numeric(lq),
          col='black', pch='.') 
   points(positions, 
          as.numeric(uq),
          col='black',pch='.') 
}
legend("bottomright", legend=iids.ca,
       col=COLS[1:length(iids.ca)], lty=1, cex=0.8)
dev.off()
