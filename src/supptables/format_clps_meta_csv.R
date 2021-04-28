
## LIBRARY
library(DescTools)

## PARAM
DOT.RM <- T
N.CLUSTERS <- 11
ADD.ZEROCOUNT <- T

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 4) {
    cat("format_clps_meta_csv.R <sampleped_dir> <in.unformatted.csv>",
        "<ensg_genesymbol_entrezid.tsv> <out.formatted.csv>\n")
    q()
  }
  sampleped.dir <- ARGS[1]
  in.unformatted.csv <- ARGS[2]
  symbol_ensg_entrez.tsv <- ARGS[3]
  out.formatted.csv <- ARGS[4]

  # for each cluster read in sampleped inroot, store ncases and nctrls
  nca <- list()
  nco <- list()
  for (i in 0:(N.CLUSTERS-1)) {
    sampleped.i <- paste0(sampleped.dir, "/cluster_",i,".caco.sampleped")
    ped.i <- read.table(sampleped.i, stringsAsFactors=F)
    nca.i <- nrow(ped.i[ped.i[,6]==2, , drop=F])
    nco.i <- nrow(ped.i[ped.i[,6]==1, , drop=F])
    i_str <- as.character(i)
    nca[[i_str]] <- nca.i
    nco[[i_str]] <- nco.i
  }

  # read unformatted csv to df
  df <- read.csv(in.unformatted.csv, 
                 check.names=F,
                 stringsAsFactors=F)

  # general adjustments to input unformatted dataframe
  df$gene <- gsub("'","",df$gene)

  # read table with symbol/ensg/entrezid
  symbol_ensg_entrez <- read.table(symbol_ensg_entrez.tsv, header=T, sep="\t",
                                   stringsAsFactors=F)

  # replace dots with underscores?
  if (DOT.RM == T) {
    colnames(df) <- gsub("\\.","_",colnames(df))
  }

  # merge results with ensg / enstrez id table, resort by cmh p
  colnames(symbol_ensg_entrez) <- c("gene","ensg","entrez_id")
  df <- merge(symbol_ensg_entrez, df, by="gene")
  df <- df[order(df$p_value), ]

  # derive woolf p-value for each entry
  n.clusters.0 <- N.CLUSTERS - 1
  cluster.nums <- 0:n.clusters.0
  ca.q.cols <- paste0(cluster.nums,"_Qualified Case")
  ca.nq.cols <- paste0(cluster.nums,"_Unqualified Case")
  co.q.cols <- paste0(cluster.nums,"_Qualified Ctrl")
  co.nq.cols <- paste0(cluster.nums,"_Unqualified Ctrl")
  cols.x <- c()
  cols.final <- c('gene','ensg','entrez_id','p_value','estimate','conf_low',
                  'conf_high','woolf_p_value','Qualified Case','Unqualified Case',
                  'Qualified Ctrl','Unqualified Ctrl','%QV+ Case','%QV+ Ctrl')
  for (i in cluster.nums) {
    cols.x.i <- c(paste0(i,"_Qualified Case"),
                  paste0(i,"_Unqualified Case"),
                  paste0(i,"_Qualified Ctrl"),
                  paste0(i,"_Unqualified Ctrl"))
    cols.x <- c(cols.x, cols.x.i)
    cols.final <- c(cols.final, cols.x.i, 
                    paste0(i,"_Fet P")
                   )
  }
  caco.q.nq <- df[, cols.x]
  rownames(caco.q.nq) <- df$gene
  woolf.p <- apply(caco.q.nq, 1, function(x) {
                   return(WoolfTest(array(x,dim=c(2,2,N.CLUSTERS)))$p.value)
                  })
  df[["woolf_p_value"]] <- rep(NA, nrow(df))
  rownames(df) <- df$gene
  for (gene.i in names(woolf.p)) {
    df[gene.i, "woolf_p_value"] <- woolf.p[[gene.i]]
  }

  # if odds ratio is Inf store as Inf, rather than leaving set as NA
  for (col.i in c("estimate","conf_low","conf_high")) {
    df[[col.i]] <- ifelse(is.na(df[[col.i]]), Inf, df[[col.i]])
  }

  # if p-values are NA, store as 1 instead
  for (i in cluster.nums) {
    pval.col.i <- paste0(i,"_Fet P")
    df[[pval.col.i]] <- ifelse(is.na(df[[pval.col.i]]),
                               1,
                               df[[pval.col.i]]
                              )
  }

  # if there are genes with zero-count of qualified variants, add
  # as rows to output dataframe
  if (ADD.ZEROCOUNT == T) {
    gs.missing <- setdiff(symbol_ensg_entrez$gene,
                          df$gene)
    df.qv0 <- subset(symbol_ensg_entrez, gene %in% gs.missing)
    df.qv0$p_value <- rep(1, nrow(df.qv0))
    df.qv0$estimate <- rep(0, nrow(df.qv0))
    df.qv0$conf_low <- rep(0, nrow(df.qv0))
    df.qv0$conf_high <- rep(Inf, nrow(df.qv0))
    df.qv0$woolf_p_value <- rep(1, nrow(df.qv0))
    df.qv0[['Qualified Case']] <- rep(0, nrow(df.qv0))
    df.qv0[['Unqualified Case']] <- rep(0, nrow(df.qv0))
    df.qv0[['Qualified Ctrl']] <- rep(0, nrow(df.qv0))
    df.qv0[['Unqualified Ctrl']] <- rep(0, nrow(df.qv0))
    df.qv0[['%QV+ Case']] <- rep(0, nrow(df.qv0))
    df.qv0[['%QV+ Ctrl']] <- rep(0, nrow(df.qv0))
    nca_tot <- 0
    nco_tot <- 0
    for (inum in cluster.nums) {
      i <- as.character(inum)
      df.qv0[[paste0(i,"_Qualified Case")]] <- rep(0, nrow(df.qv0))
      df.qv0[[paste0(i,"_Unqualified Case")]] <- rep(nca[[i]], nrow(df.qv0))
      df.qv0[[paste0(i,"_Qualified Ctrl")]] <- rep(0, nrow(df.qv0))
      df.qv0[[paste0(i,"_Unqualified Ctrl")]] <- rep(nco[[i]], nrow(df.qv0))
      df.qv0[[paste0(i,"_Fet P")]] <- rep(1, nrow(df.qv0))
      nca_tot <- nca_tot + nca[[i]]
      nco_tot <- nco_tot + nco[[i]]
    }
    df.qv0[['Unqualified Case']] <- rep(nca_tot, nrow(df.qv0))
    df.qv0[['Unqualified Ctrl']] <- rep(nco_tot, nrow(df.qv0))

    # after zerocount df is setup, perform row-bind to bottom of 
    # non-zerocount dataframe already existing
    df <- rbind(df, df.qv0)
  }
   
  # final reordering of columns plus column renaming if needed
  df <- df[, cols.final]

  # write to output csv
  write.csv(df,
            file=out.formatted.csv,
            row.names=F,
            quote=F)

}

if (interactive() == F) {
  Main()
}
