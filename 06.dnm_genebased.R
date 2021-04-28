#!/usr/bin/env Rscript

## DATA
GENE_MU_TSV<-"results/trio_coverage/OCDfams_2019.trios.n1.exp_dnm_rates.tsv"
OUTROOT<-"results/dnm_genebased/OCD_JHU_Cappi_n771"
ANNOTS <- c("dmg", "lof")
TRIO_COHORTS<-c("ocd_jhu_2020",
                "ocd_cappi_2019")
NTRIOS<-c(587,
          184)
TRIO_CDNM_FILES<-c(
                   "results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.inkits.cdnm",
                   "results/dnm_casectrl/Cappi_etal_2019_OCD.DNMs.ppn2hdiv.inkits.cdnm"
                  )
GENESET_FILE <- "results/input/CCDSr20.geneset"
GENES_BED <- "results/input/addjusted.CCDS.genes.index.r20.hg19.bed"

## PARAM
FUNC_SETS <- list("dmg"=c("non","splice","frameshift","misD"),
                  "dmg_snv"=c("non","splice","misD"),
                  "lof"=c("non","splice","frameshift"),
                  "lof_snv"=c("non","splice"))
Main <- function() {
  
  # read gene mutation rate table
  gene.mu <- read.table(GENE_MU_TSV,header=T,sep="\t",row.names=1)

  # damaging snvs : LOF snvs + misD
  gene.mu$dmg_snv <- rowSums(gene.mu[,c("lof_snv","misD")])

  # read test geneset, subset mutation rate table on it
  test.geneset <- scan(GENESET_FILE, what=character(), quiet=T)
  gene.mu <- gene.mu[rownames(gene.mu) %in% test.geneset, , drop=F]

  # define total number of tests as the total number of genes in the 
  # gene-based mutation rate table
  ntests <- nrow(gene.mu)

  # read dnm callsets and store to memory
  ntrios_total <- 0
  cdnm <- NULL
  cdnm_cohorts <- list()
  for (i in 1:length(TRIO_COHORTS)) {
    cohort_i <- TRIO_COHORTS[i]
    ntrios_i <- as.numeric(NTRIOS[i])
    cdnm_i <- read.table(TRIO_CDNM_FILES[i], stringsAsFactors=F)
    ntrios_total <- ntrios_total + ntrios_i
    if (is.null(cdnm)) {
      cdnm <- cdnm_i
    } else {
      cdnm <- rbind(cdnm, cdnm_i)
    }
    cdnm_cohorts[[cohort_i]] <- cdnm_i
  }

  # iterate through exome-wide analysis for each single input annotation set
  for (i in 1:length(ANNOTS)) {

    # get annotation string
    annot_i <- ANNOTS[i]

    # based on annotation, get effect classifications to keep
    effs_keep_i <- FUNC_SETS[[annot_i]]

    # init output df
    df <- data.frame(gene=rownames(gene.mu), 
                     rate.per.sample=gene.mu[[annot_i]],
                     count.exp=gene.mu[[annot_i]] * ntrios_total * 2,
                     count.obs=rep(0,nrow(gene.mu)), 
                     p.value=rep(1,nrow(gene.mu)))
    rownames(df) <- df$gene

    # only keep variants in gene.mu with applicable effs.keep
    cdnm_i <- cdnm[(cdnm[,7] %in% rownames(gene.mu)) &
                   (cdnm[,8] %in% effs_keep_i), , drop=F]
    genes_i <- cdnm_i[,7]
    genes_i <- unique(genes_i)
    rep.genes_i <- unique(sort(genes_i[duplicated(genes_i)==T]))
    
    # get results df
    for (gene_j in genes_i) {
      count.exp <- df[gene_j, "count.exp"]
      count.obs <- nrow(cdnm_i[cdnm_i[,7] == gene_j, , drop=F])
      pval <- ppois(count.obs - 1, lambda=count.exp, lower.tail=F)
      df[gene_j, "count.obs"] <- count.obs
      df[gene_j,"p.value"] <- pval
    }
    df <- df[order(df$p.value), ]

    # write all gene-based results, all recurrent hit results
    write.table(df, file=paste0(OUTROOT, ".", annot_i, ".all.tsv"),
                row.names=F, quote=F, sep="\t")
    df.recurr <- subset(df, count.obs > 1)
    write.table(df.recurr, file=paste0(OUTROOT, ".", annot_i, ".recurr.tsv"),
                row.names=F, quote=F, sep="\t")

    genes <- read.table(GENES_BED, header=F, sep="\t", stringsAsFactors=F)
    colnames(genes)<-c("chrom",   "start",   "end",     "gene")
    genes <- genes[ duplicated(genes$gene) == F, ]
    genes$midpt <- (genes$end + genes$start) / 2
    genes <- genes[, c("gene","chrom","midpt")]
    genes$chrom <- gsub("chr","", genes$chrom)
    df <- merge(df, genes, by="gene")
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

    df$logp <- -log10( df$p.value )


    pdf(paste0(OUTROOT,".",annot_i,".mannhattan.pdf"))

    df <- df[order(df$midpt.rank), ]

    add.to.plot <- F
    for (chrom in seq(1,23,by=1)) {
      if (chrom %% 2 == 1) {
        colchrom = rgb(0.0,0.3,0.2,0.5)
      } else {
        colchrom = rgb(0.0,0.0,0.8,0.5)
      }
     
      df.c <- df[ df$chrom == chrom, ]
      if (chrom == 1) {
        plot(df.c$midpt.rank, df.c$logp, type="s",
             xlim=c(1,max(df$midpt.rank)), xaxt='n',
             ylim=c(0,6),
             col=colchrom,
             xlab="", ylab="-log10(P)", main="")
      
      } else {
        lines(df.c$midpt.rank, df.c$logp, type="s",
              col=colchrom)
      }
    }

    axis(1, at=chroms.ticks, labels=seq(1,23,by=1), cex.axis=0.6, las=3)

    abline(h= -log10( 0.05 / (ntests*length(ANNOTS)) ), col=rgb(0.7,0,0,0.5))
    abline(h= -log10( 0.0001 ), col=rgb(0,0,0,0.5))

    df.p <- df[ df$logp > 2 & df$count.obs > 1, , drop=F]
    text(df.p$midpt.rank, df.p$logp+0.1, label=df.p$gene)

    dev.off()

  }

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

if (interactive() == F) {
  Main()
}
