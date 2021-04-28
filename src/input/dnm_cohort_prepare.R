#!/usr/bin/env Rscript

TRIO.FIDS.PRUNE <- c(

                     # treated like a trio but actually are quartets
                     "ocda","ocdg","ocdk","ocdv"

                     )

QUARTET.FIDS.PRUNE <- c()

SAMPLEVAR.PRUNE <- c(

                     # No data returned in sanger validation. In concert
                     # with the call being suboptimal, have to remove.
                     "ocd4005006287aab1:2-225449678-G-A"

                    )

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 3) {
    cat("analysis_prepare.R <triosquartets.sampleped>",
        "<triosquartets.callset.dnm.csv> <outroot>\n")
    q()
  }

  # get ARGS
  triosquartets.sampleped <- ARGS[1]
  triosquartets.callset.dnm.csv <- ARGS[2]
  outroot <- ARGS[3]

  # read input files
  ped <- read.table(triosquartets.sampleped,stringsAsFactors=F)
  colnames(ped) <- c("FID","IID","PID","MID","SEX","PHE",
                     "SEQTYPE","EXOMEPREPKIT")
  callset <- read.csv(triosquartets.callset.dnm.csv, stringsAsFactors=F)
  
  # add samplevar column to callset
  callset$samplevar <- paste(callset$Sample.Name, callset$Variant.ID, sep=":")

  # manually prune failing samples/vars
  callset <- subset(callset,
                    (samplevar %in% SAMPLEVAR.PRUNE)==F)

  # get trios and quartet FIDs
  fid.counts <- table(ped$FID)
  trio.fids <- names(fid.counts[fid.counts==3])
  quartet.fids <- names(fid.counts[fid.counts>3])
  trio.iids <- ped[(ped$FID %in% trio.fids) & (ped$PHE==2), "IID"]
  quartet.iids <- ped[(ped$FID %in% quartet.fids) & (ped$PHE==2), "IID"]
  proband.iids <- unique(c(trio.iids,quartet.iids))
  cat("Total trios (pre-filtering) :", length(trio.fids), "\n")
  cat("Total quartets (pre-filtering) :", length(quartet.fids), "\n")
  
  # ID probands (and corresponding FIDs) which have > 5 DNM calls
  dnm.counts <- sort(table(callset$Sample.Name))
  dnm.count.outliers.iids <- names(dnm.counts[dnm.counts>5])
  dnm.count.outliers.fids <- unique(ped[ped$IID %in% dnm.count.outliers.iids,
                                    "FID"])
  dnm.count.outliers.iids.full <- unique(ped[ped$FID %in%
                                         dnm.count.outliers.fids, "IID"])

  # remove IIDs and FIDs where at least one proband has an outlier dnm count
  trio.fids <- setdiff(trio.fids, dnm.count.outliers.fids)
  quartet.fids <- setdiff(quartet.fids, dnm.count.outliers.fids)
  trio.iids <- setdiff(trio.iids, dnm.count.outliers.iids)
  quartet.iids <- setdiff(quartet.iids, dnm.count.outliers.iids)
  callset <- subset(callset, (Sample.Name %in% dnm.count.outliers.iids.full)==F)
  cat("N probands in trio or quartet fams with n_dnm>5 :",
      length(dnm.count.outliers.iids), "\n")
  cat("N fams removed (proband or sib with n_dnm>5) :",
      length(dnm.count.outliers.fids), "\n")
  cat("N trio fams remaining :", length(trio.fids), "\n")
  cat("N quartet fams remaining :", length(quartet.fids), "\n")

  # remove trios and quartets which are manually labelled for pruning
  iids.prune <- unique(ped[(ped$PHE == 2) & (ped$FID %in% 
                           c(TRIO.FIDS.PRUNE, QUARTET.FIDS.PRUNE)), "IID"])
  trio.fids <- setdiff(trio.fids, TRIO.FIDS.PRUNE)
  quartet.fids <- setdiff(quartet.fids, QUARTET.FIDS.PRUNE)
  callset <- subset(callset, (Sample.Name %in% iids.prune)==F)
  cat("Total trios labelled for manual pruning :", 
      length(TRIO.FIDS.PRUNE), "\n")
  cat("Total quartets labelled for manual pruning :", 
      length(QUARTET.FIDS.PRUNE), "\n")
  cat("N trio fams remaining :", length(trio.fids), "\n")
  cat("N quartet fams remaining :", length(quartet.fids), "\n")

  # create dnm file df
  chroms <- c()
  pos <- c()
  refs <- c()
  alts <- c()
  for (i in 1:nrow(callset)) {
    varid <- callset[i,"Variant.ID"]
    varid.info <- strsplit(varid,"-")[[1]]
    chrom.i <- varid.info[1]
    pos.i <- varid.info[2]
    ref.i <- varid.info[3]
    alt.i <- varid.info[4]
    chroms <- c(chroms, chrom.i)
    pos <- c(pos, pos.i)
    refs <- c(refs, ref.i)
    alts <- c(alts, alt.i)
  }
  dnm <- callset[,c("Sample.Name","Variant.ID","samplevar")]
  dnm$CHR <- chroms
  dnm$POS <- pos
  dnm$REF <- refs
  dnm$ALT <- alts
  dnm <- dnm[,c("Sample.Name","CHR","POS","Variant.ID","REF","ALT","samplevar")]
  dnm <- dnm[order(dnm$CHR, dnm$POS),]

  # split ped, dnm and callset into trios and quartets
  ped <- subset(ped, FID %in% c(trio.fids, quartet.fids))
  ped.t <- subset(ped, FID %in% trio.fids)
  ped.q <- subset(ped, FID %in% quartet.fids)
  callset <- subset(callset, (Sample.Name %in% c(trio.iids, quartet.iids)))
  callset.t <- subset(callset, Sample.Name %in% trio.iids)
  callset.q <- subset(callset, Sample.Name %in% quartet.iids)
  dnm <- subset(dnm, (Sample.Name %in% c(trio.iids, quartet.iids)))
  dnm.t <- subset(dnm, Sample.Name %in% trio.iids)
  dnm.q <- subset(dnm, Sample.Name %in% quartet.iids)

  # get subset of calls from callset where jointcov >= 20x
  callset.20 <- subset(callset,
                       (DP.Bin >= 20) & 
                       (DP.Bin..mother. >= 20) &
                       (DP.Bin..father. >= 20) )
  samplevar.20 <- callset.20$samplevar
  callset.t.20 <- subset(callset.t, samplevar %in% samplevar.20)
  callset.q.20 <- subset(callset.q, samplevar %in% samplevar.20)
  dnm.20 <- subset(dnm, samplevar %in% samplevar.20)

  # print some call stats to stdout
  n.trios <- length(trio.fids)
  n.quartets <- length(quartet.fids)
  n.t.probands <- n.trios
  n.q.probands <- length(ped[(ped$PHE == 2) & (ped$FID %in% quartet.fids),
                             "IID"])
  n.t.snv <- nrow(callset.t[callset.t$Variant.Type == "snv",,drop=F])
  n.t.indel <- nrow(callset.t[callset.t$Variant.Type == "indel",,drop=F])
  n.t.20.snv <- nrow(callset.t.20[callset.t.20$Variant.Type == "snv",
                                  ,drop=F])
  n.t.20.indel <- nrow(callset.t.20[callset.t.20$Variant.Type == "indel",
                                    ,drop=F])
  n.t.snv <- nrow(callset.t[callset.t$Variant.Type == "snv",,drop=F])
  n.t.indel <- nrow(callset.t[callset.t$Variant.Type == "indel",,drop=F])
  n.t.20.snv <- nrow(callset.t.20[callset.t.20$Variant.Type == "snv",
                                  ,drop=F])
  n.t.20.indel <- nrow(callset.t.20[callset.t.20$Variant.Type == "indel",
                                    ,drop=F])
  n.q.snv <- nrow(callset.q[callset.q$Variant.Type == "snv",,drop=F])
  n.q.indel <- nrow(callset.q[callset.q$Variant.Type == "indel",,drop=F])
  n.q.20.snv <- nrow(callset.q.20[callset.q.20$Variant.Type == "snv",
                                  ,drop=F])
  n.q.20.indel <- nrow(callset.q.20[callset.q.20$Variant.Type == "indel",
                                    ,drop=F])
  n.q.snv <- nrow(callset.q[callset.q$Variant.Type == "snv",,drop=F])
  n.q.indel <- nrow(callset.q[callset.q$Variant.Type == "indel",,drop=F])
  n.q.20.snv <- nrow(callset.q.20[callset.q.20$Variant.Type == "snv",
                                  ,drop=F])
  n.q.20.indel <- nrow(callset.q.20[callset.q.20$Variant.Type == "indel",
                                    ,drop=F])

  cat("Total SNVs (trios, 10x) :", n.t.snv, "\n")
  cat("Total indels (trios, 10x) :", n.t.indel, "\n")
  cat("Mean SNVs per proband (trios, 10x) :", n.t.snv/n.t.probands,"\n")
  cat("Mean indels per proband (trios, 10x) :",n.t.indel/n.t.probands,"\n") 
  cat("Total SNVs (trios, 20x) :", n.t.20.snv, "\n")
  cat("Total indels (trios, 20x) :", n.t.20.indel, "\n")
  cat("Mean SNVs per proband (trios, 20x) :", 
      n.t.20.snv/n.t.probands,"\n")
  cat("Mean indels per proband (trios, 20x) :",
      n.t.20.indel/n.t.probands,"\n") 
  cat("Total SNVs (quartets, 10x) :", n.q.snv, "\n")
  cat("Total indels (quartets, 10x) :", n.q.indel, "\n")
  cat("Mean SNVs per proband (quartets, 10x) :", 
      n.q.snv/n.q.probands,"\n")
  cat("Mean indels per proband (quartets, 10x) :",
      n.q.indel/n.q.probands,"\n") 
  cat("Total SNVs (quartets, 20x) :", n.q.20.snv, "\n")
  cat("Total indels (quartets, 20x) :", n.q.20.indel, "\n")
  cat("Mean SNVs per proband (quartets, 20x) :", 
      n.q.20.snv/n.q.probands,"\n")
  cat("Mean indels per proband (quartets, 20x) :",
      n.q.20.indel/n.q.probands,"\n") 








  # remove samplevar columns
  dnm$samplevar <- NULL  
  dnm.20$samplevar <- NULL 
  callset.t$samplevar <- NULL
  callset.q$samplevar <- NULL
  callset.t.20$samplevar <- NULL
  callset.q.20$samplevar <- NULL

  # write pruned samplepeds and dnm callsets to file
  write.table(ped.t, file=paste0(outroot,".trios.sampleped"),
              row.names=F, sep="\t", quote=F, col.names=F)
  write.table(ped.q, file=paste0(outroot,".quartets.sampleped"),
              row.names=F, sep="\t", quote=F, col.names=F)
  write.table(ped, file=paste0(outroot,".triosquartets.sampleped"),
              row.names=F, sep="\t", quote=F, col.names=F)
  write.csv(callset.t, file=paste0(outroot,".trios.DNMs.10x.csv"),
            row.names=F, quote=F)
  write.csv(callset.q, file=paste0(outroot,".quartets.DNMs.10x.csv"),
            row.names=F, quote=F)
  write.table(dnm, file=paste0(outroot, ".triosquartets.10x.dnm"),
              col.names=F, row.names=F, quote=F, sep="\t")
  write.csv(callset.t.20, file=paste0(outroot,".trios.DNMs.20x.csv"),
            row.names=F, quote=F)
  write.csv(callset.q.20, file=paste0(outroot,".quartets.DNMs.20x.csv"),
            row.names=F, quote=F)
  write.table(dnm.20, file=paste0(outroot, ".triosquartets.20x.dnm"),
              col.names=F, row.names=F, quote=F, sep="\t")


}



if (interactive() == F) {
  Main()
}
