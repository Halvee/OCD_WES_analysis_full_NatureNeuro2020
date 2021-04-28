
library(denovolyzeR)
library(epitools)

## GLOBAL VARS
EFF.SETS <- list("Synon"=c("syn"),
                 "Mis"=c("misB","misP","misD","misU"),
                 "MisND"=c("misB","misP","misU"),
                 "MisD"=c("misD"),
                 "LoF"=c("non","splice","frameshift"),
                 "Mis/LoF"=c("non","splice","frameshift","misD",
                              "misP","misB","misU"))
EFF.SETS.TRANS <- list("syn"="Synon",
                       "mis"="Mis",
                       "misB"="MisB",
                       "misP"="MisP",
                       "misU"="MisU",
                       "misND"="MisND",
                       "misD"="MisD",
                       "lof"="LoF",
                       "misD_hc"="hcMisD")
EFF.SETS.CACO <- EFF.SETS

N.CO.TRIOS <- 1911

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T) 
  if (length(ARGS) != 13) {
    cat("run_rate_caco_analysis.R <gnomad.tsv> <exac.tsv> <ccds.geneset> <trios.sampleped>",
        "<cov.tsv> <samplecov.tsv> <genecov.tsv>",
        "<cases.cdnm> <ctrls.cdnm> <gene.mu.tsv>",
        "<cvg.genic.stats.tsv> <hcmisd.varlist>",
        "<outroot>\n")
    q()
  }

  # get ARGS
  gnomad.tsv <- ARGS[1]
  exac.tsv <- ARGS[2]
  ccds.geneset <- ARGS[3]
  trios.sampleped <- ARGS[4]
  cov.tsv <- ARGS[5]
  samplecov.tsv <- ARGS[6]
  genecov.tsv <- ARGS[7]
  cases.cdnm <- ARGS[8]
  ctrls.cdnm <- ARGS[9]
  gene.mu.tsv <- ARGS[10]
  cvg.genic.stats.tsv <- ARGS[11]
  mpcge2.varlist <- ARGS[12]
  outroot <- ARGS[13]

  # read input files
  gnomad<-read.table(gnomad.tsv,stringsAsFactors=F,
                     header=T,sep="\t")
  exac <- read.table(gnomad.tsv,header=T,sep="\t")
  ccds <- scan(ccds.geneset, what=character(),quiet=T)
  ped <- read.table(trios.sampleped, stringsAsFactors=F)
  colnames(ped) <- c("FID","IID","PID","MID","SEX","PHE","SEQTYPE","EXOMEPREPKIT")
  cov <- read.table(cov.tsv, sep="\t",
                    header=T, stringsAsFactors=F)
  samplecov <- read.table(samplecov.tsv, sep="\t",
                          header=T,stringsAsFactors=F)
  genecov <- read.table(genecov.tsv, sep="\t",
                        header=T, stringsAsFactors=F)
  ca.cdnm <- read.table(cases.cdnm, stringsAsFactors=F)
  co.cdnm <- read.table(ctrls.cdnm, stringsAsFactors=F)
  cdnm.cols <- c("IID","CHR","POS","VARID","REF","ALT",
                 "GENE","EFF","INKITS")
  colnames(ca.cdnm) <- cdnm.cols
  colnames(co.cdnm) <- cdnm.cols
  gene.mu <- read.table(gene.mu.tsv, header=T, sep="\t", stringsAsFactors=F)
  cvg.genic.stats <- read.table(cvg.genic.stats.tsv,header=T,sep="\t",stringsAsFactors=F,check.names=F)
  mpcge2.varids <- scan(mpcge2.varlist, what=character())

  # switch out gene.mu columns
  for (i in 1:ncol(gene.mu)) {
    if (colnames(gene.mu)[i] %in% names(EFF.SETS.TRANS)) {
      colnames(gene.mu)[i] <- EFF.SETS.TRANS[[ colnames(gene.mu)[i] ]]
    }
  }

  # collect info on ped trios
  fid.counts <- table(ped$FID)
  trio.fids <- names(fid.counts[fid.counts==3])
  n.ca.trios <- length(trio.fids)

  # collect fids for roche probands
  roche.iids <- ped[(ped$EXOMEPREPKIT == "Roche") & (ped$PHE==2), "IID"]
  roche.fids <- ped[(ped$EXOMEPREPKIT == "Roche") & (ped$PHE==2), "FID"]
  trio.fids.roche <- intersect(trio.fids, roche.fids)
  n.ca.trios.roche <- length(trio.fids.roche)
  #ca.cdnm <- subset(ca.cdnm, IID %in% roche.iids)
  #n.ca.trios <- n.ca.trios.roche

  ## ANALYSIS PART 1 : coverage summary stats
  samplecov <- read.table(samplecov.tsv, sep="\t",
                          header=T,stringsAsFactors=F)
  genecov <- read.table(genecov.tsv, sep="\t",
                        header=T, stringsAsFactors=F)
  cat("Average percent of CCDS bases jointly covered across trios :",
      mean(samplecov$Percent.Covered..Y.exc.), "\n")
  cat("Median percent of CCDS bases jointly covered across trios :",
      median(samplecov$Percent.Covered..Y.exc.), "\n")
  cat("Number of samples with > 90% of CCDS bases jointly covered :",
      nrow(subset(samplecov, Percent.Covered..Y.exc. > 0.9)), "\n")
  cat("Number of samples with <= 90% of CCDS bases jointly covered :",
      nrow(subset(samplecov, Percent.Covered..Y.exc. <= 0.9)), "\n")
  cat("Min percent of CCDS bases jointly covered in a trio :",
      min(samplecov$Percent.Covered..Y.exc.), "\n")
  cat("Percent of samples with > 90% of CCDS bases jointly covered :",
      nrow(subset(samplecov, Percent.Covered..Y.exc. > 0.9)) / nrow(samplecov),
      "\n")
  cat("Total number of genes included :",
      nrow(genecov), "\n")
  cat("Average percent of CCDS bases jointly covered across genes :",
      mean(genecov$X.Mean), "\n")
  cat("Median percent of CCDS bases jointly covered across genes :",
      median(genecov$X.Mean), "\n")
  cat("Number of genes with mean percent of CCDS bases cov > 90% :",
      nrow(subset(genecov, X.Mean > 90)), "\n")
  cat("Percent of genes with mean percent of CCDS bases cov > 90% :",
      nrow(subset(genecov, X.Mean > 90)) / nrow(genecov), "\n")
  cat("Number of genes with mean percent of CCDS bases cov <= 90% :",
      nrow(subset(genecov, X.Mean <= 90)), "\n")
  cat("Percent of genes with mean percent of CCDS bases cov <= 99% :",
      nrow(subset(genecov, X.Mean <= 99)) / nrow(genecov), "\n") 
  cat("Number of genes with mean percent of CCDS bases cov <= 99% :",
      nrow(subset(genecov, X.Mean <= 99)), "\n")

  ## ANALYSIS PART 2 : SNV and indel rate across trio samples
  fids_iids <- strsplit(samplecov[["Trio.SampleID"]],"_")
  iids <- c()
  for (i in 1:length(fids_iids)) { iids <- c(iids, fids_iids[[i]][2]) } 
  samplecov$IID <- iids
  rownames(samplecov) <- samplecov$IID
  samplecov$n_snv <- rep(0, nrow(samplecov))
  samplecov$n_indel <- rep(0, nrow(samplecov))
  for (i in 1:nrow(ca.cdnm)) {
    iid.i <- ca.cdnm[i, "IID"]
    if ((nchar(ca.cdnm[i,"REF"])==1) &  (nchar(ca.cdnm[i,"ALT"])==1)) {
      samplecov[iid.i, "n_snv"] <- samplecov[iid.i, "n_snv"] + 1
    } else {
      samplecov[iid.i, "n_indel"] <- samplecov[iid.i, "n_indel"] + 1
    }
  }
  samplecov$rate_snv <- samplecov$n_snv/(2*samplecov$Covered.Bases..Y.exc.)
  samplecov$rate_indel <- samplecov$n_indel/(2*samplecov$Covered.Bases..Y.exc.)
  cat("Mean number of SNVs per offspring :",
      mean(samplecov$n_snv), "\n")
  cat("Mean rate of indels per offspring :",
      mean(samplecov$n_indel), "\n")
  cat("Mean rate of SNVs per jointly covered nucleotide :",
      mean(samplecov$rate_snv), "\n")
  cat("Mean rate of indels per jointly covered nucleotide :",
      mean(samplecov$rate_indel), "\n")

  ## ANALYSIS PART 3 : CONSTRUCTION OF RATE, CACO TESTS

  # dfs for storing all rate.internal, rate.denovolyzeR and caco results
  rate.internal.all <- NULL
  rate.denovolyzeR.all <- NULL
  caco.all <- NULL

  # load rate tables
  denovolyzer.ptbl <- viewProbabilityTable(format="long")
 
  # get subsets of cdnm sets for case/control analysis
  ca.cdnm.caco.raw <- ca.cdnm
  co.cdnm.caco.raw <- co.cdnm
  ca.cdnm.caco <- subset(ca.cdnm.caco.raw, INKITS==TRUE)
  co.cdnm.caco <- subset(co.cdnm.caco.raw, INKITS==TRUE)

  # subset caco files on ccds
  ca.cdnm <- subset(ca.cdnm, GENE %in% ccds)
  ca.cdnm.caco <- subset(ca.cdnm.caco, GENE %in% ccds)
  co.cdnm.caco <- subset(co.cdnm.caco, GENE %in% ccds)

  # get subset of case cdnm set for rate analysis (internal and denovolyzeR)
  ca.cdnm.rate <- ca.cdnm
  ca.cdnm.rate.denovolyzer <- ca.cdnm.rate   
  for (mis.str in c("misB","misP","misD","misU")) {
    ca.cdnm.rate.denovolyzer$EFF <- gsub(mis.str,"mis", 
                                         ca.cdnm.rate.denovolyzer$EFF)
  }
  ca.cdnm.rate.denovolyzer <- subset(ca.cdnm.rate.denovolyzer, 
                                     EFF %in% c("syn","mis","non",
                                                "splice","frameshift"))
  ca.cdnm.rate.denovolyzer <- subset(ca.cdnm.rate.denovolyzer, 
                                     GENE %in% denovolyzer.ptbl$hgncSymbol)
 
  # split pli into genesets
  exome <- genecov$HGNC.CCDS
  gnomad <- subset(gnomad, gene %in% exome)
  gnomad.x <- subset(gnomad, is.na(oe_mis_pphen)==F)
  # intol <- gnomad.x[(gnomad.x$oe_lof_upper_bin <= 4) &
  #                  (gnomad.x$oe_mis_pphen <= 0.767230), "gene"]
  intol <- gnomad[(gnomad$oe_lof_upper_bin <= 4) &
                  (gnomad$oe_mis_upper <= 0.995), "gene"]
  #intol <- gnomad[gnomad$oe_lof_upper_bin <= 4, "gene"]
  #intol <- gnomad[(gnomad$oe_lof_upper < 0.35) | (gnomad$oe_mis_upper<0.70), "gene"]
  cat("N INTOLERANT GENES:", length(intol), "\n")
  # tol <- gnomad.x[(gnomad.x$oe_lof_upper_bin > 4) &
  #                 (gnomad.x$oe_mis_pphen > 0.767230), "gene"]
  tol <- gnomad[(gnomad$oe_lof_upper_bin > 4) |
                (gnomad$oe_mis_upper > 0.995), "gene"]
  #tol <- gnomad[gnomad$oe_lof_upper_bin > 4, "gene"]
  #tol <- gnomad[(gnomad$oe_lof_upper >= 0.35) & (gnomad$oe_mis_upper>=0.70), "gene"]
  cat("N TOLERANT GENES:", length(tol), "\n")
  ultraintol <- gnomad[gnomad$oe_lof_upper_bin == 0, "gene"]
  geneset.names <- c("tolerant","intolerant")
  genesets <- list("exome"=exome,"tolerant"=tol, "intolerant"=intol)

  # store LOEUF bins
  loeuf.deciles <- list()
  gnomad.x <- gnomad
  for (i in 0:9) {
    genes.i <- gnomad.x[gnomad.x$oe_lof_upper_bin == i, "gene"]
    loeuf.deciles[[paste0("LOEUF_",as.character(i+1))]] <- genes.i
  }

  # store mis_pphen bins
  gnomad.x <- subset(gnomad, is.na(oe_mis_pphen) == F)
  gnomad.x <- subset(gnomad.x, (oe_mis_pphen <= 0.767230) &
                               (oe_lof_upper_bin <= 4))
  mis.pphen.deciles <- list()
  for (i in 0:4) {
    genes.i <- gnomad.x[gnomad.x$oe_lof_upper_bin == i, "gene"]
    mis.pphen.deciles[[paste0("LOEUF_",as.character(i+1))]] <- genes.i 
  }

  # add misDhc annotations
  #misd.hc.i <- which((ca.cdnm$VARID %in% misd.hc.vars)
  #                   & 
  #                   (ca.cdnm$EFF == "misD"))
  #ca.cdnm[misd.hc.i,"EFF"] <- "misD_hc"
  #misd.hc.i <- which((ca.cdnm.caco$VARID %in% misd.hc.vars)
  #                   &
  #                   (ca.cdnm.caco$EFF == "misD"))
  #ca.cdnm.caco[misd.hc.i,"EFF"] <- "misD_hc"
  #misd.hc.i <- which((co.cdnm.caco$VARID %in% misd.hc.vars)
  #                  &
  #                  (co.cdnm.caco$EFF == "misD"))
  #co.cdnm.caco[misd.hc.i,"EFF"] <- "misD_hc"
  
  # form combined expected rate for LoF+hcMisD
  gene.mu[["MisD/LoF"]] <- gene.mu[["MisD"]] + gene.mu[["LoF"]]
  gene.mu[["Mis/LoF"]] <- gene.mu[["Mis"]] + gene.mu[["LoF"]]

  #### ANALYSIS PART 5 : exome / RVIS>=50% / RVIS < 50%

  ## ANALYSIS PART 5A : INTERNAL RATE TESTS 
  rate.df <-RateTestsInternal(ca.cdnm, gene.mu, n.ca.trios,
                              EFF.SETS, genesets)
  rate.internal.all <- rate.df
  write.csv(rate.df, file=paste0(outroot,".exome_rate.internal.csv"),
            row.names=F, quote=F)

  ## ANALYSIS PART 5B : DENVOLYZER RATE TESTS
  denovolyzer.res <- RateTestsDenovolyzeR(ca.cdnm, n.ca.trios, exome, genesets)
  rate.denovolyzeR.all <- denovolyzer.res
  write.csv(denovolyzer.res, 
            file=paste0(outroot,".exome_rate.denovolyzeR.csv"),
            row.names=F, quote=F)

  ## ANALYSIS PART 5C : CASE/CONTROL
  caco.df <- CacoTests(ca.cdnm.caco, co.cdnm.caco,
                       n.ca.trios, N.CO.TRIOS, 
                       genesets, EFF.SETS.CACO)
  caco.all <- caco.df
  write.csv(caco.df, file=paste0(outroot, ".exome_caco.csv"),
            row.names=F, quote=F)
  
  # do rate-based, case/control tests of LoF per LOEUF decile
  EFF.SETS <- list(
                   "LoF"=c("non","splice","frameshift"))

  rate.df.x <- RateTestsInternal(ca.cdnm, gene.mu, n.ca.trios,
                                 EFF.SETS, loeuf.deciles)
  caco.df.x <- CacoTests(ca.cdnm.caco, co.cdnm.caco,
                         n.ca.trios, N.CO.TRIOS,
                         loeuf.deciles, EFF.SETS)
  denovolyzer.res.x <- RateTestsDenovolyzeR(ca.cdnm, n.ca.trios, 
                                            exome, loeuf.deciles)
  rate.internal.all <- rbind(rate.internal.all, rate.df.x)
  rate.denovolyzeR.all <- rbind(rate.denovolyzeR.all, denovolyzer.res.x)
  caco.all <- rbind(caco.all, caco.df.x)
  write.csv(rate.df.x,
            file=paste0(outroot,".lof_LOEUF_deciles_rate.internal.csv"),
            row.names=F, quote=F)
  write.csv(denovolyzer.res.x,
            file=paste0(outroot,".lof_LOEUF_deciles_rate.denovolyzeR.csv"),
            row.names=F, quote=F)
  write.csv(caco.df.x,
            file=paste0(outroot,".lof_LOEUF_deciles_caco.csv"),
            row.names=F, quote=F)

  # do rate-based, case/control tests of MisD per 'MisD' decile
  EFF.SETS <- list("MisD"=c("misD"))
  rate.df.x <- RateTestsInternal(ca.cdnm, gene.mu, n.ca.trios,
                                 EFF.SETS,  mis.pphen.deciles)
  caco.df.x <- CacoTests(ca.cdnm.caco, co.cdnm.caco,
                         n.ca.trios, N.CO.TRIOS,
                         mis.pphen.deciles, EFF.SETS)
  write.csv(rate.df.x,
            file=paste0(outroot,".MisD_LOEUF_deciles_rate.internal.csv"),
            row.names=F, quote=F)
  write.csv(caco.df.x,
            file=paste0(outroot,".MisD_LOEUF_deciles_caco.csv"),
            row.names=F, quote=F)

  # do rate-based, case/control tests of hcMisD 
  EFF.SETS.HCMISD <- list(
                          "hcMisD"=c("misD_hc")
                          )
  misd.hc.i <- which((ca.cdnm$VARID %in% mpcge2.varids)
                     & 
                     (grepl("misD",ca.cdnm$EFF)))
  ca.cdnm[misd.hc.i,"EFF"] <- "misD_hc"
  misd.hc.i <- which((ca.cdnm.caco.raw$VARID %in% mpcge2.varids)
                     &
                     (grepl("misD",ca.cdnm.caco.raw$EFF)))
  ca.cdnm.caco.raw[misd.hc.i,"EFF"] <- "misD_hc"
  misd.hc.i <- which((co.cdnm.caco.raw$VARID %in% mpcge2.varids)
                    &
                     (grepl("misD",co.cdnm.caco.raw$EFF)))
                    #&
                    #(co.cdnm.caco$EFF == "misD"))
  co.cdnm.caco.raw[misd.hc.i,"EFF"] <- "misD_hc"
  rate.df.x <- RateTestsInternal(ca.cdnm, gene.mu, n.ca.trios,
                                 EFF.SETS.HCMISD, genesets)
  caco.df.x <- CacoTests(ca.cdnm.caco.raw, co.cdnm.caco.raw,
                         n.ca.trios, N.CO.TRIOS,
                         genesets, EFF.SETS.HCMISD)
  rate.internal.all <- rbind(rate.internal.all, rate.df.x)
  caco.all <- rbind(caco.all, caco.df.x)
  write.csv(rate.df.x,
            file=paste0(outroot,".MisD_hc_rate.internal.csv"),
            row.names=F, quote=F)
  write.csv(caco.df.x,
            file=paste0(outroot,".MisD_hc_caco.csv"),
            row.names=F, quote=F)
 
  # take genes with at least one LoF. Is there an excess of missense mutation
  # within these genes?
  lof.genes <- ca.cdnm[(ca.cdnm$EFF %in% c("non","splice","frameshift"))
                       &
                       (ca.cdnm$GENE %in% intol),
                       "GENE"]
  genesets <- list("LoF_Genes"=lof.genes)
  EFF.SETS <- list("Synon"=c("syn"),
                   "MisND"=c("misB","misP"),
                   "MisD"=c("misD","misD_hc")
                  )
  rate.df.x <- RateTestsInternal(ca.cdnm, gene.mu, n.ca.trios,
                                 EFF.SETS, genesets)
  caco.df.x <- CacoTests(ca.cdnm.caco, co.cdnm.caco,
                         n.ca.trios, N.CO.TRIOS,
                         genesets, EFF.SETS.CACO)
  rate.internal.all <- rbind(rate.internal.all, rate.df.x)
  caco.all <- rbind(caco.all, caco.df.x)
  write.csv(rate.df.x,
            file=paste0(outroot,".lof_genes.internal.csv"),
            row.names=F, quote=F)
  write.csv(caco.df.x,
            file=paste0(outroot,".lof_genes.caco.csv"),
            row.names=F, quote=F)

  # function for FET odds ratio table
  BuildTbl <- function(df, classif.str, res=NULL) {
    fet <- fisher.test(df)
    out.df <- data.frame(Annotation=classif.str,
                         ca_carrier=df[1,1],
                         ca_noncarrier=df[2,1],
                         ca_perc_carrier=df[1,1]/sum(df[,1]),
                         co_carrier=df[1,2],
                         co_noncarrier=df[2,2],
                         co_perc_carrier=df[1,2]/sum(df[,2]),
                         caco_perc_carrier_diff=df[1,1]/sum(df[,1])-df[1,2]/sum(df[,2]),
                         or_95ci_l=fet$conf.int[1],
                         or_95ci_u=fet$conf.int[2],
                         or_est=fet$estimate,
                         fet_p=fet$p.value)
    if (is.null(res)) { res <- out.df }
    else { res <- rbind(res, out.df) }
    return(res)
  }
  BuildTblVsSyn <- function(df, classif.str, res=NULL) {
    fet <- fisher.test(df, alternative='greater')
    out.df <- data.frame(Annotation=classif.str,
                         n_ca_var=df[1,1],
                         n_ca_syn=df[2,1],
                         n_co_var=df[1,2],
                         n_co_syn=df[2,2],
                         or_95ci_l=fet$conf.int[1],
                         or_95ci_u=fet$conf.int[2],
                         or_est=fet$estimate,
                         fet_p=fet$p.value)
    if (is.null(res)) { res <- out.df }
    else { res <- rbind(res, out.df) }
    return(res)
  }


  # do caco comparison of samples with an LoF in LOEUF<50% & oe_mis_pphen<0.76723
  gnomad.x <- subset(gnomad, is.na(oe_mis_pphen)==F)
  intol <- gnomad.x[(gnomad.x$oe_lof_upper_bin<=4) & 
                    (gnomad.x$oe_mis_upper<=0.995), "gene"]
  MISD.EFFS <- c("misD_hc","misD","misP","misB","misU")
  ca.syn <- ca.cdnm.caco[(ca.cdnm.caco$EFF == "syn") &
                         (ca.cdnm.caco$GENE %in% intol), "IID"]
  ca.misd <- ca.cdnm.caco[(ca.cdnm.caco$EFF %in% MISD.EFFS) &
                          (ca.cdnm.caco$GENE %in% intol), "IID"]
  ca.lof <- ca.cdnm.caco[(ca.cdnm.caco$EFF %in% c("non","splice","frameshift")) &
                         (ca.cdnm.caco$GENE %in% intol), "IID"]
  co.syn <- co.cdnm.caco[(co.cdnm.caco$EFF == "syn") &
                         (co.cdnm.caco$GENE %in% intol), "IID"]
  co.misd <- co.cdnm.caco[(co.cdnm.caco$EFF %in% MISD.EFFS) &
                          (co.cdnm.caco$GENE %in% intol), "IID"]
  co.lof <- co.cdnm.caco[(co.cdnm.caco$EFF %in% c("non","splice","frameshift")) &
                         (co.cdnm.caco$GENE %in% intol), "IID"]
  n.ca.syn <- length(unique(ca.syn))
  n.ca.misd <- length(unique(ca.misd))
  n.ca.lof <- length(unique(ca.lof))
  n.ca.lofmisd <- length(unique(c(ca.misd, ca.lof)))
  n.co.syn <- length(unique(co.syn))
  n.co.misd <- length(unique(co.misd))
  n.co.lof <- length(unique(co.lof))
  n.co.lofmisd <- length(unique(c(co.misd, co.lof)))
  syn.df <- data.frame(ca=c(n.ca.syn, n.ca.trios - n.ca.syn),
                       co=c(n.co.syn, N.CO.TRIOS - n.co.syn))
  misd.df <- data.frame(ca=c(n.ca.misd, n.ca.trios - n.ca.misd),
                        co=c(n.co.misd, N.CO.TRIOS - n.co.misd))
  lof.df <- data.frame(ca=c(n.ca.lof, n.ca.trios - n.ca.lof),
                       co=c(n.co.lof, N.CO.TRIOS - n.co.lof))
  lofmisd.df <- data.frame(ca=c(n.ca.lofmisd, n.ca.trios - n.ca.lofmisd),
                           co=c(n.co.lofmisd, N.CO.TRIOS - n.co.lofmisd))
  res.full <- NULL
  res.full <- BuildTbl(syn.df, "Synon")
  res.full <- BuildTbl(lof.df, "LoF, intol genes", res=res.full)
  res.full <- BuildTbl(misd.df, "Mis, intol genes", res=res.full) 
  res.full <- BuildTbl(lofmisd.df, "LoF/Mis, intol genes", res=res.full)
  write.csv(res.full, file=paste0(outroot, ".Mis_LoF_intol.fet.csv"),
            row.names=F, quote=F)

  # do caco comparison of samples with an LoF in LOEUF<10% or 
  # a MisD variant with MPC >= 2
  gnomad.x <- subset(gnomad, is.na(oe_mis_pphen)==F)
  lofintol <- gnomad[(gnomad$oe_lof_upper_bin==0), "gene"]
  print(table(ca.cdnm.caco$EFF))
  print(ca.cdnm[ca.cdnm$EFF == "misD_hc", ])
  mis.annots <- c("misU","misB","misP","misD","misD_hc")
  ca.syn <- ca.cdnm.caco.raw[(ca.cdnm.caco.raw$EFF == "syn"),
                             "IID"]
  ca.mis <- ca.cdnm.caco.raw[(ca.cdnm.caco.raw$EFF %in% mis.annots),
                             "IID"]
  ca.lof <- ca.cdnm.caco.raw[(ca.cdnm.caco.raw$EFF %in%
                              c("non","splice","frameshift")), "IID"]
  ca.misdhc <- ca.cdnm.caco.raw[(ca.cdnm.caco.raw$EFF %in% c("misD_hc")),
                                "IID"]
  ca.lofhc <- ca.cdnm.caco.raw[(ca.cdnm.caco.raw$EFF %in% c("non","splice","frameshift")) &
                               (ca.cdnm.caco.raw$GENE %in% lofintol), 
                               "IID"]
  co.syn <- co.cdnm.caco.raw[(co.cdnm.caco.raw$EFF == "syn"),
                             "IID"]
  co.mis <- co.cdnm.caco.raw[(co.cdnm.caco.raw$EFF %in% mis.annots),
                             "IID"]
  co.lof <- co.cdnm.caco.raw[co.cdnm.caco.raw$EFF %in% 
                             c("non","splice","frameshift"), "IID"]
  co.misdhc <- co.cdnm.caco.raw[(co.cdnm.caco.raw$EFF %in% c("misD_hc")),
                                "IID"]
  co.lofhc <- co.cdnm.caco.raw[(co.cdnm.caco.raw$EFF %in% c("non","splice","frameshift")) &
                               (co.cdnm.caco.raw$GENE %in% lofintol), 
                               "IID"]
  n.ca.syn <- length(unique(ca.syn))
  n.ca.misdhc <- length(unique(ca.misdhc))
  n.ca.misnd <- length(ca.mis) - n.ca.misdhc
  n.ca.lof <- length(unique(ca.lof))
  n.ca.lofhc <- length(unique(ca.lofhc))
  n.ca.synmisnd <- n.ca.syn + n.ca.misnd
  n.ca.lofmisdhc <- length(unique(c(ca.misdhc, ca.lofhc)))
  n.co.syn <- length(unique(co.syn))
  n.co.misdhc <- length(unique(co.misdhc))
  n.co.misnd <- length(co.mis) - n.co.misdhc 
  n.co.lof <- length(unique(co.lof))
  n.co.lofhc <- length(unique(co.lofhc))
  n.co.synmisnd <- n.co.syn + n.co.misnd
  n.co.lofmisdhc <- length(unique(c(co.misdhc, co.lofhc)))
  syn.df <- data.frame(ca=c(n.ca.syn, n.ca.trios - n.ca.syn),
                       co=c(n.co.syn, N.CO.TRIOS - n.co.syn))
  misdhc.df <- data.frame(ca=c(n.ca.misdhc, n.ca.trios - n.ca.misdhc),
                          co=c(n.co.misdhc, N.CO.TRIOS - n.co.misdhc))
  lofhc.df <- data.frame(ca=c(n.ca.lofhc, n.ca.trios - n.ca.lofhc),
                         co=c(n.co.lofhc, N.CO.TRIOS - n.co.lofhc))
  lofmisdhc.df <- data.frame(ca=c(n.ca.lofmisdhc, n.ca.trios - n.ca.lofmisdhc),
                             co=c(n.co.lofmisdhc, N.CO.TRIOS - n.co.lofmisdhc))
  res.full <- NULL
  res.full <- BuildTbl(syn.df, "Synon")
  res.full <- BuildTbl(lofhc.df, "hcLoF (LOEUF < 10%)", res=res.full)
  res.full <- BuildTbl(misdhc.df, "hcMisD (MPC>=2)", res=res.full) 
  res.full <- BuildTbl(lofmisdhc.df, "hcLoF/hcMisD", res=res.full)
  write.csv(res.full, file=paste0(outroot, ".hc_damaging.fet.csv"),
            row.names=F, quote=F)

  # do caco comparison of samples with an LoF in LOEUF<10% or 
  # a MisD variant with MPC >= 2, this time relative to synon var calls per
  # cohort
  nvar.ca.syn <- length(ca.syn)
  nvar.co.syn <- length(co.syn)
  nvar.ca.misdhc <- length(ca.misdhc)
  nvar.co.misdhc <- length(co.misdhc)
  nvar.ca.mis <- length(ca.mis)
  nvar.co.mis <- length(co.mis)
  nvar.ca.misnd <- length(ca.mis) - nvar.ca.misdhc
  nvar.co.misnd <- length(co.mis) - nvar.co.misdhc
  nvar.ca.synmisnd <- nvar.ca.syn + nvar.ca.misnd
  nvar.co.synmisnd <- nvar.co.syn + nvar.co.misnd
  nvar.ca.lof <- length(ca.lof)
  nvar.co.lof <- length(co.lof)
  nvar.ca.lofhc <- length(ca.lofhc) 
  nvar.co.lofhc <- length(co.lofhc) 
  nvar.ca.ns <- nvar.ca.mis + nvar.ca.lof
  nvar.co.ns <- nvar.co.mis + nvar.co.lof
  nvar.ca.lofmisdhc <- nvar.ca.misdhc + nvar.ca.lofhc
  nvar.co.lofmisdhc <- nvar.co.misdhc + nvar.co.lofhc
  mis.df <- data.frame(ca=c(nvar.ca.mis, nvar.ca.syn),
                       co=c(nvar.co.mis, nvar.co.syn))
  lof.df <- data.frame(ca=c(nvar.ca.lof, nvar.ca.syn),
                       co=c(nvar.co.lof, nvar.co.syn))
  ns.df <- data.frame(ca=c(nvar.ca.ns, nvar.ca.syn),
                      co=c(nvar.co.ns, nvar.co.syn))
  misdhc.df <- data.frame(ca=c(nvar.ca.misdhc, nvar.ca.syn),
                          co=c(nvar.co.misdhc, nvar.co.syn))
  lofhc.df <- data.frame(ca=c(nvar.ca.lofhc, nvar.ca.syn),
                         co=c(nvar.co.lofhc, nvar.co.syn))
  lofmisdhc.df <- data.frame(ca=c(nvar.ca.lofmisdhc, nvar.ca.syn),
                             co=c(nvar.co.lofmisdhc, nvar.co.syn))
  res.full <- NULL
  res.full <- BuildTblVsSyn(mis.df, "Missense", res=res.full)
  res.full <- BuildTblVsSyn(lof.df, "LoF", res=res.full)
  res.full <- BuildTblVsSyn(ns.df, "Nonsynon", res=res.full)
  res.full <- BuildTblVsSyn(lofhc.df, "hcLoF (LOEUF < 10%)", res=res.full)
  res.full <- BuildTblVsSyn(misdhc.df, "hcMisD (MPC>=2)", res=res.full) 
  res.full <- BuildTblVsSyn(lofmisdhc.df, "hcLoF/hcMisD", res=res.full)
  write.csv(res.full, file=paste0(outroot, ".hc_damaging.fet.rel_to_n_synon.csv"),
            row.names=F, quote=F)

  # write full set of results to tables
  write.csv(rate.internal.all,
            file=paste0(outroot,".tests.rate.internal.csv"),
            row.names=F, quote=F)
  write.csv(rate.denovolyzeR.all,
            file=paste0(outroot,".tests.rate.denovolyzeR.csv"),
            row.names=F, quote=F)
  write.csv(caco.all,
            file=paste0(outroot,".tests.caco.csv"),
            row.names=F, quote=F)


  q()

}

RatePerSample <- function(cdnm, genes, effs, n.trios, ...) {
  cdnm.x <- subset(cdnm, (GENE %in% genes) & (EFF %in% effs))
  n.dnm <- nrow(cdnm.x)
  iids <- cdnm.x$IID
  iid.counts <- table(iids)
  iids.0 <- n.trios - length(iid.counts)
  counts <- rep(0, iids.0)
  for (iid in names(iid.counts)) {
    counts <- c(counts, iid.counts[[iid]])
  }
  dat <- t.test(counts, ...)
  pois.dat <- pois.exact(n.dnm, n.trios)
  res <- list("count"=nrow(cdnm.x),
              "counts"=counts,
              "rate"=pois.dat$rate,
              "rate_95ci_l"=dat$conf.int[1],
              "rate_95ci_u"=dat$conf.int[2],
              "pois_95ci_l"=pois.dat$lower,
              "pois_95ci_u"=pois.dat$upper)
              #"pois_95ci_l"=cdnm.x - sqrt(n.dnm/n.trios),
              #"pois_95ci_u"=cdnm.x + sqrt(n.dnm/n.trios))
  return(res)
}

RatePerNt <- function(cdnm, geneset, effs, cov.df) {
  cov.df.i <- subset(cov.df, HGNC.CCDS %in% geneset)
  cov.df.i <- cov.df.i[, -(1:3)]
  rates.df <- data.frame(FID_IID=colnames(cov.df.i))
  iids <- c()
  ncovs <-c()
  for (fid_iid in colnames(cov.df.i)) {
    n.nt <- sum(cov.df.i[[fid_iid]])
    iid <- strsplit(fid_iid, "_")[[1]][2]
    iids <- c(iids, iid)
    ncovs <- c(ncovs, n.nt)
  }
  rates.df$IID <- iids
  rates.df$NCOV <- ncovs
  rownames(rates.df) <- rates.df$IID
  rates.df$COUNT <- rep(0,nrow(rates.df))

  cdnm.i <- subset(cdnm, (GENE %in% geneset) & (EFF %in% effs))
  counts <- table(cdnm.i$IID)     
  for (iid in names(counts)) {
    rates.df[iid, "COUNT"] <- counts[[iid]]
  }
  rates.df$RATE <- rates.df$COUNT / rates.df$NCOV
  stats <- t.test(rates.df$RATE)
  total.rate <- sum(rates.df$RATE)
  return(total.rate)
}

RateTestsInternal <- function(ca.cdnm, gene.mu, n.ca.trios,
                              eff.sets, genesets.list) {
  rate.df <- data.frame(geneset=character(),
                        eff=character(),
                        count_obs=numeric(),
                        rate_obs=numeric(),
                        count_exp=numeric(),
                        rate_exp=numeric(),
                        rate_ratio=numeric(),
                        pois_95ci_l=numeric(),
                        pois_95ci_u=numeric(),
                        rate_95ci_l=numeric(),
                        rate_95ci_u=numeric(),
                        ptest_p=numeric())
  for (geneset.name in names(genesets.list)) {
    ca.cdnm.i <- subset(ca.cdnm, GENE %in% genesets.list[[geneset.name]])
    gene.mu.i <- subset(gene.mu, Gene %in% genesets.list[[geneset.name]])
    for (eff.set.name in names(eff.sets)) {
      count.exp <- sum(gene.mu.i[[eff.set.name]]) 
      rate.exp <- count.exp / n.ca.trios 
      obs <- RatePerSample(ca.cdnm, genesets.list[[geneset.name]],
                           eff.sets[[eff.set.name]], n.ca.trios)
      ptest.res <- poisson.test(obs$count, T=n.ca.trios, 
                                r=rate.exp, alternative="two.sided")$p.value
      rate.df<-rbind(rate.df,
                     data.frame(geneset=geneset.name,
                                eff=eff.set.name,
                                count_exp=round(count.exp,1),
                                count_obs=round(obs$count,1),
                                rate_exp=signif(rate.exp,3),
                                rate_obs=signif(obs$rate,3),
                                rate_ratio=signif(obs$rate/rate.exp,3), 
                                pois_95ci_l=signif(obs$pois_95ci_l,3),
                                pois_95ci_u=signif(obs$pois_95ci_u,3),
                                rate_95ci_l=signif(obs$rate_95ci_l,3),
                                rate_95ci_u=signif(obs$rate_95ci_u,3),
                                ptest_p=signif(ptest.res,4))
                     )
    }
  }
  rownames(rate.df) <- NULL
  return(rate.df)
}

RateTestsDenovolyzeR <- function(ca.cdnm, n.trios,
                                 exome.genes, genesets.list) {

  # replace all missense annotations to missense
  ca.cdnm$EFF <- gsub("mis[A-Z]*","mis",ca.cdnm$EFF)
  denovolyzer.ptbl <- viewProbabilityTable(format="long")
  denovolyzer.genes <- intersect(exome.genes, denovolyzer.ptbl$hgncSymbol)
  denovolyzer.ptbl <- subset(denovolyzer.ptbl, 
                             hgncSymbol %in% denovolyzer.genes)
  denovolyzer.res <- NULL
  dnm.genes.effs <- ca.cdnm[,c("GENE","EFF")]
  dnm.genes.effs <- subset(dnm.genes.effs, GENE %in% denovolyzer.genes)
  dnm.genes.effs <- subset(dnm.genes.effs, EFF %in% c("syn","mis","non",
                                                      "splice","frameshift"))
  for (geneset.name in names(genesets.list)) {
    geneset <- genesets.list[[geneset.name]]
    dnm.genes.effs.i <- subset(dnm.genes.effs, GENE %in% geneset)
    denovolyzer.ptbl.i <- subset(denovolyzer.ptbl, hgncSymbol %in% geneset)
    res <- denovolyze(as.character(dnm.genes.effs.i$GENE),
                      as.character(dnm.genes.effs.i$EFF),
                      nsamples=n.trios,
                      includeClasses = c("syn", "mis","lof","prot"),
                      probTable=denovolyzer.ptbl.i)
    res <- cbind(data.frame("geneset"=geneset.name), res)
    if (is.null(denovolyzer.res)) { denovolyzer.res <- res }
    else { denovolyzer.res <- rbind(denovolyzer.res, res) }
  }
  return(denovolyzer.res)
}

CacoTests <- function(ca.cdnm, co.cdnm, 
                      n.ca.trios, n.co.trios, 
                      genesets.list, eff.sets,
                      by.carriers=F) {
  caco.df <- data.frame(geneset=character(),
                        eff=character(),
                        ca_rate=numeric(),
                        co_rate=numeric(),
                        ca_rate_95ci_l=numeric(),
                        ca_rate_95ci_u=numeric(),
                        co_rate_95ci_l=numeric(),
                        co_rate_95ci_u=numeric(),
                        ca_pois_95ci_l=numeric(), 
                        ca_pois_95ci_u=numeric(), 
                        n_ca_dnm1=numeric(),
                        n_ca_dnm0=numeric(),
                        n_co_dnm1=numeric(),
                        n_co_dnm0=numeric(),
                        caco_or_est=numeric(),
                        caco_or_95ci_l=numeric(),
                        caco_or_95ci_u=numeric(),
                        caco_p=numeric(),
                        pois_p=numeric())
  if (by.carriers == T) {
    n.ca.trios <- length(unique(ca.cdnm[,1]))
    n.co.trios <- length(unique(co.cdnm[,1]))
  }
  for (geneset.name in names(genesets.list)) {
    ca.cdnm.i <- subset(ca.cdnm, GENE %in% genesets.list[[geneset.name]])
    co.cdnm.i <- subset(co.cdnm, GENE %in% genesets.list[[geneset.name]])
    for (eff.set.name in names(eff.sets)) {
      ca.obs <- RatePerSample(ca.cdnm, genesets.list[[geneset.name]],
                              eff.sets[[eff.set.name]], n.ca.trios)
      co.obs <- RatePerSample(co.cdnm, genesets.list[[geneset.name]],
                              eff.sets[[eff.set.name]], n.co.trios)
      ca.cdnm.ij <- subset(ca.cdnm.i, EFF %in% eff.sets[[eff.set.name]])
      co.cdnm.ij <- subset(co.cdnm.i, EFF %in% eff.sets[[eff.set.name]])
      ca.ij.dnm1 <- unique(ca.cdnm.ij$IID)
      co.ij.dnm1 <- unique(co.cdnm.ij$IID)
      n.ca.dnm1 <- length(ca.ij.dnm1)
      n.co.dnm1 <- length(co.ij.dnm1)
      n.ca.dnm0 <- n.ca.trios - n.ca.dnm1
      n.co.dnm0 <- n.co.trios - n.co.dnm1
      fet.dnm.tbl <- data.frame(ca=c(n.ca.dnm1,n.ca.dnm0),
                                co=c(n.co.dnm1,n.co.dnm0))
      fet.dnm.res <- fisher.test(fet.dnm.tbl,alternative="two.sided")
      ptest.res <- poisson.test(nrow(ca.cdnm.ij),
                                T=n.ca.trios, r=co.obs$rate,
                                alternative='two.sided')$p.value
      caco.df <- rbind(caco.df, 
                       data.frame(geneset=geneset.name,
                                  eff=eff.set.name,
                                  ca_rate=signif(ca.obs$rate,3),
                                  co_rate=signif(co.obs$rate,3),
                                  ca_rate_95ci_l=signif(ca.obs$rate_95ci_l,3),
                                  ca_rate_95ci_u=signif(ca.obs$rate_95ci_u,3),
                                  co_rate_95ci_l=signif(co.obs$rate_95ci_l,3),
                                  co_rate_95ci_u=signif(co.obs$rate_95ci_u,3),
                                  ca_pois_95ci_l=signif(ca.obs$pois_95ci_l,3),
                                  ca_pois_95ci_u=signif(ca.obs$pois_95ci_u,3),
                                  n_ca_dnm1=n.ca.dnm1,
                                  n_ca_dnm0=n.ca.dnm0,
                                  n_co_dnm1=n.co.dnm1,
                                  n_co_dnm0=n.co.dnm0,
                                  caco_or_est=round(fet.dnm.res$estimate,3),
                                  caco_or_95ci_l=round(fet.dnm.res$conf.int[1],3),
                                  caco_or_95ci_u=round(fet.dnm.res$conf.int[2],3),
                                  caco_p=signif(fet.dnm.res$p.value,4),
                                  pois_p=signif(ptest.res,3))      
                      )
    }
  }
  return(caco.df)
}

if (interactive() == F) {
  Main()
}

