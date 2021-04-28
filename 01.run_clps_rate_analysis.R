#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(gridExtra)

## DATA
IN_DRAGENDB_CA_CSV<-"data/clps_rate/OCD_clps_rate.cases.dragendb.csv"
IN_DRAGENDB_CO_CSV<-"data/clps_rate/OCD_clps_rate.ctrls.dragendb.csv"
CCDS_GENESET<-"results/input/CCDSr20.geneset"
GNOMAD_LOEUF_TSV<-"results/input/gnomad.v2.1.1.loeuf_bins.tsv"
GNOMAD_PLI_TSV<-"results/input/gnomad.v2.1.1.pLI.tsv"
NONPSYCH_PLI_TSV<-"data/input/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"

## PARAM
LOF.COL <- "firebrick4"

Main <- function() {

  # get user args defining input cohort
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 5) {
    cat("03.run_clps_rate_analysis.R <in.sampleped> ",
        "<in.evec> <in.syn.genotypes.csv>",
        "<in.lof.genotypes.csv> <outroot>\n")
    q()
  }
  SAMPLEPED<-ARGS[1]
  EVEC<-ARGS[2]
  LVG_SYN_CSV<-ARGS[3]
  LVG_LOF_CSV<-ARGS[4]
  OUTROOT<-ARGS[5]

  # read sampleped, prepare columns for analysis
  ped <- read.table(SAMPLEPED, stringsAsFactors=F)
  colnames(ped) <- c("FID","IID","PID","MID","SEX","PHE","SEQTYPE","PREPKIT")
  samples <- ped$IID
  cases <- ped[ped$PHE == 2, "IID"]
  ctrls <- ped[ped$PHE == 1, "IID"]
  ped$SEX <- ped$SEX - 1
  ped$PHE <- ped$PHE - 1

  # read dragendb csvs, merge with ped
  ddb.ca <- read.csv(IN_DRAGENDB_CA_CSV, stringsAsFactors=F)
  ddb.co <- read.csv(IN_DRAGENDB_CO_CSV, stringsAsFactors=F)
  ddb <- rbind(ddb.ca, ddb.co)
  ddb <- ddb[, c("CHGVID","CCDSCoverage","CCDSBasesCov10X",
                 "African_prob","Caucasian_prob","EastAsian_prob",
                 "Hispanic_prob", "MiddleEastern_prob", "SouthAsian_prob")]
  colnames(ddb)[1] <- "IID"
  ped <- merge(ped, ddb, by="IID")

  # read evec
  evec <-read.table(EVEC, stringsAsFactors=F)
  n.pc <- ncol(evec) - 2
  colnames(evec) <- c("IID",paste0("PC",1:n.pc),"PHENO")
  evec$PHENO <- NULL
  ped <- merge(ped, evec, by="IID")

  # assign IID as rownames of ped
  rownames(ped) <- ped$IID

  # read genesets form, intolerant genes
  ccds <- scan(CCDS_GENESET, quiet=T, what=character())
  gnomad <- read.table(GNOMAD_LOEUF_TSV, sep="\t",
                       header=T, stringsAsFactors=F)
  gnomad <- subset(gnomad, gene %in% ccds)
  gnomad <- subset(gnomad, is.na(oe_lof_upper_bin)==F)

  # store gene symbols for each individual LOEUF decile to list
  lof.deciles <- list()
  for (i in 0:9) {
    gnomad.d <- subset(gnomad, oe_lof_upper_bin == i)
    d.genes <- gnomad.d$gene
    j <- i + 1
    lof.deciles[[paste0("LOEUFd",j)]] <- d.genes
  }

  # read genotypes, apply general formatting to dataframe
  df.syn <- fread(LVG_SYN_CSV, stringsAsFactors=F, data.table=F)
  df.lof <- fread(LVG_LOF_CSV, stringsAsFactors=F, data.table=F)
  
  # let's make sure that only synon variants are in synon set,
  # and that only LoF variants are in LoF set 
  df.syn <- subset(df.syn, grepl("synonymous_variant", Effect))
  lof.annots <-  c("frameshift_variant",
                   "splice_donor_variant",
                   "splice_acceptor_variant",
                   "stop_gained")
  df.lof <- subset(df.lof, Effect %in% lof.annots)
 
  # combine syn and lof dataframes
  df <- rbind(df.syn, df.lof)

  # format the merged dataframe
  colnames(df) <- gsub(" ",".", colnames(df))
  df$Gene.Name <- gsub("'","",df$Gene.Name)
  
  # get variant counts (synon, lof)
  count.df <- data.frame(step=character(), 
                         n_syn=numeric(),
                         n_lof=numeric())
  n_syn <- nrow(subset(df, grepl("synonymous_variant", Effect)))
  lof.annots <-  c("frameshift_variant",
                   "splice_donor_variant",
                   "splice_acceptor_variant",
                   "stop_gained")
  n_lof <- nrow(subset(df, Effect %in% lof.annots))
  count.df <-rbind(count.df,
                   data.frame(step="input",
                              n_syn=n_syn,n_lof=n_lof)
                  )

  ## requirements for LOF rate tests :

  # 1. heterozygous variant calls in autosomes
  df <- subset(df, 
               (GT == "het")
               &
               (grepl("X-", Variant.ID)==F)
              )
  n_syn <- nrow(subset(df, grepl("synonymous_variant", Effect)))
  n_lof <- nrow(subset(df, Effect %in% lof.annots))
  count.df <-rbind(count.df,
                   data.frame(step="het_auto_only",
                              n_syn=n_syn,n_lof=n_lof)
                  )  

  # 2. variants are of high quality (PASS only in case/control,
  #    no QC-fail variants in cases or controls, no non-PASS in gnomAD)
  df <- subset(df, FILTER == "PASS")
  df <- subset(df, QC.Fail.Case == 0)
  df <- subset(df, QC.Fail.Ctrl == 0)
  for (g.filter in c("gnomAD.Exome.FILTER","gnomAD.Genome.FILTER")) {
    df[[g.filter]] <- ifelse(is.na(df[[g.filter]]), "PASS", df[[g.filter]])
    df <- df[df[[g.filter]] == "PASS", , drop=F]
  }
  n_syn <- nrow(subset(df, grepl("synonymous_variant", Effect)))
  n_lof <- nrow(subset(df, Effect %in% lof.annots))
  count.df <-rbind(count.df,
                   data.frame(step="quality_extensive_pruning",
                              n_syn=n_syn,n_lof=n_lof)
                  )

  # 3. variants are near-absent (MAF < 0.00001) from full gnomAD 
  #    and ExAC exome cohorts. Since we're focused on LoF in genes
  #    that may be more intolerant, in many of those cases we
  #    wouldn't expect these variants to be found in more than 1-3 
  #    individuals in a cohort of 60K-150K
  df <- subset(df, LOO.AF < 0.00001)
  for (af.col in c("gnomAD.Exome.global_AF","ExAC.global.af")) {
    df[[af.col]] <- ifelse(is.na(df[[af.col]]), 0, df[[af.col]])
    df <- df[df[[af.col]] < 0.00001, , drop=F]
  }
  n_syn <- nrow(subset(df, grepl("synonymous_variant", Effect)))
  n_lof <- nrow(subset(df, Effect %in% lof.annots))
  count.df <-rbind(count.df,
                   data.frame(step="MAF_lt_1e-5",
                              n_syn=n_syn,n_lof=n_lof)
                  )

  # write pruning count table to df
  write.csv(count.df,
            file=paste0(OUTROOT,".var_pruning.counts.csv"),
            row.names=F, quote=F)


  ## split by variant effect :

  # synonymous
  df.syn <- subset(df, grepl("synonymous_variant", Effect))

  # loss of function
  lof.annots <-  c("frameshift_variant",
                   "splice_donor_variant",
                   "splice_acceptor_variant",
                   "stop_gained")
  df.lof <- subset(df, Effect %in% lof.annots) 

  # establish list of df sets
  df.sets <- list(
                  "syn"=df.syn, 
                  "LoF"=df.lof
                 )

  # load genesets
  genesets <- list("exome"=ccds)
  for (d.name in names(lof.deciles)) {
    genesets[[d.name]] <- lof.deciles[[d.name]]
  }

  # for each geneset, get number of variants per annot for each sample
  for (geneset.name in names(genesets)) {
    for (annot.name in names(df.sets)) {
      g.a <- paste0(geneset.name, "_", annot.name)
      ped <- CountsFromDf(df.sets[[annot.name]], 
                          genesets[[geneset.name]], 
                          g.a, 
                          ped) 
    }
  }

  # plot density of variants across all samples
  pdf(paste0(OUTROOT,".exome_syn.density.all_samples.pdf"))
  plot(density(ped$exome_syn), main="all samples", 
       xlab="synonymous variants",
       ylab="density")
  dev.off()
  pdf(paste0(OUTROOT,".exome_syn.density.cases_vs_ctrls.pdf"))
  plot(density(ped$exome_syn), main="cases vs controls", 
       xlab="synonymous variants",
       ylab="density")
  lines(density(subset(ped, PHE==0)$exome_syn), col='blue')
  lines(density(subset(ped, PHE==1)$exome_syn), col='red')
  dev.off()

  # what is the mean synonymous variant count in cases and in controls?
  cat("Mean exome synon rate (cases) :",
      mean(subset(ped,PHE==1)$exome_syn),"\n")
  cat("Mean exome synon rate (controls) :",
      mean(subset(ped,PHE==0)$exome_syn),"\n")

  # in general, how much of a predictor for OCD status is total exome count?
  print(summary(glm(PHE~exome_syn, family=binomial, data=ped)))
  print(summary(lm(exome_syn~PHE, data=ped)))

  # load variable for sample exome kit represented in both cases and controls:
  # Roche v3, IDTERPv1
  ped$SHAREDKITS <- ifelse(ped$PREPKIT %in% c("Roche","IDTERPv1"), 1, 0)

  ## get counts per annotation, ID which covars are predictors
  ## of synon var count
  covars <- c(
              "CCDSBasesCov10X",
              "SEX",
              "SHAREDKITS",
              paste0("PC",1:n.pc)
             )

  # if only ocd samples in ped, then we are doing a supplemental analysis
  # on LoF rate in cases versus controls requested by reviewers, and 
  # need to only use a subset of the covariates above
  n.non_ocd <- nrow(subset(ped, grepl("ocd",IID)==F))
  if (n.non_ocd == 0) {
    
    # case-only supplemental analysis. remove SEX and SHAREDKITS
    covars <- c(
                "CCDSBasesCov10X",
                paste0("PC",1:n.pc)
               )
  
  }
  
  covars.keep <- c()
  for (covar in covars) {

    # linear 
    mdl.str <- paste0("exome_syn~",covar)
    mdl <- lm(as.formula(mdl.str), data=ped)
    res <- summary(mdl)
    estimate <- res$coefficients[covar,1]
    pval.linear <- res$coefficients[covar,4]
    cat("LINEAR :",covar,":",estimate," (p=",pval.linear,")\n")

    # logistic
    mdl.str <- paste0("PHE~",covar)
    mdl <- glm(as.formula(mdl.str),family=binomial,data=ped)
    res <- summary(mdl)
    estimate <- res$coefficients[covar,1]
    pval.logistic <- res$coefficients[covar,4]
    cat("LOGISTIC :",covar,":",estimate," (p=",pval.logistic,")\n")
  
    # if nom signif in linear or logistic models add as a covar
    if ((pval.linear < 0.01) | (pval.logistic < 0.01)) { 
      covars.keep <- c(covars.keep, covar)
    }
  }
  cat("Selecting covariates to keep best on pval_linear < 0.01 .. \n")
  cat("Covariates to keep : ",paste(covars.keep, collapse=", "),"\n")

  # test mdl with predictive covars
  covars.str <- paste(covars.keep, collapse="+")
  mdl.str<-paste0("PHE~",covars.str,"+exome_syn")
  mdl <- glm(as.formula(mdl.str),family=binomial,data=ped)  
  print(summary(mdl))
  mdl.str <- paste0("exome_syn~",covars.str,"+PHE")
  mdl <- lm(as.formula(mdl.str), data=ped)
  print(summary(mdl))

  # test model on exome-wide LoF rate
  covars.str <- paste(covars.keep, collapse="+")
  mdl.str<-paste0("PHE~",covars.str,"+exome_syn+exome_LoF")
  mdl <- glm(as.formula(mdl.str),family=binomial,data=ped)  
  print(summary(mdl))
  mdl.str <- paste0("exome_LoF~",covars.str,"+exome_syn+PHE")
  mdl <- lm(as.formula(mdl.str), data=ped)
  print(summary(mdl))

  ## setup and execute linear regression
  annot.names <- c("LoF")
  res.df <- NULL
  geneset.names <- names(genesets)
  for (annot in annot.names) {
    for (geneset in names(lof.deciles)) {
      cat(geneset,"/",annot,"\n")
      res.df <- LogisticReg(ped, covars.str, geneset, annot, res.df=res.df)    
    }
  }

  ## setup and execute complementary FET analyses of LoF vs syn var
  annot.names <- c("LoF")
  fet.df <- NULL
  geneset.names <- names(genesets)
  for (annot in annot.names) {
   for (geneset in names(lof.deciles)) {
      cat(geneset,"/",annot,"\n")
      fet.df <- FetRelTo(ped, 
                         geneset, 
                         annot.1=annot,
                         annot.0="syn",
                         res.df=fet.df)
    }
  }

  # store full unformatted dfs for writing to file later
  res.full.df <- res.df
  fet.full.df <- fet.df 

  # get subset of rows with LoF in LOUEF bins
  res.lof <- subset(res.df, grepl("_LoF", res.df$predictor))
  res.lof <- subset(res.lof, grepl("LOEUFd", res.lof$predictor))
  res.lof$predictor <- gsub("LOEUFd","",res.lof$predictor)
  res.lof$predictor <- gsub("_LoF","",res.lof$predictor)
  res.lof$predictor <- factor(res.lof$predictor, levels=res.lof$predictor)  

  # establish pd
  pd <- position_dodge(width = 0.4)

  # plot LoF in LOEUF deciles (logistic regression)
  gg.rate.loeuf <- ggplot(res.lof, 
                          aes(y=odds_ratio_logistic, x=predictor))
  gg.rate.loeuf <- gg.rate.loeuf + geom_linerange(aes(ymin=ci95_low_logistic,
                                                      ymax=ci95_high_logistic),
                                      position=pd, color=LOF.COL )
  gg.rate.loeuf <- gg.rate.loeuf + geom_point(shape = 19, size  = 4, 
                                              color=LOF.COL,
                                              position=pd, shape=19)
  
  gg.rate.loeuf <- gg.rate.loeuf + ylab("Case/control odds ratio")
  gg.rate.loeuf <- gg.rate.loeuf + xlab("LOEUF decile")
  gg.rate.loeuf <- gg.rate.loeuf + geom_hline(yintercept=1)
  gg.rate.loeuf <- gg.rate.loeuf + theme(legend.position = 'none')
  gg.rate.loeuf <- gg.rate.loeuf + theme(axis.title.x = element_text(face = "bold",
                                                                     size = 14),
                                         axis.title.y = element_text(face = "bold",
                                                                     size = 14),
                                         axis.text.x = element_text(face = "bold",
                                                                    size = 14, angle = 0),
                                         axis.text.y = element_text(face = "bold",
                                                                    size = 14, angle = 45)
                                        )
  ggsave(paste0(OUTROOT, ".LOF_LOEUF.logistic.pdf"))

  # plot LoF in LOEUF deciles (linear regression)
  gg.rate.loeuf <- ggplot(res.lof,
                          aes(y=estimate_linear, x=predictor))
  gg.rate.loeuf <- gg.rate.loeuf + geom_linerange(aes(ymin=ci95_low_linear,
                                                      ymax=ci95_high_linear),
                                      position=pd, color=LOF.COL)
  gg.rate.loeuf <- gg.rate.loeuf + geom_point(shape = 19, size  = 4, 
                                              color=LOF.COL,
                                              position=pd, shape=19)
  gg.rate.loeuf <- gg.rate.loeuf + theme(axis.text.x = element_text(face = "bold", 
                                                                    size = 14, angle = 0),
                                         axis.text.y = element_text(face = "bold",
                                                                    size = 14, angle = 45)
                                         )
  gg.rate.loeuf <- gg.rate.loeuf + ylab("Case rate diff. versus controls")
  gg.rate.loeuf <- gg.rate.loeuf + xlab("LOEUF decile")
  gg.rate.loeuf <- gg.rate.loeuf + geom_hline(yintercept=0)
  gg.rate.loeuf <- gg.rate.loeuf + theme(legend.position = 'none')
  gg.rate.loeuf <- gg.rate.loeuf + theme(axis.title.x = element_text(face = "bold",
                                                                     size = 14),
                                         axis.title.y = element_text(face = "bold",
                                                                     size = 14),
                                         axis.text.x = element_text(face = "bold",
                                                                    size = 14, angle = 45),
                                         axis.text.y = element_text(face = "bold",
                                                                    size = 14, angle = 45)
                                        )
  ggsave(paste0(OUTROOT, ".LOF_LOEUF.linear.pdf"))

  # format fet input
  print(fet.df)
  res.lof <- subset(fet.df, grepl("LOEUFd", geneset))
  print(res.df)
  res.lof$geneset <- gsub("LOEUFd","",res.lof$geneset)
  res.lof$geneset <- factor(res.lof$geneset, levels=res.lof$geneset)  

  print(res.lof)

  # plot fet results across LOEUF deciles
  gg.fet.loeuf <- ggplot(res.lof,
                          aes(y=fet_or, x=geneset))
  gg.fet.loeuf <- gg.fet.loeuf + geom_linerange(aes(ymin=fet_or_95ci_l,
                                                    ymax=fet_or_95ci_u),
                                                position=pd, color=LOF.COL)
  gg.fet.loeuf <- gg.fet.loeuf + geom_point(shape = 19, size  = 4, 
                                            color=LOF.COL,
                                            position=pd, shape=19)
  gg.fet.loeuf <- gg.fet.loeuf + theme(axis.text.x = element_text(face = "bold", 
                                                                   size = 14, angle = 0),
                                        axis.text.y = element_text(face = "bold",
                                                                   size = 14, angle = 45)
                                       )
  gg.fet.loeuf <- gg.fet.loeuf + ylab("LoF/synon case/control odds ratio")
  gg.fet.loeuf <- gg.fet.loeuf + xlab("LOEUF decile")
  gg.fet.loeuf <- gg.fet.loeuf + geom_hline(yintercept=1)
  gg.fet.loeuf <- gg.fet.loeuf + theme(legend.position = 'none')
  gg.fet.loeuf <- gg.fet.loeuf + theme(axis.title.x = element_text(face = "bold",
                                                                     size = 14),
                                         axis.title.y = element_text(face = "bold",
                                                                     size = 14),
                                         axis.text.x = element_text(face = "bold",
                                                                    size = 14, angle = 45),
                                         axis.text.y = element_text(face = "bold",
                                                                    size = 14, angle = 45)
                                        )
  ggsave(paste0(OUTROOT, ".LOF_LOEUF.fet.pdf"), plot=gg.fet.loeuf)

  # combine LOEUF linear and fet plots into one figure, for supp
  pdf(paste0(OUTROOT, ".LOF_LOEUF.linear_fet.pdf"), width=6, height=7)
  grid.arrange(gg.rate.loeuf, gg.fet.loeuf, ncol=1)
  dev.off()

  # addition : assess LoF burden in genes with pLI > 0.995 from gnomAD
  pli <- read.table(GNOMAD_PLI_TSV, header=T, stringsAsFactors=F, sep="\t")
  pli <- subset(pli, is.na(pLI)==F)
  pli_genes_0 <- subset(pli, pLI<=0.5)$gene
  pli_genes_1 <- subset(pli, (pLI > 0.5) & (pLI <= 0.995))$gene
  pli_genes_2 <- subset(pli, pLI > 0.995)$gene
  pli_genes.loeuf1 <- intersect(pli_genes_2, lof.deciles[["LOEUFd1"]])

  # assess LoF burden in genes with pLI > 0.995 from a nonpsychiatric
  # cohort of samples (ExAC nonpsych)
  pli <- read.table(NONPSYCH_PLI_TSV, header=T, stringsAsFactors=F,
                    sep="\t")
  pli <- subset(pli, is.na(pLI)==F)
  pli_nonpsych <- subset(pli, pLI > 0.995)$gene
  pli_nonpsych.loeuf1 <- intersect(pli_nonpsych, 
                                   lof.deciles[["LOEUFd1"]])

  # get counts
  genesets = list("pLIle0.5"=pli_genes_0,
                  "pLI0.5to0.995"=pli_genes_1,
                  "pLIgt0.995"=pli_genes_2,
                  "pLIgt0.995_LOEUF1"=pli_genes.loeuf1,
                  "nonpsych_pLIgt0.995_LOEUF1"=pli_nonpsych.loeuf1)
  for (geneset.name in names(genesets)) {
    for (annot.name in names(df.sets)) {
      g.a <- paste0(geneset.name, "_", annot.name)
      ped <- CountsFromDf(df.sets[[annot.name]], 
                          genesets[[geneset.name]], 
                          g.a, 
                          ped) 
    }
  }
  ## setup and execute linear regression
  annot.names <- c("LoF")
  res.df <- NULL
  geneset.names <- names(genesets)
  for (annot in annot.names) {
    for (geneset in geneset.names) {
      cat(geneset,"/",annot,"\n")
      res.df <- LogisticReg(ped, covars.str, 
                            geneset, annot, 
                            res.df=res.df)    
    }
  }

  ## setup and execute fisher's exact tests
  annot.names <- c("LoF")
  fet.df <- NULL
  geneset.names <- names(genesets)
  for (annot in annot.names) {
    for (geneset in geneset.names) {
      cat(geneset, "/", annot, "\n")
      fet.df <- FetRelTo(ped, 
                         geneset, 
                         annot.1=annot,
                         annot.0="syn",
                         res.df=fet.df)
    }
  }
 
  # write results to csv
  write.csv(res.df,
            file=paste0(OUTROOT, ".pLI.csv"),
            row.names=F,
            quote=F)
  write.csv(fet.df,
            file=paste0(OUTROOT,".pLI.fet.csv"),
            row.names=F,
            quote=F)

  # add pli results to other dfs from LOEUF
  res.full.df <- rbind(res.full.df, res.df)
  fet.full.df <- rbind(fet.full.df, fet.df)

  # format geneset names, p-values
  res.df$predictor <- c("pLI<=0.5",
                        "pLI>0.5,\n<=0.995", 
                        "pLI>0.995",
                        "pLI>0.995,\nLOEUF<10%",
                        "nonpsych\npLI>0.995,\nLOEUF<10%")
  res.df$predictor <- factor(res.df$predictor,
                             levels=res.df$predictor)
  res.df$p_value_logistic <- PvalFix(res.df$p_value_logistic)
  res.df$p_value_linear <- PvalFix(res.df$p_value_linear)

  # plot results
  gg.logistic <- ggplot(res.df, 
                        aes(y=odds_ratio_logistic, x=predictor, 
                            label=p_value_logistic))
  gg.logistic <- gg.logistic + geom_linerange(aes(ymin=ci95_low_logistic,
                                                      ymax=ci95_high_logistic),
                                              position=pd, color=LOF.COL )
  gg.logistic <- gg.logistic + geom_point(shape = 19, size  = 4, 
                                            color=LOF.COL,
                                            position=pd, shape=19)
  
  gg.logistic <- gg.logistic + ylab("Case/control odds ratio")
  gg.logistic <- gg.logistic + geom_hline(yintercept=1)
  gg.logistic <- gg.logistic + theme(legend.position = 'none')
  gg.logistic <- gg.logistic + theme(axis.title.x=element_blank(),
                                     axis.title.y = element_text(face = "bold",
                                                                 size = 11),
                                     axis.text.x = element_text(face = "bold",
                                                                size = 11, angle = 0),
                                     axis.text.y = element_text(face = "bold",
                                                                size = 11, angle = 45)
                                    )
  gg.logistic <- gg.logistic + geom_text(color="black", angle=90, size=4,                               
                                         position=position_nudge(x=0.2))
 
  # plot LoF in pLI bins (linear regression)
  gg.linear <- ggplot(res.df,
                      aes(y=estimate_linear, x=predictor, label=p_value_linear))
  gg.linear <- gg.linear + geom_linerange(aes(ymin=ci95_low_linear,
                                                      ymax=ci95_high_linear),
                                      position=pd, color=LOF.COL)
  gg.linear <- gg.linear + geom_point(shape = 19, size  = 4, 
                                              color=LOF.COL,
                                              position=pd, shape=19)
  gg.linear <- gg.linear + ylab("Case rate diff. versus controls")
  gg.linear <- gg.linear + geom_hline(yintercept=0)
  gg.linear <- gg.linear + theme(legend.position = 'none')
  gg.linear <- gg.linear + theme(axis.title.x=element_blank(),
                                         axis.title.y = element_text(face = "bold",
                                                                     size = 11),
                                         axis.text.x = element_text(face = "bold",
                                                                    size = 11,
                                                                    angle = 0),
                                         axis.text.y = element_text(face = "bold",
                                                                    size = 11, 
                                                                    angle = 45)
                                        )
  gg.linear <- gg.linear + geom_text(color="black", angle=90, size=4,
                                     position=position_nudge(x=0.2))

  # format fet input
  print(fet.df)
  res.lof <- subset(fet.df, grepl("pLI", geneset))
  print(res.lof)
  res.lof$geneset <- factor(res.lof$geneset, levels=res.lof$geneset)
  res.lof$geneset <- c("pLI<=0.5",
                       "pLI>0.5,\n<=0.995",
                       "pLI>0.995",
                       "pLI>0.995,\nLOEUF<10%",
                       "nonpsych\npLI>0.995,\nLOEUF<10%")
  res.lof$geneset <- factor(res.lof$geneset, levels=res.lof$geneset)
  res.lof$fet_p <- PvalFix(res.lof$fet_p)

  # plot LoF in pLI bins (FET)
  gg.fet <- ggplot(res.lof,
                   aes(y=fet_or, x=geneset, label=fet_p))
  gg.fet <- gg.fet + geom_linerange(aes(ymin=fet_or_95ci_l,
                                        ymax=fet_or_95ci_u),
                                      position=pd, color=LOF.COL)
  gg.fet <- gg.fet + geom_point(shape = 19, size  = 4, 
                                              color=LOF.COL,
                                              position=pd, shape=19)
  gg.fet <- gg.fet + ylab("LoF/synon odds ratio")
  gg.fet <- gg.fet + geom_hline(yintercept=1)
  gg.fet <- gg.fet + theme(legend.position = 'none')
  gg.fet <- gg.fet + theme(axis.title.x=element_blank(),
                                         axis.title.y = element_text(face = "bold",
                                                                     size = 11),
                                         axis.text.x = element_text(face = "bold",
                                                                    size = 11,
                                                                    angle = 0),
                                         axis.text.y = element_text(face = "bold",
                                                                    size = 11, 
                                                                    angle = 45)
                                        )
  gg.fet <- gg.fet + geom_text(color="black", angle=90, size=4,
                                     position=position_nudge(x=0.2))

  # make one figure with logistic + linear res for pLI
  pdf(paste0(OUTROOT, ".pLI.pdf"), width=6, height=7)
  grid.arrange(gg.logistic, gg.linear, gg.fet, ncol=1)
  dev.off()

  # save the pLI figures by themselves (logistic only, linear only)
  ggsave(paste0(OUTROOT,".logistic.pLI.pdf"),
         plot=gg.logistic, height=6, width=7)
  ggsave(paste0(OUTROOT,".linear.pLI.pdf"),
         plot=gg.linear, height=6, width=7)
  ggsave(paste0(OUTROOT,".fet.pLI.pdf"),
         plot=gg.fet, height=6, width=7)

  # LOEUF<10% + pLI>0.995 per LoF annotation? (stop-gain snv, splice snv, framshfit indel)
  
  # split into annots
  df.non <- subset(df, Effect == "stop_gained")
  df.splice <- subset(df, Effect %in% c("splice_donor_variant","splice_acceptor_variant"))
  df.fs <- subset(df, Effect == "frameshift_variant")
  df.sets <- list("non"=df.non, "splice"=df.splice, "fs"=df.fs)

  # counts per annot
  genesets = list(
                  "pLIgt0.995_LOEUF1"=pli_genes.loeuf1
                 )
  for (geneset.name in names(genesets)) {
    for (annot.name in names(df.sets)) {
      g.a <- paste0(geneset.name, "_", annot.name)
      ped <- CountsFromDf(df.sets[[annot.name]],
                          genesets[[geneset.name]],
                          g.a,
                          ped)
    }
  }

  ## setup and execute regression
  annot.names <- c("non","splice","fs")
  res.df <- NULL
  geneset.names <- names(genesets)
  for (annot in annot.names) {
    for (geneset in geneset.names) {
      cat(geneset,"/",annot,"\n")
      res.df <- LogisticReg(ped, covars.str, geneset, annot, res.df=res.df)
    }
  }

  # write results per LoF classification
  write.csv(res.df,
            file=paste0(OUTROOT, ".LOEUF1_pLI0.995.per_LoF_annot.csv"),
            row.names=F, quote=F)
  print(res.df)

  # write full complete results to files
  write.table(res.full.df, file=paste0(OUTROOT, ".rate_tests.res.tsv"),
              row.names=F, col.names=T, sep="\t", quote=F)
  write.table(fet.full.df, file=paste0(OUTROOT, ".fet_tests.res.tsv"),
              row.names=F, col.names=T, sep="\t", quote=F)
  write.table(ped, file=paste0(OUTROOT, ".counts.tsv"),
              row.names=F, col.names=T, sep="\t", quote=F)

}

InitCountMatrix <- function(df, ped, genes) {
  df.count <- df[, c("Gene.Name","Sample.Name")]
  genes.qual <- unique(sort(df.count$Gene.Name))
  iids.qual <- unique(sort(df.count$Sample.Name))
  mat <- as.data.frame.matrix(table(df.count))

  # fill in columns with samples that have no calls
  iids0 <- setdiff(ped[,2], iids.qual)
  for (iid in iids0) { mat[[iid]] <- rep(0, nrow(mat)) }

  return(mat)
}

CountsFromDf <- function(df, genes, classif, res.tbl) {
  df <- df[, c("Gene.Name","Sample.Name")]
  df <- df[df$Gene.Name %in% genes, , drop=F]
  res.tbl[[classif]] <- rep(0, nrow(res.tbl))
  counts <- table(df$Sample.Name)
  for (iid.i in names(counts)) { res.tbl[iid.i, classif] <- counts[[iid.i]] }
  return(res.tbl)
}

MakeDensityPlot<-function(ped, classif, out.pdf) {
  pdf(out.pdf)
  plot(density(ped[[classif]]), main=classif)
  lines(density(ped[ped$PHE==0, classif]), col='blue')
  lines(density(ped[ped$PHE==1, classif]), col='red')
  dev.off()
}

LogisticReg <- function(df, covars.str, geneset.name, 
                        annot, fet.relto=NULL, res.df=NULL) {
  # define pred, geneset_syn strings
  pred <- paste0(geneset.name, "_", annot)
  gs.syn <- paste0(geneset.name, "_syn")
  
  # make exome counts seperate from geneset counts
  if (geneset.name != "exome") {
    gs.syn <- paste0(geneset.name, "_syn")
    eqn.str.lg <- paste0("PHE~",covars.str,
                         "+exome_syn+",gs.syn,"+",pred)
    eqn.str.lm <- paste0(pred,"~",covars.str,
                         "+exome_syn+",gs.syn,"+PHE")
  } else {
    eqn.str.lg <- paste0("PHE~",covars.str,"+",pred)
    eqn.str.lm <- paste0(pred,"~",covars.str,"+PHE")
  }
  print(eqn.str.lg)
  print(eqn.str.lm)

  # linear regression
  mdl.lm <- lm(as.formula(eqn.str.lm), data=df)
  res.lm <- summary(mdl.lm)

  # logistic regression
  mdl.lg <- glm(as.formula(eqn.str.lg),data=df, family=binomial)
  res.lg <- summary(mdl.lg)
  print(res.lg$coefficients)

  # store model results
  res.df.i <- data.frame(predictor=pred,
                         control_rate=mean(subset(df,PHE==0)[[pred]]),
                         estimate_linear=res.lm$coefficients["PHE",1],
                         ci95_low_linear=confint(mdl.lm)["PHE",1],
                         ci95_high_linear=confint(mdl.lm)["PHE",2], 
                         p_value_linear=res.lm$coefficients["PHE",4],
                         odds_ratio_logistic=exp(coef(mdl.lg)[pred]),
                         ci95_low_logistic=exp(confint(mdl.lg)[pred,1]),
                         ci95_high_logistic=exp(confint(mdl.lg)[pred,2]),
                         p_value_logistic=res.lg$coefficients[pred,4]
                         )
  if (is.null(res.df)) {
    res.df <- res.df.i
  } else {
    res.df <- rbind(res.df, res.df.i)
  }
  return(res.df)
}

# do two-sided fishers exact test of LoF variation counts relative to
# counts of synonymous variation
FetRelTo <- function(df, geneset.name, 
                     annot.1="LoF",
                     annot.0="syn",
                     res.df=NULL) {
  ann.1.col <- paste0(geneset.name, "_", annot.1)
  ann.0.col <- paste0(geneset.name, "_", annot.0)
  n.1.ca <- sum(subset(df, PHE==1)[[ann.1.col]])
  n.0.ca <- sum(subset(df, PHE==1)[[ann.0.col]])
  n.1.co <- sum(subset(df, PHE==0)[[ann.1.col]])
  n.0.co <- sum(subset(df, PHE==0)[[ann.0.col]])
  fet.df <- data.frame(ca=c(n.1.ca, n.0.ca),
                       co=c(n.1.co, n.0.co))
  fet.res <- fisher.test(fet.df)
  res.df.i <- data.frame(geneset=geneset.name,
                         ann_1=annot.1,
                         ann_0=annot.0,
                         n_1_ca=n.1.ca,
                         n_0_ca=n.0.ca,
                         n_1_co=n.1.co,
                         n_0_co=n.0.co,
                         fet_or=fet.res$estimate,
                         fet_or_95ci_l=fet.res$conf.int[1],
                         fet_or_95ci_u=fet.res$conf.int[2],
                         fet_p=fet.res$p.value)
  if (is.null(res.df)) {
    res.df <- res.df.i
  } else {
    res.df <- rbind(res.df, res.df.i)
  }
  return(res.df)
}

# fix pval cols
PvalFix <- function(pvals) {
  pvals <- ifelse(pvals <= 1,
                  paste0("p = ",formatC(pvals,format='e',digits=2)),
                  "")
  return(pvals)
}

if (interactive() == F){
  Main()
}
