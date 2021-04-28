#!/usr/bin/env Rscript

library(denovolyzeR)
library(ggplot2)

## DATA
GNOMAD_TSV <- "results/input/gnomad.v2.1.1.loeuf_bins.tsv"
CCDS_GENESET <- "results/input/CCDSr20_autosomal.geneset"
TRIOS_SAMPLEPED <- "results/dnm_cohort/OCDfams_2019.cohort.trios.sampleped"
COV_TSV <- "results/trio_coverage/OCDfams_2019.trios.cvg_stats.10x.tsv"
SAMPLECOV_TSV <- "results/trio_coverage/OCDfams_2019.trios.cvg_stats.10x.sampleProfileSummary.tsv"
GENECOV_TSV <- "results/trio_coverage/OCDfams_2019.trios.cvg_stats.10x.genicProfileSummary.tsv"
CASES_CDNM <- "results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.inkits.cdnm"
CTRLS_CDNM <- "results/dnm_casectrl/Iossifov_etal_2014_healthysib_DNMs.ppn2hdiv.inkits.cdnm"
GENE_MU_TSV <- "results/trio_coverage/OCDfams_2019.trios.10x.exp_dnm_rates.tsv"
CVG_GENIC_STATS_TSV <- "results/trio_coverage/OCDfams_2019.trios.cvg_stats.10x.genicProfileSummary.tsv"
MPC0TO1_VARLIST <- "results/input/MPC0to1.varlist.gz"
MPC1TO2_VARLIST <- "results/input/MPC1to2.varlist.gz"
MPCGT2_VARLIST <- "results/input/MPCgt2.varlist.gz"
GNOMAD_NONNEURO_VARLIST <- "data/input/gnomAD_exome_non_neuro_popmaxMAFgt0.varlist"
NDD_GENESET <- "results/input/NDD.n187.geneset"
TS_GE2DNM_GENESET <- "results/input/TS.ge2_dmg_DNMs_q_lt30.geneset"
TS_1DNM_GENESET <- "results/input/TS.ge1_dmg_DNM_q_ge30.geneset"
OUTROOT <- "results/dnm_rate_caco/OCDfams_2019.cohort.trios.10x"

## PARAM
N.CO.TRIOS <- 1911
COLORBLIND_ORANGE=rgb(0.9,0.6,0)
COLORBLIND_BLUE=rgb(0.0,0.45,0.70)
COLORBLIND_GREEN=rgb(0.0,0.6,0.5)
LOF_COLOR="firebrick4"
DAMAGING_COLOR <- rgb(0.80,0.40,0,1)

Main <- function() {

  # read input files
  gnomad<-read.table(GNOMAD_TSV,stringsAsFactors=F,
                     header=T,sep="\t")
  ccds <- scan(CCDS_GENESET, what=character(),quiet=T)
  ped <- read.table(TRIOS_SAMPLEPED, stringsAsFactors=F)
  colnames(ped) <- c("FID","IID","PID","MID","SEX","PHE","SEQTYPE","EXOMEPREPKIT")
  cov <- read.table(COV_TSV, sep="\t",
                    header=T, stringsAsFactors=F)
  samplecov <- read.table(SAMPLECOV_TSV, sep="\t",
                          header=T,stringsAsFactors=F)
  genecov <- read.table(GENECOV_TSV, sep="\t",
                        header=T, stringsAsFactors=F)
  ca.cdnm <- read.table(CASES_CDNM, stringsAsFactors=F)
  co.cdnm <- read.table(CTRLS_CDNM, stringsAsFactors=F)
  cdnm.cols <- c("IID","CHR","POS","VARID","REF","ALT",
                 "GENE","EFF","INKITS")
  colnames(ca.cdnm) <- cdnm.cols
  colnames(co.cdnm) <- cdnm.cols
  gene.mu <- read.table(GENE_MU_TSV, header=T, sep="\t", stringsAsFactors=F)
  cvg.genic.stats <- read.table(CVG_GENIC_STATS_TSV,header=T,sep="\t",
                                stringsAsFactors=F,check.names=F)

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
  samplecov <- read.table(SAMPLECOV_TSV, sep="\t",
                          header=T,stringsAsFactors=F)
  genecov <- read.table(GENECOV_TSV, sep="\t",
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

  # load denovolyzeR rate table
  denovolyzer.ptbl <- viewProbabilityTable(format="long")
  denovolyzer.ptbl <- subset(denovolyzer.ptbl, hgncSymbol %in% ccds)
 
  # subset caco files on ccds
  ca.cdnm <- subset(ca.cdnm, GENE %in% ccds)
  co.cdnm <- subset(co.cdnm, GENE %in% ccds)

  # get subset of case cdnm set for rate analysis (internal and denovolyzeR)
  ca.cdnm.rate.denovolyzer <- ca.cdnm
  for (mis.str in c("misB","misP","misD","misU")) {
    ca.cdnm.rate.denovolyzer$EFF <- gsub(mis.str,"mis", 
                                         ca.cdnm.rate.denovolyzer$EFF)
  }
  ca.cdnm.rate.denovolyzer <- subset(ca.cdnm.rate.denovolyzer, 
                                     EFF %in% c("syn","mis","non",
                                                "splice","frameshift"))
  ca.cdnm.rate.denovolyzer <- subset(ca.cdnm.rate.denovolyzer, 
                                     GENE %in% denovolyzer.ptbl$hgncSymbol)
  res.dr <- denovolyze(ca.cdnm.rate.denovolyzer$GENE, ca.cdnm.rate.denovolyzer$EFF,
                       nsamples=587, includeClasses=c("syn","mis","lof"),
                       probTable=denovolyzer.ptbl)
  write.csv(res.dr, row.names=F, quote=F,
            file=paste0(OUTROOT,"exome.annotations.rate_denovolyzeR.csv"))

  # base callset : remove all non-autosomal calls from case and control callsets
  ca.cdnm <- subset(ca.cdnm, (CHR %in% c("X","Y"))==F)
  co.cdnm <- subset(co.cdnm, (CHR %in% c("X","Y"))==F)

  # analysis 1 : comparison of proportion of cases and controls that 
  # carry a qualifying de novo based on coding annotation. Do this using
  # a logistic regression model where phenotype is outcome, main predictor
  # is dnm count per annotation_x, and covariate is number of dnms in sample
  # that are outside of capture kit intersect.
  # annotation 1 : synonymous
  # annotation 2 : missense non-damaging (ppn2 hdiv < 0.957)
  # annotation 3 : missense damaging (ppn2 hdiv >= 0.957)
  # annotation 4 : LoF (nonsense, splice donor/acceptor, frameshift)
  # annotation 5 : missense damaging or LOF

  # establish callset to use
  ca.cdnm.caco <- CdnmMultiMutRemove(ca.cdnm, dist.thresh=200)
  co.cdnm.caco <- CdnmMultiMutRemove(co.cdnm, dist.thresh=200) 

  # get list of total dnm counts (non-zero) in cases and in controls
  n_dnm_ca<-table(ca.cdnm.caco$IID)
  n_dnm_co<-table(co.cdnm.caco$IID)

  # store phenotypes of sample IDs with nonzero amounts of dnms
  phenotypes <- list()
  for (iid_i in names(n_dnm_ca)) { phenotypes[[iid_i]] <- 1 }
  for (iid_i in names(n_dnm_co)) { phenotypes[[iid_i]] <- 0 }

  # form a single case/control data frame
  caco.cdnm <- rbind(ca.cdnm.caco, co.cdnm.caco)
  caco.iids <- unique(sort(caco.cdnm$IID))

  # init output data frame
  a1_count_df <- data.frame(iid=character(), phe=numeric(), n_dnm=numeric(),
                            n_outofkits_dnm=numeric(), n_inkits_dnm=numeric(),
                            n_snv_dnm=numeric(),
                            n_indel_dnm=numeric(), n_syn_dnm=numeric(),
                            n_misnd_dnm=numeric(), n_misd_dnm=numeric(),
                            n_fsindel_dnm=numeric(), n_lofsnv_dnm=numeric())

  intol_genes <- gnomad[gnomad$oe_lof_upper_bin==0, "gene"]  
  # for each case with a non-zero amount of dnms, store the number of dnms
  # per individual coding annotation
  for (iid_i in caco.iids) {
    phe_i <- phenotypes[[iid_i]]
    n_dnm_i <- nrow(subset(caco.cdnm, IID == iid_i))
    n_snv_dnm_i <- nrow(subset(caco.cdnm,
                               (IID == iid_i) & (nchar(REF)==1) & (nchar(ALT)==1)))
    n_indel_dnm_i <- nrow(subset(caco.cdnm,
                                 (IID == iid_i) & ((nchar(REF)==1) | (nchar(ALT)!=1))))
    n_outofkits_dnm_i <- nrow(subset(caco.cdnm, 
                                     IID==iid_i & INKITS == FALSE))
    n_inkits_dnm_i <- nrow(subset(caco.cdnm,
                                  IID==iid_i & INKITS == TRUE))
    n_syn_dnm_i <- nrow(subset(caco.cdnm, 
                              IID==iid_i & EFF=="syn"))
    n_misnd_dnm_i <- nrow(subset(caco.cdnm, 
                                IID==iid_i & EFF %in% c("misB","misP")))
    n_misd_dnm_i <- nrow(subset(caco.cdnm, 
                               IID==iid_i & EFF %in% c("misD")))
    n_fsindel_dnm_i <- nrow(subset(caco.cdnm, 
                                  IID==iid_i & EFF %in% c("frameshift")))
    n_lofsnv_dnm_i <- nrow(subset(caco.cdnm, 
                                 IID==iid_i & EFF %in% c("non","splice")))   

    # add values to count df
    a1_count_df <- rbind(a1_count_df,
                         data.frame(iid=iid_i, 
                                    phe=phe_i, 
                                    n_dnm=n_dnm_i,
                                    n_outofkits_dnm=n_outofkits_dnm_i, 
                                    n_inkits_dnm=n_inkits_dnm_i,
                                    n_snv_dnm=n_snv_dnm_i,
                                    n_indel_dnm=n_indel_dnm_i, 
                                    n_syn_dnm=n_syn_dnm_i,
                                    n_misnd_dnm=n_misnd_dnm_i, 
                                    n_misd_dnm=n_misd_dnm_i,
                                    n_fsindel_dnm=n_fsindel_dnm_i, 
                                    n_lofsnv_dnm=n_lofsnv_dnm_i)
                        )

  }

  # get set of case trio IIDs that don't have any DNM calls
  ca.iids <- subset(ped, PHE==2)$IID
  ca.iids.dnm0 <- ca.iids[(ca.iids %in% a1_count_df$iid)==F]  

  # make set of synthetic control trio IIDs for samples that 
  # don't have DNM calls based on the total number of healthy
  # control trios in the published cohort (1911)
  n.co.dnm1 <- nrow(subset(a1_count_df, phe==0))
  n.co.dnm0 <- N.CO.TRIOS - n.co.dnm1
  co.iids.dnm0 <- paste0("co_",(n.co.dnm1+1):N.CO.TRIOS)

  # get total number of samples with zero DNMs
  n.dnm0 <- length(ca.iids.dnm0) + length(co.iids.dnm0)

  # add all samples with zero DNMs to full count table
  a1_count_df <- rbind(a1_count_df,
                       data.frame(iid=c(ca.iids.dnm0, co.iids.dnm0), 
                                  phe=c(rep(1,length(ca.iids.dnm0)), 
                                        rep(0, length(co.iids.dnm0))), 
                                  n_dnm=rep(0,n.dnm0),                            
                                  n_outofkits_dnm=rep(0,n.dnm0), 
                                  n_inkits_dnm=rep(0,n.dnm0),
                                  n_snv_dnm=rep(0,n.dnm0),           
                                  n_indel_dnm=rep(0,n.dnm0), 
                                  n_syn_dnm=rep(0,n.dnm0),
                                  n_misnd_dnm=rep(0,n.dnm0), 
                                  n_misd_dnm=rep(0,n.dnm0),
                                  n_fsindel_dnm=rep(0,n.dnm0), 
                                  n_lofsnv_dnm=rep(0,n.dnm0)) 
                      )

  # get count sums and subsets of interest (misD/LOF, non_frameshift indels)
  a1_count_df$n_lof_dnm <- a1_count_df$n_lofsnv_dnm + a1_count_df$n_fsindel_dnm
  a1_count_df$n_misdlofsnv_dnm <- a1_count_df$n_misd_dnm + a1_count_df$n_lofsnv_dnm
  a1_count_df$n_misdlof_dnm <- a1_count_df$n_misd_dnm + a1_count_df$n_lof_dnm

  # remove samples with ndnm > 5
  # a1_count_df <- subset(a1_count_df, n_dnm <= 5)

  # do some basic logistic regression models ..
  # 1. does outside-kit DNM count predict case/control status?
  # 2. does in-kit DNM count predict case/control status?
  outofkits_mdl <- glm(phe~n_outofkits_dnm, data=a1_count_df, family=binomial)
  inkits_mdl <- glm(phe~n_inkits_dnm, data=a1_count_df, family=binomial)
  print("OUTOFKITS : ")
  print(summary(outofkits_mdl))
  print("INKITS : ")
  print(summary(inkits_mdl))
  print("OUTOFKITS/INKITS : ")
  inoutofkits_mdl <- glm(phe~n_outofkits_dnm+n_inkits_dnm,
                         data=a1_count_df, family=binomial)
  print(summary(inoutofkits_mdl))

  # form result dataframe for analysis 1
  a1_df <- data.frame(Annotation=character(),
                      mean_case=numeric(),
                      mean_ctrl=numeric(),
                      linear_estimate=numeric(),
                      linear_ci95_low=numeric(),
                      linear_ci95_high=numeric(),
                      linear_p_value=numeric(),
                      logistic_odds_ratio=numeric(),
                      logistic_ci95_low=numeric(),
                      logistic_ci95_high=numeric(),
                      logistic_p_value=numeric())

  # formal tests on each of the following annotations : synon, misND, misD, LOF, misD_or_LOF
  # for (ann_i in c("syn","misnd","misd","lof","misdlof")) {
  for (ann_i in c("syn","misnd","misd","lof","misdlof")) {                      
    pred_i <- paste0("n_",ann_i,"_dnm")
    logistic_mdl_str <- paste0("phe~n_outofkits_dnm+",pred_i)
    linear_mdl_str <- paste0(pred_i,"~n_outofkits_dnm+phe")
    logistic_mdl_i <- glm(as.formula(logistic_mdl_str), data=a1_count_df, family=binomial)
    linear_mdl_i <- lm(as.formula(linear_mdl_str), data=a1_count_df)
    logistic_res_i <- summary(logistic_mdl_i)
    linear_res_i <- summary(linear_mdl_i)
    oddsratios_i <- exp(coef(logistic_mdl_i))

    a1_df <- rbind(a1_df,
                   data.frame(Annotation=ann_i, 
                              mean_case=mean(subset(a1_count_df, phe==1)[[pred_i]]),
                              mean_ctrl=mean(subset(a1_count_df, phe==0)[[pred_i]]),
                              linear_estimate=linear_res_i$coefficients["phe",1],
                              linear_ci95_low=confint(linear_mdl_i)["phe",1],
                              linear_ci95_high=confint(linear_mdl_i)["phe",2],
                              linear_p_value=linear_res_i$coefficients["phe",4],
                              logistic_odds_ratio=exp(coef(logistic_mdl_i)[pred_i]),
                              logistic_ci95_low=exp(confint(logistic_mdl_i)[pred_i,1]),
                              logistic_ci95_high=exp(confint(logistic_mdl_i)[pred_i,2]),
                              logistic_p_value=logistic_res_i$coefficients[pred_i,4])
                   )
  }

  # write analysis 1 results to csv
  write.csv(a1_df, quote=F, row.names=F,
            file=paste0(OUTROOT,".exome.annotations.caco.csv"))

  # plot results of analysis 1 (exome, case/control, per annotation)
  a1_df <- subset(a1_df, Annotation %in% c("syn","misnd","misd","lof"))
  a1_df$Annotation <- c("synon","misND","misD","LoF")
  pd <- position_dodge(width = 0.9)
  a1_df$logistic_odds_ratio <- as.numeric(a1_df$logistic_odds_ratio)
  a1_df$Annotation <- factor(a1_df$Annotation, 
                              levels=a1_df$Annotation)
  a1_df$logistic_p_value <- PvalFix(a1_df$logistic_p_value)
  gg.a1 <- ggplot(a1_df, aes(x=Annotation, y=logistic_odds_ratio, 
                             label=logistic_p_value))
  gg.a1 <- gg.a1 + geom_hline(yintercept=1,color='black')
  gg.a1 <- gg.a1 + geom_linerange(aes(ymin=logistic_ci95_low,                   
                                      ymax=logistic_ci95_high),             
                                      position=pd,
                                      color=c(rep(COLORBLIND_BLUE, 2),
                                              COLORBLIND_GREEN, LOF_COLOR),
                                 )
  gg.a1 <- gg.a1 + geom_point(shape = 19, size  = 4, position=pd, shape=19,
                              color=c(rep(COLORBLIND_BLUE, 2),
                                      COLORBLIND_GREEN, LOF_COLOR),
                             )
  gg.a1 <- gg.a1 + theme(text=element_text(size=12))
  gg.a1 <- gg.a1 + theme(legend.position = 'none',
                         axis.title.x = element_text(face = "bold", size = 13),
                         axis.title.y = element_text(face = "bold", size = 13),
                         axis.text.x = element_text(face = "bold", size = 13,
                                                    angle = 0),
                         axis.text.y = element_text(hjust=1,
                                                    face = "bold", size = 13, 
                                                    angle = 45)
                        )
  gg.a1 <- gg.a1 + geom_text(color="black", angle=90, size=4,
                             position=position_nudge(x=0.2))
  gg.a1 <- gg.a1 + ylab("case/control odds ratio")
  gg.a1 <- gg.a1 + xlab("Annotation")
  ggsave(paste0(OUTROOT, ".exome.annotations.caco.pdf"))

  # make additional figure from analysis 1 that represents differences in
  # case rate per trio and control rate per trio
  pd <- position_dodge(width = 0.9)
  a1_df$linear_estimate <- as.numeric(a1_df$linear_estimate)
  a1_df$linear_p_value <- PvalFix(a1_df$linear_p_value)
  gg.a1.r <- ggplot(a1_df, aes(x=Annotation, y=linear_estimate, 
                               label=linear_p_value))
  gg.a1.r <- gg.a1.r + geom_hline(yintercept=0,color='black')
  gg.a1.r <- gg.a1.r + geom_linerange(aes(ymin=linear_ci95_low,                   
                                          ymax=linear_ci95_high),             
                                      position=pd,
                                      color=c(rep(COLORBLIND_BLUE, 2),
                                              COLORBLIND_GREEN, LOF_COLOR),
                                     )
  gg.a1.r <- gg.a1.r + geom_point(shape = 19, size  = 4, position=pd, shape=19,
                                  color=c(rep(COLORBLIND_BLUE, 2),
                                          COLORBLIND_GREEN, LOF_COLOR),
                                 )
  gg.a1.r <- gg.a1.r + theme(text=element_text(size=12))
  gg.a1.r <- gg.a1.r + theme(legend.position = 'none',
                             axis.title.x = element_text(face = "bold", size = 13),
                             axis.title.y = element_text(face = "bold", size = 13),
                             axis.text.x = element_text(face = "bold", size = 13,
                                                        angle = 0),
                             axis.text.y = element_text(hjust=1,
                                                        face = "bold", size = 13, 
                                                        angle = 45)
                            )
  gg.a1.r <- gg.a1.r + geom_text(color="black", angle=90, size=4,
                                 position=position_nudge(x=0.2))
  gg.a1.r <- gg.a1.r + ylab("rate difference (case-ctrl)")
  gg.a1.r <- gg.a1.r + xlab("Annotation")
  ggsave(paste0(OUTROOT, ".exome.annotations.caco_rate.pdf"))

  # analysis 2 : repeat analysis 1 (de novo SNV analysis across the exome), this time
  # comparing case de novo mutation rate against sequence-based context 
  
  # establish master set of dnms : general set of de novo mutations from cases
  ca.cdnm.rate <- ca.cdnm

  # get observed and expected dnm counts per annotation
  o_e_syn <- ObsExpMutRate(ca.cdnm.rate, gene.mu, c("syn"))
  o_e_misnd <- ObsExpMutRate(ca.cdnm.rate, gene.mu, c("misB","misP"))
  o_e_misd <- ObsExpMutRate(ca.cdnm.rate, gene.mu, c("misD"))
  o_e_lof <- ObsExpMutRate(ca.cdnm.rate, gene.mu, c("non","splice","frameshift"))
  o_e_misdlof <- ObsExpMutRate(ca.cdnm.rate, gene.mu, c("misD","non","splice","frameshift"))

  # get poisson p-values for observed vs expected
  o_e_syn_p <- PoissonP(o_e_syn[1], o_e_syn[2])
  o_e_misnd_p <- PoissonP(o_e_misnd[1], o_e_misnd[2])
  o_e_misd_p <- PoissonP(o_e_misd[1], o_e_misd[2])
  o_e_lof_p <- PoissonP(o_e_lof[1], o_e_lof[2])
  o_e_misdlof_p <- PoissonP(o_e_misdlof[1], o_e_misdlof[2])

  # get confidence intervals and p-values from poisson test (one-sided)
  o_e_syn_res <- poisson.test(o_e_syn[1], n.ca.trios,
                              r=o_e_syn[2]/n.ca.trios,
                              alternative='greater')
  o_e_misnd_res <- poisson.test(o_e_misnd[1], n.ca.trios,
                                r=o_e_misnd[2]/n.ca.trios,
                                alternative='greater')
  o_e_misd_res <- poisson.test(o_e_misd[1], n.ca.trios,
                               r=o_e_misd[2]/n.ca.trios,
                               alternative='greater')
  o_e_lof_res <- poisson.test(o_e_lof[1], n.ca.trios,
                                 r=o_e_lof[2]/n.ca.trios,
                                 alternative='greater')
  o_e_misdlof_res <- poisson.test(o_e_misdlof[1], n.ca.trios,
                                  r=o_e_misdlof[2]/n.ca.trios,
                                  alternative='greater')

  # store analysis 2 results to table
  a2_df <- data.frame(Annotation=c("synon","misND","misD","LoF","misD/LoF"),
                      count_obs=c(o_e_syn[1], o_e_misnd[1], o_e_misd[1],
                                  o_e_lof[1], 
                                  o_e_misdlof[1]),
                      count_exp=c(o_e_syn[2], o_e_misnd[2], o_e_misd[2],
                                  o_e_lof[2], 
                                  o_e_misdlof[2]),
                      rate_obs=c(o_e_syn[1]/n.ca.trios,
                                 o_e_misnd[1]/n.ca.trios,
                                 o_e_misd[1]/n.ca.trios,
                                 o_e_lof[1]/n.ca.trios,
                                 o_e_misdlof[1]/n.ca.trios),
                      rate_obs_ci95_low=c(o_e_syn_res$conf.int[1],
                                          o_e_misnd_res$conf.int[1],
                                          o_e_misd_res$conf.int[1],
                                          o_e_lof_res$conf.int[1],
                                          o_e_misdlof_res$conf.int[1]),
                      rate_obs_ci95_high=c(o_e_syn_res$conf.int[2],
                                           o_e_misnd_res$conf.int[2],
                                           o_e_misd_res$conf.int[2], 
                                           o_e_lof_res$conf.int[2],
                                           o_e_misdlof_res$conf.int[2]),
                      rate_exp=c(o_e_syn[2]/n.ca.trios,
                                 o_e_misnd[2]/n.ca.trios,
                                 o_e_misd[2]/n.ca.trios,
                                 o_e_lof[2]/n.ca.trios,
                                 o_e_misdlof[2]/n.ca.trios),
                      rate_ratio=c(o_e_syn[1]/o_e_syn[2],
                                   o_e_misnd[1]/o_e_misnd[2],
                                   o_e_misd[1]/o_e_misd[2],
                                   o_e_lof[1]/o_e_lof[2],
                                   o_e_misdlof[1]/o_e_misdlof[2]),
                      poisson_p=c(o_e_syn_res$p.value,
                                  o_e_misnd_res$p.value,
                                  o_e_misd_res$p.value,
                                  o_e_lof_res$p.value,
                                  o_e_misdlof_res$p.value)
                     )

  # write analysis 2 results to csv
  write.csv(a2_df, quote=F, row.names=F,
            file=paste0(OUTROOT,".exome.annotations.rate.csv"))

  # plot analysis 2 results
  pd <- position_dodge(width = 0.4)
  a2_df$p_rr <- paste0(a2_df$count_obs, " obs / ", 
                       formatC(a2_df$count_exp, format='f', digits=2), 
                       " exp",
                       "\n",
                       RateRatioFix(a2_df$rate_ratio),
                       "\n",
                       PvalFix(a2_df$poisson_p))
  a2_df$Annotation <- factor(a2_df$Annotation, levels=a2_df$Annotation)
  a2_df$rate_obs <- a2_df$rate_obs - a2_df$rate_exp
  a2_df$rate_obs_ci95_low <- a2_df$rate_obs_ci95_low - a2_df$rate_exp
  a2_df$rate_obs_ci95_high <- a2_df$rate_obs_ci95_high - a2_df$rate_exp
  a2_df$poisson_p <- PvalFix(a2_df$poisson_p)
  gg.a2 <- ggplot(a2_df, aes(y=rate_obs, x=Annotation, label=p_rr))
  gg.a2 <- gg.a2 + geom_linerange(aes(ymin=rate_obs_ci95_low,
                                      ymax=rate_obs),
                                      position=pd )
  gg.a2 <- gg.a2 + geom_point(shape = 19, size  = 4, position=pd, shape=19)
  # position=position_nudge(x=0.4),
  gg.a2 <- gg.a2 + geom_text(position=position_nudge(x=0.1),
                             hjust=1, vjust=1,
                             angle=90,size=4)
  gg.a2 <- gg.a2 + ylab("rate (obs - exp, mutability)")
  gg.a2 <- gg.a2 + xlab("Annotation")
  gg.a2 <- gg.a2 + geom_hline(yintercept=0)
  gg.a2 <- gg.a2 + theme(legend.position = 'none',
                         axis.title.x = element_text(face = "bold", size = 13),
                         axis.title.y = element_text(face = "bold", size = 13),
                         axis.text.x = element_text(face = "bold", size = 13,
                                                    angle = 0),
                         axis.text.y = element_text(hjust=1,
                                                    face = "bold", size = 13, 
                                                    angle = 45)
                        )
  gg.a2 <- gg.a2 + ylim(-0.0825, 0.05)
  ggsave(paste0(OUTROOT, ".exome.annotations.rate.pdf"))

  # analysis 3 : case/control tests of LoF burden per LOEUF bin using the
  # framework from analysis 1.

  # establish callset to use
  ca.cdnm.caco <- CdnmMultiMutRemove(ca.cdnm, dist.thresh=200)
  co.cdnm.caco <- CdnmMultiMutRemove(co.cdnm, dist.thresh=200) 
 
  # form a single case/control data frame
  caco.cdnm <- rbind(ca.cdnm.caco, co.cdnm.caco)

  # analysis 3 count dataframe is a derivative of the analysis 1 count dataframe
  a3_count_df <- a1_count_df[,c("iid","phe","n_dnm","n_outofkits_dnm",
                                "n_snv_dnm", "n_indel_dnm","n_misd_dnm")]

  # gather LOEUF bins
  loeuf_bins <- list("1"=c(0),
                     "2 - 3"=c(1,2),
                     "4 - 6"=c(3,4,5),
                     "7 - 10"=c(6,7,8,9)
                    )
  
  # iterate through each LOEUF decile and get sample-level counts for the 
  # following in genes found within LOEUF decile :
  # 1. total num of DNMs that are in geneset but not in capture intersect loci
  # 2. total num of damaging DNMs in geneset
  #    a. loss of function (non, splice, frameshift)
  #    b. missense damaging (misD)
  annot_list <- list("lof"=c("non","splice","frameshift"),
                     "misd"=c("misD")
                    )
  for (annot in names(annot_list)) {
    annots <- annot_list[[annot]]
    for (bin_name in names(loeuf_bins)) {
      bin_name_x <- gsub(" - ", "_", bin_name)
      loeuf_deciles <- loeuf_bins[[bin_name]]
      loeuf_bin_genes <- gnomad[gnomad$oe_lof_upper_bin %in% loeuf_deciles, "gene"]
      n_outofkits_syn_loeuf_bin_i <- c()
      n_outofkits_dnm_loeuf_bin_i <- c()
      n_dnm_loeuf_bin_i <- c()
      for (iid_j in a3_count_df$iid) {
        n_loeuf_bin_i_iid_j_0 <- nrow(subset(caco.cdnm,
                                             (IID == iid_j) &
                                             (GENE %in% loeuf_bin_genes) &
                                             (INKITS == FALSE)
                                            )
                                     )
        n_loeuf_bin_i_iid_j_syn <- nrow(subset(caco.cdnm,
                                               (IID == iid_j) & 
                                               (GENE %in% loeuf_bin_genes) &
                                               (EFF %in% c("syn"))
                                              )
                                       )
        n_loeuf_bin_i_iid_j <- nrow(subset(caco.cdnm,
                                           (IID == iid_j) & 
                                           (GENE %in% loeuf_bin_genes) &
                                           (EFF %in% annots)
                                          )
                                   )
        
        n_outofkits_dnm_loeuf_bin_i <- c(n_outofkits_dnm_loeuf_bin_i,
                                         n_loeuf_bin_i_iid_j_0)
        n_outofkits_syn_loeuf_bin_i <- c(n_outofkits_syn_loeuf_bin_i,
                                         n_loeuf_bin_i_iid_j_syn)
        n_dnm_loeuf_bin_i <- c(n_dnm_loeuf_bin_i, n_loeuf_bin_i_iid_j)    
      }
      outofkits_pred_i <- paste0("n_outofkits_dnm_LOEUF_",bin_name_x)
      loeuf_pred_i <- paste0("n_",annot,"_dnm_LOEUF_",bin_name_x)
      outofkits_pred_i <- paste0("n_outofkits_dnm_LOEUF_",bin_name_x)
      a3_count_df[[outofkits_pred_i]] <- n_outofkits_dnm_loeuf_bin_i
      loeuf_pred_i <- paste0("n_syn_dnm_LOEUF_",bin_name_x)
      a3_count_df[[loeuf_pred_i]] <- n_outofkits_syn_loeuf_bin_i 
      loeuf_pred_i <- paste0("n_",annot,"_dnm_LOEUF_",bin_name_x) 
      a3_count_df[[loeuf_pred_i]] <- n_dnm_loeuf_bin_i
    }
  }

  # form result dataframe for analysis 3 (lof)
  a3_df <- data.frame(Annotation=character(),
                      mean_case=numeric(),
                      mean_ctrl=numeric(),
                      linear_estimate=numeric(),
                      linear_ci95_low=numeric(),
                      linear_ci95_high=numeric(),
                      linear_p_value=numeric(),
                      logistic_odds_ratio=numeric(),
                      logistic_ci95_low=numeric(),
                      logistic_ci95_high=numeric(),
                      logistic_p_value=numeric())
 
  # iterate through each LOEUF bin and perform case/control testing (lof)
  for (bin_name in names(loeuf_bins)) {
    bin_name_x <- gsub(" - ", "_", bin_name)
    covar_i <- paste0("n_outofkits_dnm+n_outofkits_dnm_LOEUF_", bin_name_x)
    pred_i <- paste0("n_lof_dnm_LOEUF_",bin_name_x)
    logistic_mdl_str <- paste0("phe~",covar_i,"+",pred_i)
    linear_mdl_str <- paste0(pred_i,"~",covar_i,"+phe")
    logistic_mdl_i <- glm(as.formula(logistic_mdl_str), data=a3_count_df, family=binomial)
    linear_mdl_i <- lm(as.formula(linear_mdl_str), data=a3_count_df)
    logistic_res_i <- summary(logistic_mdl_i)
    linear_res_i <- summary(linear_mdl_i)
    oddsratios_i <- exp(coef(logistic_mdl_i))

    a3_df <- rbind(a3_df,
                   data.frame(Annotation=bin_name, 
                              mean_case=mean(subset(a3_count_df, phe==1)[[pred_i]]),
                              mean_ctrl=mean(subset(a3_count_df, phe==0)[[pred_i]]),
                              linear_estimate=linear_res_i$coefficients["phe",1],
                              linear_ci95_low=confint(linear_mdl_i)["phe",1],
                              linear_ci95_high=confint(linear_mdl_i)["phe",2],
                              linear_p_value=linear_res_i$coefficients["phe",4],
                              logistic_odds_ratio=exp(coef(logistic_mdl_i)[pred_i]),
                              logistic_ci95_low=exp(confint(logistic_mdl_i)[pred_i,1]),
                              logistic_ci95_high=exp(confint(logistic_mdl_i)[pred_i,2]),
                              logistic_p_value=logistic_res_i$coefficients[pred_i,4])
                   )
 
  }

  # form result dataframe for analysis 3 (lof)
  a3_misd_df <- data.frame(Annotation=character(),
                           mean_case=numeric(),
                           mean_ctrl=numeric(),
                           linear_estimate=numeric(),
                           linear_ci95_low=numeric(),
                           linear_ci95_high=numeric(),
                           linear_p_value=numeric(),
                           logistic_odds_ratio=numeric(),
                           logistic_ci95_low=numeric(),
                           logistic_ci95_high=numeric(),
                           logistic_p_value=numeric())
  
  # iterate through each LOEUF bin and perform case/control testing (lof)
  for (bin_name in names(loeuf_bins)) {
    bin_name_x <- gsub(" - ", "_", bin_name)
    covar_i <- paste0("n_outofkits_dnm+n_outofkits_dnm_LOEUF_", bin_name_x)
    pred_i <- paste0("n_misd_dnm_LOEUF_",bin_name_x)
    logistic_mdl_str <- paste0("phe~",covar_i,"+",pred_i)
    linear_mdl_str <- paste0(pred_i,"~",covar_i,"+phe")
    logistic_mdl_i <- glm(as.formula(logistic_mdl_str), data=a3_count_df, family=binomial)
    linear_mdl_i <- lm(as.formula(linear_mdl_str), data=a3_count_df)
    logistic_res_i <- summary(logistic_mdl_i)
    linear_res_i <- summary(linear_mdl_i)
    oddsratios_i <- exp(coef(logistic_mdl_i))

    a3_misd_df <- rbind(a3_misd_df,
                        data.frame(Annotation=bin_name, 
                                   mean_case=mean(subset(a3_count_df, phe==1)[[pred_i]]),
                                   mean_ctrl=mean(subset(a3_count_df, phe==0)[[pred_i]]),
                                   linear_estimate=linear_res_i$coefficients["phe",1],
                                   linear_ci95_low=confint(linear_mdl_i)["phe",1],
                                   linear_ci95_high=confint(linear_mdl_i)["phe",2],
                                   linear_p_value=linear_res_i$coefficients["phe",4],
                                   logistic_odds_ratio=exp(coef(logistic_mdl_i)[pred_i]),
                                   logistic_ci95_low=exp(confint(logistic_mdl_i)[pred_i,1]),
                                   logistic_ci95_high=exp(confint(logistic_mdl_i)[pred_i,2]),
                                   logistic_p_value=logistic_res_i$coefficients[pred_i,4])
                       )
 
  }


  # write results of analysis 3 to csv (lof)
  write.csv(a3_df, row.names=F, quote=F,
            file=paste0(OUTROOT,".loeuf_bins.lof.caco.csv"))
  
  # write results of analysis 3 to csv (misd)
  write.csv(a3_misd_df, row.names=F, quote=F,
            file=paste0(OUTROOT,".loeuf_bins.misd.caco.csv"))

  # plot results of analysis 3 using ggplot2
  # establish pd
  pd <- position_dodge(width = 0.4)

  # plot LoF in LOEUF deciles (logistic regression)
  a3_df$logistic_p_value <- PvalFix(a3_df$logistic_p_value)
  a3_df$Annotation <- factor(a3_df$Annotation,
                             levels=a3_df$Annotation )
  gg.a3 <- ggplot(a3_df, 
                  aes(y=logistic_odds_ratio, x=Annotation,
                      label=logistic_p_value))
  gg.a3 <- gg.a3 + geom_linerange(aes(ymin=logistic_ci95_low,
                                      ymax=logistic_ci95_high),
                                  position=pd, color=LOF_COLOR )
  gg.a3 <- gg.a3 + geom_point(shape = 19, size  = 4, 
                              color=LOF_COLOR,
                              position=pd, shape=19)
  gg.a3 <- gg.a3 + geom_text(position=position_nudge(x=0.2),
                             angle=90,size=4)
  gg.a3 <- gg.a3 + ylab("case/control odds ratio")
  gg.a3 <- gg.a3 + xlab("LOEUF decile")
  gg.a3 <- gg.a3 + geom_hline(yintercept=1)
  gg.a3 <- gg.a3 + theme(legend.position = 'none')
  gg.a3 <- gg.a3 + theme(axis.title.x = element_text(face = "bold", size = 14),
                         axis.title.y = element_text(face = "bold", size = 14),
                         axis.text.x = element_text(face = "bold", size = 14,
                                                    angle = 0),
                         axis.text.y = element_text(face = "bold", size = 14, 
                                                    angle = 45)
                        )
  ggsave(paste0(OUTROOT, ".loeuf_bins.lof.caco.pdf"))

  # make additional figure from analysis depicting case/control rate diffs
  a3_df$linear_p_value <- PvalFix(a3_df$linear_p_value)
  gg.a3.r <- ggplot(a3_df, 
                    aes(y=linear_estimate, x=Annotation,
                        label=linear_p_value))
  gg.a3.r <- gg.a3.r + geom_linerange(aes(ymin=linear_ci95_low,
                                          ymax=linear_ci95_high),
                                  position=pd, color=LOF_COLOR )
  gg.a3.r <- gg.a3.r + geom_point(shape = 19, size  = 4, 
                                  color=LOF_COLOR,
                                  position=pd, shape=19)
  gg.a3.r <- gg.a3.r + geom_text(position=position_nudge(x=0.2),
                                 angle=90,size=4)
  gg.a3.r <- gg.a3.r + ylab("rate difference (case-ctrl)")
  gg.a3.r <- gg.a3.r + xlab("LOEUF decile")
  gg.a3.r <- gg.a3.r + geom_hline(yintercept=0)
  gg.a3.r <- gg.a3.r + theme(legend.position = 'none')
  gg.a3.r <- gg.a3.r + theme(axis.title.x = element_text(face = "bold", size = 14),
                             axis.title.y = element_text(face = "bold", size = 14),
                             axis.text.x = element_text(face = "bold", size = 14,
                                                        angle = 0),
                             axis.text.y = element_text(face = "bold", size = 14, 
                                                        angle = 45)
                            )
  ggsave(paste0(OUTROOT, ".loeuf_bins.lof.caco_rate.pdf"))



  # roll analysis 1 and analysis 3 odds ratio figures into one figure
  pdf(paste0(OUTROOT, ".caco.pdf"))
  multiplot(gg.a1, gg.a3, cols=1)
  dev.off()

  # roll analysis 1 and analysis rate diff figures into one figure
  pdf(paste0(OUTROOT, ".caco_rate.pdf"))
  multiplot(gg.a1.r, gg.a3.r, cols=1)
  dev.off()  

  # do the above for misd variants and roll into into one figure as well
  a3_misd_df$logistic_p_value <- PvalFix(a3_misd_df$logistic_p_value)
  a3_misd_df$Annotation <- factor(a3_misd_df$Annotation,
                                  levels=a3_misd_df$Annotation )
  gg.a3 <- ggplot(a3_misd_df, 
                  aes(y=logistic_odds_ratio, x=Annotation,
                      label=logistic_p_value))
  gg.a3 <- gg.a3 + geom_linerange(aes(ymin=logistic_ci95_low,
                                      ymax=logistic_ci95_high),
                                  position=pd, color=COLORBLIND_GREEN )
  gg.a3 <- gg.a3 + geom_point(shape = 19, size  = 4, 
                              color=COLORBLIND_GREEN,
                              position=pd, shape=19)
  gg.a3 <- gg.a3 + geom_text(position=position_nudge(x=0.2),
                             angle=90,size=4)
  gg.a3 <- gg.a3 + ylab("case/control odds ratio")
  gg.a3 <- gg.a3 + xlab("LOEUF decile")
  gg.a3 <- gg.a3 + geom_hline(yintercept=1)
  gg.a3 <- gg.a3 + theme(legend.position = 'none')
  gg.a3 <- gg.a3 + theme(axis.title.x = element_text(face = "bold", size = 14),
                         axis.title.y = element_text(face = "bold", size = 14),
                         axis.text.x = element_text(face = "bold", size = 14,
                                                    angle = 0),
                         axis.text.y = element_text(face = "bold", size = 14, 
                                                    angle = 45)
                        )
  a3_misd_df$linear_p_value <- PvalFix(a3_misd_df$linear_p_value)
  gg.a3.r <- ggplot(a3_misd_df, 
                    aes(y=linear_estimate, x=Annotation,
                        label=linear_p_value))
  gg.a3.r <- gg.a3.r + geom_linerange(aes(ymin=linear_ci95_low,
                                          ymax=linear_ci95_high),
                                  position=pd, color=COLORBLIND_GREEN )
  gg.a3.r <- gg.a3.r + geom_point(shape = 19, size  = 4, 
                                  color=COLORBLIND_GREEN,
                                  position=pd, shape=19)
  gg.a3.r <- gg.a3.r + geom_text(position=position_nudge(x=0.2),
                                 angle=90,size=4)
  gg.a3.r <- gg.a3.r + ylab("rate difference (case-ctrl)")
  gg.a3.r <- gg.a3.r + xlab("LOEUF decile")
  gg.a3.r <- gg.a3.r + geom_hline(yintercept=0)
  gg.a3.r <- gg.a3.r + theme(legend.position = 'none')
  gg.a3.r <- gg.a3.r + theme(axis.title.x = element_text(face = "bold", size = 14),
                             axis.title.y = element_text(face = "bold", size = 14),
                             axis.text.x = element_text(face = "bold", size = 14,
                                                        angle = 0),
                             axis.text.y = element_text(face = "bold", size = 14, 
                                                        angle = 45)
                            )
  pdf(paste0(OUTROOT, ".loeuf_bins.misd.caco_caco_rate.pdf"))
  multiplot(gg.a3, gg.a3.r, cols=1)
  dev.off()  

  # analysis 4 : loss of function variants
  # observed vs expected mutation rate across LOEUF bins

  # establish master set of dnms : general set of de novo mutations from cases
  ca.cdnm.rate <- ca.cdnm

  # init analysis 4 result data frame
  a4_df <- data.frame(Annotation=character(),
                      count_obs=numeric(),
                      count_exp=numeric(),
                      rate_obs=numeric(),
                      rate_obs_ci95_low=numeric(),
                      rate_obs_ci95_high=numeric(),
                      rate_exp=numeric(),
                      rate_ratio=numeric(),
                      poisson_p=numeric())
  a4_denovolyzeR_df <- data.frame(Annotation=character(),
                                  observed=numeric(), 
                                  expected=numeric(), 
                                  enrichment=numeric(),
                                  pValue=numeric())

  # gather LOEUF bins
  loeuf_bins <- list("1"=c(0),
                     "2 - 3"=c(1,2),
                     "4 - 5"=c(3,4),
                     "6 - 7"=c(5,6),
                     "8 - 10"=c(7,8,9)
                    )

  # for each LOEUF decile, get the observed and expected number of LoF DNMs
  for (bin_name in names(loeuf_bins)) {
    bin_name_x <- gsub(" - ", "_", bin_name)
    loeuf_deciles <- loeuf_bins[[bin_name]]
    loeuf_bin_genes <- gnomad[gnomad$oe_lof_upper_bin %in% loeuf_deciles, "gene"]
 
    # subset cdnm and gene mu table on genes in loeuf bin
    ca.cdnm.rate.d <- subset(ca.cdnm.rate, GENE %in% loeuf_bin_genes)
    gene.mu.d <- subset(gene.mu, Gene %in% loeuf_bin_genes)

    # get observed and expected dnm counts
    o_e_syn <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d, c("syn")) 
    o_e_lofsnv <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d, c("non","splice"))
    o_e_lof <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d,
                             c("non","splice","frameshift"))

    # get poisson p-values for observed vs expected  
    o_e_syn_p <- PoissonP(o_e_syn[1], o_e_syn[2])
    o_e_lofsnv_p <- PoissonP(o_e_lofsnv[1], o_e_lofsnv[2])
    o_e_lof_p <- PoissonP(o_e_lof[1], o_e_lof[2])

    # one-sided poisson test to get p-value for observed vs expected
    o_e_syn_res <- poisson.test(o_e_syn[1], n.ca.trios, 
                                r=o_e_syn[2]/n.ca.trios,
                                alternative='greater')
    o_e_lofsnv_res <- poisson.test(o_e_lofsnv[1], n.ca.trios,
                                   r=o_e_lofsnv[2]/n.ca.trios,
                                   alternative='greater') 
    o_e_lof_res <- poisson.test(o_e_lof[1], n.ca.trios, 
                                r=o_e_lof[2]/n.ca.trios,
                                alternative='greater')     

    # get subset of denovolyzer table with genes from LOEUF bin
    denovolyzer.ptbl.d <- subset(denovolyzer.ptbl, 
                                 hgncSymbol %in% loeuf_bin_genes)

    # get denovolyzeR result
    ca.cdnm.rate.d$EFF <- as.character(ca.cdnm.rate.d$EFF)
    ca.cdnm.rate.d.dr <- subset(ca.cdnm.rate.d,
                                EFF %in% c("non","splice","frameshift"))
    ca.cdnm.rate.d.dr$GENE <- as.character(ca.cdnm.rate.d.dr$GENE)
    ca.cdnm.rate.d.dr$EFF <- as.character(ca.cdnm.rate.d.dr$EFF)
    if (nrow(ca.cdnm.rate.d.dr) > 0) {
      dr_res.d <- denovolyze(ca.cdnm.rate.d.dr$GENE, ca.cdnm.rate.d.dr$EFF,
                             nsamples=n.ca.trios, includeClasses=c("lof"),
                             probTable=denovolyzer.ptbl.d)
      dr_res.d$class <- NULL
      dr_res.d <- cbind(data.frame(Annotation=bin_name),
                        dr_res.d)
    
      # add denovolyzeR results to df
      a4_denovolyzeR_df <- rbind(a4_denovolyzeR_df, dr_res.d)
    } else {
      a4_denovolyzeR_df <- rbind(a4_denovolyzeR_df,
                                 data.frame(Annotation=bin_name,
                                            observed=0, 
                                            expected=0, 
                                            enrichment=0,
                                            pValue=1)
                                )
    }

    # add results to df 
    a4_df <- rbind(a4_df,
                   data.frame(Annotation=bin_name,
                              count_obs=o_e_lof[1],
                              count_exp=o_e_lof[2],
                              rate_obs=o_e_lof[1]/n.ca.trios,
                              rate_obs_ci95_low=o_e_lof_res$conf.int[1],
                              rate_obs_ci95_high=o_e_lof_res$conf.int[2], 
                              rate_exp=o_e_lof[2]/n.ca.trios,
                              rate_ratio=o_e_lof[1]/o_e_lof[2],
                              poisson_p=o_e_lof_res$p.value)
                  )
  }  

  # write results for analysis 4 to csv
  write.csv(a4_df, row.names=F, quote=F,
            file=paste0(OUTROOT,".loeuf_bins.lof.rate.csv"))
 
  # write denovolyzeR results from analysis 4 to csv
  write.csv(a4_denovolyzeR_df, row.names=F, quote=F,
            file=paste0(OUTROOT,".loeuf_bins.lof.rate_denovolyzeR.csv"))

  # plot analysis 4 results
  pd <- position_dodge(width = 0.4)  
  a4_df$p_rr <- paste0(a4_df$count_obs, " obs / ", 
                       formatC(a4_df$count_exp, format='f', digits=2), 
                       " exp",
                       "\n",
                       RateRatioFix(a4_df$rate_ratio),
                       "\n",
                       PvalFix(a4_df$poisson_p))
  a4_df$Annotation <- factor(a4_df$Annotation, levels=a4_df$Annotation)
  a4_df$rate_obs <- a4_df$rate_obs - a4_df$rate_exp
  a4_df$rate_obs_ci95_low <- a4_df$rate_obs_ci95_low - a4_df$rate_exp
  a4_df$rate_obs_ci95_high <- a4_df$rate_obs_ci95_high - a4_df$rate_exp
  a4_df$poisson_p <- PvalFix(a4_df$poisson_p)
  gg.a4 <- ggplot(a4_df, aes(y=rate_obs, x=Annotation, label=p_rr))
  gg.a4 <- gg.a4 + geom_linerange(aes(ymin=rate_obs_ci95_low,
                                      ymax=rate_obs),
                                      position=pd )
  gg.a4 <- gg.a4 + geom_point(shape = 19, size  = 4, position=pd, shape=19)
  gg.a4 <- gg.a4 + geom_text(position=position_nudge(x=0.1),
                             hjust=1, vjust=1,
                             angle=90,size=4)
  gg.a4 <- gg.a4 + ylab("rate (obs - exp, mutability)")
  gg.a4 <- gg.a4 + xlab("LOEUF deciles")
  gg.a4 <- gg.a4 + geom_hline(yintercept=0)
  gg.a4 <- gg.a4 + theme(legend.position = 'none',
                         axis.title.x = element_text(face = "bold", size = 13),
                         axis.title.y = element_text(face = "bold", size = 13),
                         axis.text.x = element_text(face = "bold", size = 13,
                                                    angle = 0),
                         axis.text.y = element_text(hjust=1,
                                                    face = "bold", size = 13, 
                                                    angle = 45)
                        )
 
  gg.a4 <- gg.a4 + ylim(-0.015,0.011)
  ggsave(paste0(OUTROOT, ".loeuf_bins.lof.rate.pdf"))

  # roll analysis 1 and analysis 3 figures into one figure
  pdf(paste0(OUTROOT, ".rate.pdf"))
  multiplot(gg.a2, gg.a4, cols=1)
  dev.off()

  # analysis 5 : run caco analysis across MPC bins (0-1, 1-2, >2)
  # does the burden of misD variants specifically partition to high-MPC?
  a5_count_df <- a3_count_df[,c("iid","phe","n_dnm","n_outofkits_dnm",
                                "n_snv_dnm", "n_indel_dnm",
                                "n_outofkits_dnm_LOEUF_1",
                                "n_lof_dnm_LOEUF_1")] 


  # define cdnm file to use
  ca.cdnm.caco <- CdnmMultiMutRemove(ca.cdnm, dist.thresh=200)
  co.cdnm.caco <- CdnmMultiMutRemove(co.cdnm, dist.thresh=200) 
  caco.cdnm <- rbind(ca.cdnm.caco, co.cdnm.caco)
  gnomad_non_neuro_varids <- scan(GNOMAD_NONNEURO_VARLIST,
                                  what=character())
  caco.cdnm <- subset(caco.cdnm, (VARID %in% gnomad_non_neuro_varids)==F)
  rm(gnomad_non_neuro_varids)

  # split into mpc bins based on MPC files that are read in
  mpcgt2.varids <- scan(MPCGT2_VARLIST, what=character())
  caco.cdnm.mpcgt2 <- subset(caco.cdnm, VARID %in% mpcgt2.varids)
  rm(mpcgt2.varids)
  mpc1to2.varids <- scan(MPC1TO2_VARLIST, what=character())
  caco.cdnm.mpc1to2 <- subset(caco.cdnm, VARID %in% mpc1to2.varids)
  rm(mpc1to2.varids)
  mpc0to1.varids <- scan(MPC0TO1_VARLIST, what=character())
  caco.cdnm.mpc0to1 <- subset(caco.cdnm, VARID %in% mpc0to1.varids)
  rm(mpc0to1.varids)

  # create list pointing to mpc cdnm bins
  mpc_bins <- list("0 - 1"=caco.cdnm.mpc0to1,
                   "1 - 2"=caco.cdnm.mpc1to2,
                   "> 2"=caco.cdnm.mpcgt2)

  # for each individual, get misD counts in each mpc bin
  syn.counts <- c()
  misnd.counts <- c()
  misd.counts <- c()
  mpc0to1.counts <- c()
  mpc1to2.counts <- c()
  mpcgt2.counts <- c()
  for (iid_i in a5_count_df$iid) {

    # get counts for this iid (synonymous, missense non-damaging)
    n_syn_dnm_i <- nrow(subset(caco.cdnm, 
                               (IID == iid_i) & (EFF=="syn")))
    n_misnd_dnm_i <- nrow(subset(caco.cdnm,
                                 (IID == iid_i) &
                                 (EFF %in% c("mis","misB","misP"))
                                )
                         )
    n_misd_dnm_i <- nrow(subset(caco.cdnm, 
                                (IID == iid_i) & (EFF=="misD")))
    
    # get counts for this iid for missense damaging partitioned into
    # MPC bins : low_constraint (0-1), medium constraint (1-2),
    # high constraint (>2)
    n_misdmpc0to1_dnm_i <- nrow(subset(caco.cdnm.mpc0to1,
                                       (IID == iid_i) & (EFF=="misD")))
    n_misdmpc1to2_dnm_i <- nrow(subset(caco.cdnm.mpc1to2,
                                       (IID == iid_i) & (EFF=="misD")))
    n_misdmpcgt2_dnm_i <- nrow(subset(caco.cdnm.mpcgt2,
                                      (IID == iid_i) & (EFF=="misD")))
    
    # add to arrays
    syn.counts <- c(syn.counts, n_syn_dnm_i)
    misnd.counts <- c(misnd.counts, n_misnd_dnm_i)
    misd.counts <- c(misd.counts, n_misd_dnm_i)
    mpc0to1.counts <- c(mpc0to1.counts, n_misdmpc0to1_dnm_i)
    mpc1to2.counts <- c(mpc1to2.counts, n_misdmpc1to2_dnm_i)
    mpcgt2.counts <- c(mpcgt2.counts, n_misdmpcgt2_dnm_i)
    
  }

  # add arrays to df
  a5_count_df[["n_syn_dnm"]]=syn.counts
  a5_count_df[["n_misnd_dnm"]]=misnd.counts  
  a5_count_df[["n_misd_dnm"]]=misd.counts  
  a5_count_df[["n_misdmpc0to1_dnm"]]=mpc0to1.counts 
  a5_count_df[["n_misdmpc1to2_dnm"]]=mpc1to2.counts
  a5_count_df[["n_misdmpcgt2_dnm"]]=mpcgt2.counts

  # form output df for analyses
  a5_df <- data.frame(Annotation=character(),
                      mean_case=numeric(),
                      mean_ctrl=numeric(),
                      linear_estimate=numeric(),
                      linear_ci95_low=numeric(),
                      linear_ci95_high=numeric(),
                      linear_p_value=numeric(),
                      logistic_odds_ratio=numeric(),
                      logistic_ci95_low=numeric(),
                      logistic_ci95_high=numeric(),
                      logistic_p_value=numeric())
 
  # assess case/control burden in presumably benign annotation (syn, misND)
  # also include misD for baseline measurement of misD burden for variants
  # absent from gnomad nonneuro.
  for (ann in c("syn", "misnd","misd")) {
    covar_i <- paste0("n_outofkits_dnm")
    pred_i <- paste0("n_",ann,"_dnm")
    logistic_mdl_str <- paste0("phe~",covar_i,"+",pred_i)
    linear_mdl_str <- paste0(pred_i,"~",covar_i,"+phe")
    logistic_mdl_i <- glm(as.formula(logistic_mdl_str), data=a5_count_df, family=binomial)
    linear_mdl_i <- lm(as.formula(linear_mdl_str), data=a5_count_df)
    logistic_res_i <- summary(logistic_mdl_i)
    linear_res_i <- summary(linear_mdl_i)
    oddsratios_i <- exp(coef(logistic_mdl_i))

    a5_df <- rbind(a5_df,
                   data.frame(Annotation=ann, 
                              mean_case=mean(subset(a5_count_df, phe==1)[[pred_i]]),
                              mean_ctrl=mean(subset(a5_count_df, phe==0)[[pred_i]]),
                              linear_estimate=linear_res_i$coefficients["phe",1],
                              linear_ci95_low=confint(linear_mdl_i)["phe",1],
                              linear_ci95_high=confint(linear_mdl_i)["phe",2],
                              linear_p_value=linear_res_i$coefficients["phe",4],
                              logistic_odds_ratio=exp(coef(logistic_mdl_i)[pred_i]),
                              logistic_ci95_low=exp(confint(logistic_mdl_i)[pred_i,1]),
                              logistic_ci95_high=exp(confint(logistic_mdl_i)[pred_i,2]),
                              logistic_p_value=logistic_res_i$coefficients[pred_i,4])
                   )
  } 

  # assess case/control burden across MPC bins
  for (bin_name in names(mpc_bins)) {
    bin_name_x <- gsub(" - ", "to", bin_name)
    bin_name_x <- gsub("> ", "gt", bin_name_x)
    covar_i <- paste0("n_outofkits_dnm")
    pred_i <- paste0("n_misdmpc",bin_name_x,"_dnm")
    logistic_mdl_str <- paste0("phe~",covar_i,"+",pred_i)
    linear_mdl_str <- paste0(pred_i,"~",covar_i,"+phe")
    logistic_mdl_i <- glm(as.formula(logistic_mdl_str), data=a5_count_df, family=binomial)
    linear_mdl_i <- lm(as.formula(linear_mdl_str), data=a5_count_df)
    logistic_res_i <- summary(logistic_mdl_i)
    linear_res_i <- summary(linear_mdl_i)
    oddsratios_i <- exp(coef(logistic_mdl_i))

    a5_df <- rbind(a5_df,
                   data.frame(Annotation=paste0("misDmpc_",bin_name_x), 
                              mean_case=mean(subset(a5_count_df, phe==1)[[pred_i]]),
                              mean_ctrl=mean(subset(a5_count_df, phe==0)[[pred_i]]),
                              linear_estimate=linear_res_i$coefficients["phe",1],
                              linear_ci95_low=confint(linear_mdl_i)["phe",1],
                              linear_ci95_high=confint(linear_mdl_i)["phe",2],
                              linear_p_value=linear_res_i$coefficients["phe",4],
                              logistic_odds_ratio=exp(coef(logistic_mdl_i)[pred_i]),
                              logistic_ci95_low=exp(confint(logistic_mdl_i)[pred_i,1]),
                              logistic_ci95_high=exp(confint(logistic_mdl_i)[pred_i,2]),
                              logistic_p_value=logistic_res_i$coefficients[pred_i,4])
                   )
 
  }

  # repeat test of LOEUF1 LOF using same covar
  covar_i <- paste0("n_outofkits_dnm+n_outofkits_dnm_LOEUF_1")
  pred_i <- "n_lof_dnm_LOEUF_1"
  logistic_mdl_str <- paste0("phe~",covar_i,"+",pred_i)
  linear_mdl_str <- paste0(pred_i,"~",covar_i,"+phe")
  logistic_mdl_i <- glm(as.formula(logistic_mdl_str), data=a5_count_df, family=binomial)
  linear_mdl_i <- lm(as.formula(linear_mdl_str), data=a5_count_df)
  logistic_res_i <- summary(logistic_mdl_i)
  linear_res_i <- summary(linear_mdl_i)
  oddsratios_i <- exp(coef(logistic_mdl_i))
  a5_df <- rbind(a5_df,
                 data.frame(Annotation="lofLOEUF1", 
                            mean_case=mean(subset(a5_count_df, phe==1)[[pred_i]]),
                            mean_ctrl=mean(subset(a5_count_df, phe==0)[[pred_i]]),
                            linear_estimate=linear_res_i$coefficients["phe",1],
                            linear_ci95_low=confint(linear_mdl_i)["phe",1],
                            linear_ci95_high=confint(linear_mdl_i)["phe",2],
                            linear_p_value=linear_res_i$coefficients["phe",4],
                            logistic_odds_ratio=exp(coef(logistic_mdl_i)[pred_i]),
                            logistic_ci95_low=exp(confint(logistic_mdl_i)[pred_i,1]),
                            logistic_ci95_high=exp(confint(logistic_mdl_i)[pred_i,2]),
                            logistic_p_value=logistic_res_i$coefficients[pred_i,4])
                 )

  # form a final single test combining misD MPC>2 and LOEUF1 LOF 
  col_x <- "n_lofloeuf1_misdmpcgt2_dnm"
  a5_count_df[[col_x]] = a5_count_df[["n_misdmpcgt2_dnm"]] +
                         a5_count_df[["n_lof_dnm_LOEUF_1"]]
  covar_i <- paste0("n_outofkits_dnm+n_outofkits_dnm_LOEUF_1")
  pred_i <- col_x
  logistic_mdl_str <- paste0("phe~",covar_i,"+",pred_i)
  linear_mdl_str <- paste0(pred_i,"~",covar_i,"+phe")
  logistic_mdl_i <- glm(as.formula(logistic_mdl_str), data=a5_count_df, family=binomial)
  linear_mdl_i <- lm(as.formula(linear_mdl_str), data=a5_count_df)
  logistic_res_i <- summary(logistic_mdl_i)
  linear_res_i <- summary(linear_mdl_i)
  oddsratios_i <- exp(coef(logistic_mdl_i))
  a5_df <- rbind(a5_df,
                 data.frame(Annotation="misDmpcgt2_lofLOEUF1", 
                            mean_case=mean(subset(a5_count_df, phe==1)[[pred_i]]),
                            mean_ctrl=mean(subset(a5_count_df, phe==0)[[pred_i]]),
                            linear_estimate=linear_res_i$coefficients["phe",1],
                            linear_ci95_low=confint(linear_mdl_i)["phe",1],
                            linear_ci95_high=confint(linear_mdl_i)["phe",2],
                            linear_p_value=linear_res_i$coefficients["phe",4],
                            logistic_odds_ratio=exp(coef(logistic_mdl_i)[pred_i]),
                            logistic_ci95_low=exp(confint(logistic_mdl_i)[pred_i,1]),
                            logistic_ci95_high=exp(confint(logistic_mdl_i)[pred_i,2]),
                            logistic_p_value=logistic_res_i$coefficients[pred_i,4])
                 )

  # write to csv
  write.csv(a5_df, row.names=F, quote=F,
            file=paste0(OUTROOT, ".exome.misD_MPC_bins.caco.csv"))

  # plot analysis 5 results
  a5_df <- subset(a5_df,
                  Annotation %in% 
                  c("syn","misnd","misDmpc_0to1","misDmpc_1to2",
                    "misDmpc_gt2",
                    "lofLOEUF1", 
                    "misDmpcgt2_lofLOEUF1")
                 ) 
  a5_df$Annotation <- c("synon","misND","misD,\nMPC 0-1",
                        "misD,\nMPC 1-2", "misD,\nMPC >2",
                        "LoF,\nLOEUF 1", "misD, MPC >2\nLoF, LOEUF 1")
  
  # remove syn and misND from table
  a5_df<-subset(a5_df, (Annotation %in% c("synon","misND"))==F)
  
  # assign colors
  colors <- c(rep(COLORBLIND_BLUE, 2),
              rep(COLORBLIND_GREEN, 3),
              LOF_COLOR, DAMAGING_COLOR)
  colors <- c(rep(COLORBLIND_GREEN, 3),
              LOF_COLOR, DAMAGING_COLOR)

  pd <- position_dodge(width = 0.9)
  a5_df$logistic_odds_ratio <- as.numeric(a5_df$logistic_odds_ratio)
  a5_df$Annotation <- factor(a5_df$Annotation, 
                              levels=a5_df$Annotation)
  a5_df$logistic_p_value <- PvalFix(a5_df$logistic_p_value)
  gg.a5 <- ggplot(a5_df, aes(x=Annotation, y=logistic_odds_ratio, 
                             label=logistic_p_value))
  gg.a5 <- gg.a5 + geom_hline(yintercept=1,color='black')
  gg.a5 <- gg.a5 + geom_linerange(aes(ymin=logistic_ci95_low,                   
                                      ymax=logistic_ci95_high),             
                                      position=pd,
                                      color=colors
                                 )
  gg.a5 <- gg.a5 + geom_point(shape = 19, size  = 4, 
                              position=pd, shape=19,
                              color=colors
                             )
  gg.a5 <- gg.a5 + theme(text=element_text(size=12))
  gg.a5 <- gg.a5 + theme(legend.position = 'none',
                         axis.title.x = element_text(face = "bold", size = 13),
                         axis.title.y = element_text(face = "bold", size = 13),
                         axis.text.x = element_text(face = "bold", size = 13,
                                                    angle = 0),
                         axis.text.y = element_text(hjust=1,
                                                    face = "bold", size = 13, 
                                                    angle = 45)
                        )
  gg.a5 <- gg.a5 + geom_text(color="black", angle=90, size=4,
                             position=position_nudge(x=0.2))
  gg.a5 <- gg.a5 + ylab("case/control odds ratio")
  gg.a5 <- gg.a5 + xlab("Annotation")
  ggsave(paste0(OUTROOT, ".exome.misD_MPC_bins.caco.pdf"))

  # analysis 6 : is there a difference between observed and expected rate of
  # lof and missense damaging de novos in 1) 124 NDD genes from Coe et al 2019,
  # 2) 9 TS genes from Wang 2018 that have q<0.3, 3) 199 genes from TS Wang
  # 2018 that have q>=0.3?

  # read genesets
  ndd <- scan(NDD_GENESET, what=character())
  ts2 <- scan(TS_GE2DNM_GENESET, what=character())
  ts1 <- scan(TS_1DNM_GENESET, what=character())

  # base cdnm
  ca.cdnm.rate <- ca.cdnm

  # get observed / expected lof + missense damaging counts per geneset
  o_e_ndd <- ObsExpMutRate(subset(ca.cdnm.rate, GENE %in% ndd),
                           subset(gene.mu, Gene %in% ndd), 
                           c("misD","non","splice","frameshift"))
  o_e_ts2 <- ObsExpMutRate(subset(ca.cdnm.rate, GENE %in% ts2),
                           subset(gene.mu, Gene %in% ts2), 
                           c("misD","non","splice","frameshift"))
  o_e_ts1 <- ObsExpMutRate(subset(ca.cdnm.rate, GENE %in% ts1),
                           subset(gene.mu, Gene %in% ts1),
                           c("misD","non","splice","frameshift"))

  # one-sided poisson test to get p-value for observed vs expected
  o_e_ndd_res <- poisson.test(o_e_ndd[1], n.ca.trios, 
                              r=o_e_ndd[2]/n.ca.trios,
                              alternative='greater')
  o_e_ts2_res <- poisson.test(o_e_ts2[1], n.ca.trios,
                              r=o_e_ts2[2]/n.ca.trios,
                              alternative='greater') 
  o_e_ts1_res <- poisson.test(o_e_ts1[1], n.ca.trios, 
                              r=o_e_ts1[2]/n.ca.trios,
                              alternative='greater')     

  # form result df
  genesets <- c("NDD (124 genes)", "TS, >=2 damaging DNMs (6 genes)",
              "TS, 1 damaging DNM (199 genes")
  a6_df <- data.frame(Geneset=genesets,
                      count_obs=c(o_e_ndd[1], o_e_ts2[1], o_e_ts1[1]),
                      count_exp=c(o_e_ndd[2], o_e_ts2[2], o_e_ts1[2]),
                      rate_obs=c(o_e_ndd[1], o_e_ts2[1], o_e_ts1[1])/n.ca.trios,
                      rate_obs_ci95_low=c(o_e_ndd_res$conf.int[1],
                                          o_e_ts2_res$conf.int[1], 
                                          o_e_ts1_res$conf.int[1]),
                      rate_obs_ci95_high=c(o_e_ndd_res$conf.int[2],
                                           o_e_ts2_res$conf.int[2],
                                           o_e_ts1_res$conf.int[2]),
                      rate_exp=c(o_e_ndd[2], o_e_ts2[2], o_e_ts1[2])/n.ca.trios,
                      rate_ratio=c(o_e_ndd[1]/o_e_ndd[2],
                                   o_e_ts2[1]/o_e_ts2[2],
                                   o_e_ts1[1]/o_e_ts1[2]),
                      poisson_p=c(o_e_ndd_res$p.value,
                                  o_e_ts2_res$p.value,  
                                  o_e_ts1_res$p.value)
                     )
          
  # write results to csv
  write.csv(a6_df, row.names=F, quote=F,
            file=paste0(OUTROOT, ".genesets.misdlof.rate.csv"))

  # plot analysis 5 results
  pd <- position_dodge(width = 0.4)  
  a6_df$Geneset <- c('Neurodev.\nDisorders\n(187 genes)', 
                     'Tourettes,\n>1 damaging DNM\n(6 genes)',
                     'Tourettes,\n1 damaging DNM\n(199 genes)')
  a6_df$p_rr <- paste0(a6_df$count_obs, " obs vs. ", 
                       formatC(a6_df$count_exp, format='f', digits=2), 
                       " exp",
                       "\n",
                       RateRatioFix(a6_df$rate_ratio),
                       "\n",
                       PvalFix(a6_df$poisson_p))
  a6_df$Geneset <- factor(a6_df$Geneset, levels=a6_df$Geneset)
  a6_df$rate_obs <- a6_df$rate_obs - a6_df$rate_exp
  a6_df$rate_obs_ci95_low <- a6_df$rate_obs_ci95_low - a6_df$rate_exp
  a6_df$rate_obs_ci95_high <- a6_df$rate_obs_ci95_high - a6_df$rate_exp
  a6_df$poisson_p <- PvalFix(a6_df$poisson_p)
  gg.a6 <- ggplot(a6_df, aes(y=rate_obs, x=Geneset, label=p_rr))
  gg.a6 <- gg.a6 + geom_linerange(aes(ymin=rate_obs_ci95_low,
                                      ymax=rate_obs),
                                      position=pd )
  gg.a6 <- gg.a6 + geom_point(shape = 19, size  = 4, position=pd, shape=19)
  gg.a6 <- gg.a6 + geom_text(hjust=1, vjust=1,
                             position=position_nudge(x=0.2),
                             angle=90,size=4)
  gg.a6 <- gg.a6 + ylab("rate (obs - exp, mutability)")
  gg.a6 <- gg.a6 + xlab("geneset")
  gg.a6 <- gg.a6 + geom_hline(yintercept=0)
  gg.a6 <- gg.a6 + theme(legend.position = 'none',
                         axis.title.x = element_blank(),
                         axis.title.y = element_text(face = "bold", size = 13),
                         axis.text.x = element_text(face = "bold", size = 13,
                                                    angle = 0),
                         axis.text.y = element_text(hjust=1,
                                                    face = "bold", size = 13, 
                                                    angle = 45)
                        )
 
  ggsave(paste0(OUTROOT, ".genesets.misdlof.rate.pdf"))

  # analysis 7 : loss of function variants
  # observed vs expected mutation rate across LOEUF bins
  # (essentially a repeat of analysis 4 but stratify by LoF SNVs, 
  #  LoF frameshift indels)

  # establish master set of dnms : general set of de novo mutations from cases
  ca.cdnm.rate <- ca.cdnm

  # init analysis 4 result data frame
  a7_df <- data.frame(loeuf_bin=character(),
                      annotation=character(),
                      count_obs=numeric(),
                      count_exp=numeric(),
                      rate_obs=numeric(),
                      rate_obs_ci95_low=numeric(),
                      rate_obs_ci95_high=numeric(),
                      rate_exp=numeric(),
                      rate_ratio=numeric(),
                      poisson_p=numeric())

  # gather LOEUF bins
  loeuf_bins <- list("1"=c(0),
                     "2 - 10"=c(1,2,3,4,5,6,7,8,9)
                    )

  # for each LOEUF decile, get the observed and expected number of LoF DNMs
  for (bin_name in names(loeuf_bins)) {
    bin_name_x <- gsub(" - ", "_", bin_name)
    loeuf_deciles <- loeuf_bins[[bin_name]]
    loeuf_bin_genes <- gnomad[gnomad$oe_lof_upper_bin %in% loeuf_deciles, "gene"]
 
    # subset cdnm and gene mu table on genes in loeuf bin
    ca.cdnm.rate.d <- subset(ca.cdnm.rate, GENE %in% loeuf_bin_genes)
    gene.mu.d <- subset(gene.mu, Gene %in% loeuf_bin_genes)

    # get observed and expected dnm counts
    o_e_syn <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d, c("syn")) 
    o_e_lofsnv <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d, c("non","splice"))
    o_e_non <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d, c("non"))
    o_e_splice <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d, c("splice"))
    o_e_fs <- ObsExpMutRate(ca.cdnm.rate.d, gene.mu.d, c("frameshift"))

    # get poisson p-values for observed vs expected  
    o_e_syn_p <- PoissonP(o_e_syn[1], o_e_syn[2])
    o_e_lofsnv_p <- PoissonP(o_e_lofsnv[1], o_e_lofsnv[2])
    o_e_non_p <- PoissonP(o_e_non[1], o_e_non[2])
    o_e_splice_p <- PoissonP(o_e_splice[1], o_e_splice[2])
    o_e_fs_p <- PoissonP(o_e_fs[1], o_e_fs[2])

    # one-sided poisson test to get p-value for observed vs expected
    o_e_syn_res <- poisson.test(o_e_syn[1], n.ca.trios, 
                                r=o_e_syn[2]/n.ca.trios,
                                alternative='greater')
    o_e_lofsnv_res <- poisson.test(o_e_lofsnv[1], n.ca.trios,
                                   r=o_e_lofsnv[2]/n.ca.trios,
                                   alternative='greater')
    o_e_non_res <- poisson.test(o_e_non[1], n.ca.trios, 
                                r=o_e_non[2]/n.ca.trios,
                                alternative='greater')
    o_e_splice_res <- poisson.test(o_e_splice[1], n.ca.trios, 
                                   r=o_e_splice[2]/n.ca.trios,
                                   alternative='greater')
    o_e_fs_res <- poisson.test(o_e_fs[1], n.ca.trios, 
                               r=o_e_fs[2]/n.ca.trios,
                               alternative='greater')

    # add results to df (lof snvs)
    a7_df <- rbind(a7_df,
                   data.frame(loeuf_bin=bin_name,
                              annotation="stopgain SNV",
                              count_obs=o_e_non[1],
                              count_exp=o_e_non[2],
                              rate_obs=o_e_non[1]/n.ca.trios,
                              rate_obs_ci95_low=o_e_non_res$conf.int[1],
                              rate_obs_ci95_high=o_e_non_res$conf.int[2], 
                              rate_exp=o_e_non[2]/n.ca.trios,
                              rate_ratio=o_e_non[1]/o_e_non[2],
                              poisson_p=o_e_non_res$p.value)
                  )
    a7_df <- rbind(a7_df,
                   data.frame(loeuf_bin=bin_name,
                              annotation="splice donor/acceptor SNV",
                              count_obs=o_e_splice[1],
                              count_exp=o_e_splice[2],
                              rate_obs=o_e_splice[1]/n.ca.trios,
                              rate_obs_ci95_low=o_e_splice_res$conf.int[1],
                              rate_obs_ci95_high=o_e_splice_res$conf.int[2], 
                              rate_exp=o_e_splice[2]/n.ca.trios,
                              rate_ratio=o_e_splice[1]/o_e_splice[2],
                              poisson_p=o_e_splice_res$p.value)
                  )

    # add results to df (frameshift indels)
    a7_df <- rbind(a7_df,
                   data.frame(loeuf_bin=bin_name,
                              annotation="frameshift indel",
                              count_obs=o_e_fs[1],
                              count_exp=o_e_fs[2],
                              rate_obs=o_e_fs[1]/n.ca.trios,
                              rate_obs_ci95_low=o_e_fs_res$conf.int[1],
                              rate_obs_ci95_high=o_e_fs_res$conf.int[2], 
                              rate_exp=o_e_fs[2]/n.ca.trios,
                              rate_ratio=o_e_fs[1]/o_e_fs[2],
                              poisson_p=o_e_fs_res$p.value)
                  )

  }
  
  # write a7 results to csv
  write.csv(a7_df, row.names=F, quote=F,
            file=paste0(OUTROOT, ".loeuf_bins.non_splice_frameshift.rate.csv")
           )

  # turn loeuf bins into factor
  a7_df$loeuf_bin <- factor(a7_df$loeuf_bin,
                            levels=c("1","2 - 10")
                           )
 
  # get p / rr strings
  a7_df$p_rr <- paste0(a7_df$count_obs, " obs / ", 
                       formatC(a7_df$count_exp, format='f', digits=2), 
                       " exp",
                       "\n",
                       RateRatioFix(a7_df$rate_ratio),
                       "\n",
                       PvalFix(a7_df$poisson_p))
 
  # get adjust observed rates based on exp 
  a7_df$rate_obs <- a7_df$rate_obs - a7_df$rate_exp
  a7_df$rate_obs_ci95_low <- a7_df$rate_obs_ci95_low - a7_df$rate_exp
  a7_df$rate_obs_ci95_high <- a7_df$rate_obs_ci95_high - a7_df$rate_exp
  a7_df$poisson_p <- PvalFix(a7_df$poisson_p)
 
  # split into LoF SNV and Frameshift indel dfs
  a7_df.1 <- subset(a7_df, annotation=="stopgain SNV")
  a7_df.2 <- subset(a7_df, annotation=="splice donor/acceptor SNV")
  a7_df.3 <- subset(a7_df, annotation=="frameshift indel")

  # plot analysis 7 results
  pd <- position_dodge(width = 0.4)  
  gg.a7.1 <- ggplot(a7_df.1, aes(y=rate_obs, x=loeuf_bin, label=p_rr))
  gg.a7.1 <- gg.a7.1 + geom_linerange(aes(ymin=rate_obs_ci95_low,
                                          ymax=rate_obs),
                                      position=pd )
  gg.a7.1 <- gg.a7.1 + geom_point(shape = 19, size  = 4, position=pd, shape=19)
  gg.a7.1 <- gg.a7.1 + geom_text(hjust=1, vjust=1,
                                 position=position_nudge(x=0.2),
                                 angle=90,size=4)
  gg.a7.1 <- gg.a7.1 + ylab("rate (obs - exp, mutab.)")
  gg.a7.1 <- gg.a7.1 + geom_hline(yintercept=0)
  gg.a7.1 <- gg.a7.1 + theme(legend.position = 'none',
                             axis.title.x = element_blank(),
                             axis.title.y = element_text(face = "bold", size = 13),
                             axis.text.x = element_text(face = "bold", size = 13,
                                                        angle = 0),
                             axis.text.y = element_text(hjust=1,
                                                        face = "bold", size = 13, 
                                                        angle = 45)
                        )
  gg.a7.1 <- gg.a7.1 + ggtitle("stopgain SNVs")
  ggsave(paste0(OUTROOT, ".loeuf_bins.stopgain_snv.rate.pdf"))
  gg.a7.2 <- ggplot(a7_df.2, aes(y=rate_obs, x=loeuf_bin, label=p_rr))
  gg.a7.2 <- gg.a7.2 + geom_linerange(aes(ymin=rate_obs_ci95_low,
                                          ymax=rate_obs),
                                      position=pd )
  gg.a7.2 <- gg.a7.2 + geom_point(shape = 19, size  = 4, position=pd, shape=19)
  gg.a7.2 <- gg.a7.2 + geom_text(hjust=1, vjust=1,
                                 position=position_nudge(x=0.2),
                                 angle=90,size=4)
  gg.a7.2 <- gg.a7.2 + ylab("rate (obs - exp, mutab.)")
  gg.a7.2 <- gg.a7.2 + geom_hline(yintercept=0)
  gg.a7.2 <- gg.a7.2 + theme(legend.position = 'none',
                             axis.title.x = element_blank(),
                             axis.title.y = element_text(face = "bold", size = 13),
                             axis.text.x = element_text(face = "bold", size = 13,
                                                        angle = 0),
                             axis.text.y = element_text(hjust=1,
                                                        face = "bold", size = 13, 
                                                        angle = 45)
                        )
  gg.a7.2 <- gg.a7.2 + ggtitle("splice donor/acceptor SNVs")
  ggsave(paste0(OUTROOT, ".loeuf_bins.splice_snv.rate.pdf"))
  gg.a7.3 <- ggplot(a7_df.3, aes(y=rate_obs, x=loeuf_bin, label=p_rr))
  gg.a7.3 <- gg.a7.3 + geom_linerange(aes(ymin=rate_obs_ci95_low,
                                          ymax=rate_obs),
                                      position=pd )
  gg.a7.3 <- gg.a7.3 + geom_point(shape = 19, size  = 4, position=pd, shape=19)
  gg.a7.3 <- gg.a7.3 + geom_text(hjust=1, vjust=1,
                                 position=position_nudge(x=0.2),
                                 angle=90,size=4)
  gg.a7.3 <- gg.a7.3 + ylab("rate (obs - exp, mutab.)")
  gg.a7.3 <- gg.a7.3 + xlab("LOEUF deciles")
  gg.a7.3 <- gg.a7.3 + geom_hline(yintercept=0)
  gg.a7.3 <- gg.a7.3 + theme(legend.position = 'none',
                             axis.title.x = element_text(face = "bold", size = 13),
                             axis.title.y = element_text(face = "bold", size = 13),
                             axis.text.x = element_text(face = "bold", size = 13,
                                                        angle = 0),
                             axis.text.y = element_text(hjust=1,
                                                        face = "bold", size = 13, 
                                                        angle = 45)
                        )
  gg.a7.3 <- gg.a7.3 + ggtitle("frameshift indels")
  gg.a7.3 <- gg.a7.3 + ylim(-0.04, 0.01)
  ggsave(paste0(OUTROOT, ".loeuf_bins.fsindel.rate.pdf"))

  # roll lofsnv and fsindel results into one figure
  pdf(paste0(OUTROOT, ".loeuf_bins.non_splice_fs.rate.pdf"))
  multiplot(gg.a7.1, gg.a7.2, gg.a7.3, cols=1)
  dev.off()

  q()

}

CdnmMultiMutRemove <- function(cdnm.i, dist.thresh=200) {
  
  # sort cdnm by chrom, pos
  cdnm.i <- cdnm.i[order(cdnm.i$CHR, cdnm.i$POS), ]

  # rename rows
  rownames(cdnm.i)<-1:nrow(cdnm.i)

  # find instances where row i and row i+1 are same indiv and < 200 nt apart
  prev.iid<-"0"
  prev.chrom<-"0"
  prev.pos<-0
  rows.rm <- c()
  for (i in rownames(cdnm.i)) {
    if (cdnm.i[i,"CHR"] == prev.chrom & cdnm.i[i,"IID"] == prev.iid) {
      if (abs(cdnm.i[i,"POS"] - prev.pos) < dist.thresh) {
        rows.rm <- c(rows.rm, c(as.character(as.numeric(i)-1), i))
      }
    }
    prev.iid<-cdnm.i[i,"IID"]
    prev.chrom<-cdnm.i[i,"CHR"]
    prev.pos<-cdnm.i[i,"POS"]
  }
  rows.rm <- unique(sort(rows.rm))
  cdnm.i <- cdnm.i[setdiff(rownames(cdnm.i), rows.rm), , drop=F]
  return(cdnm.i)
}

ObsExpMutRate <- function(cdnm, exp_mut_df, effs_qual) {
  obs_mut <- 0
  exp_mut <- 0
  for (eff.i in effs_qual) {
    cdnm.i <- subset(cdnm, EFF == eff.i)
    obs_mut <- obs_mut + nrow(cdnm.i)
    exp_mut <- exp_mut + sum(exp_mut_df[[eff.i]])
  }
  return(c(obs_mut, exp_mut))
}

PoissonP <- function(obs, exp) {
  return(1-ppois(obs-1, exp))
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
  
PvalFix <- function(pvals) {
  pvals.x <- ifelse(pvals <= 1,
                    paste0("p = ",formatC(pvals,format='e',digits=2)),
                    "")
  return(pvals.x)
}

RateRatioFix <- function(rrvals) {
  return(paste0("rate ratio = ",
                formatC(rrvals, format='f', digits=2),""))
}

if (interactive() == F) {
  Main()
}

