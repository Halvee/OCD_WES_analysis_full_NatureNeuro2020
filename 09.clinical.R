#!/usr/bin/env Rscript

library(ggplot2)

## DATA
SAMPLEPED <- "results/dnm_cohort/OCDfams_2019.cohort.trios.sampleped"
CDNM.FILE <- "results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.cdnm"
LOEUF.TXT <- "results/input/gnomad.v2.1.1.loeuf_bins.tsv"
MPC.VARLIST <- "results/input/MPCgt2.varlist.gz"
AAO.CSV <- "data/JHU/AgeOfOnset.WESprobands.20161215.csv"
CLIN.CSV <- "data/JHU/VariablesforMatt_01042017_AAO.csv"
LINKER.CSV <- "data/clinical/ocd_samples_20190324.csv"
OUTROOT <- "results/clinical/OCD.clinical"

Main <- function() {

  # read files
  ped <- read.table(SAMPLEPED, stringsAsFactors=F, header=F)
  colnames(ped) <- c("FID","IID","PID","MID","SEX","PHE","SEQTYPE","EXOMEKIT")
  cdnm <- read.table(CDNM.FILE, stringsAsFactors=F, header=F)
  colnames(cdnm) <- c("IID","CHROM","POS","VARID","REF","ALT","GENE","EFF")
  loeuf <- read.table(LOEUF.TXT, header=T, sep="\t", stringsAsFactors=F)
  mpc.vars <- scan(MPC.VARLIST, what=character(), quiet=T)
  lnk <- read.csv(LINKER.CSV, stringsAsFactors=F)
  clin <- read.csv(CLIN.CSV, stringsAsFactors=F)

  # form CHGVID/OrigID linker
  lnk <- lnk[, c("CHGVID","OrigID")]
  colnames(lnk) <- c("IID","ID")

  # subset on probands in ped
  ped <- subset(ped, PHE == 2)

  # merge ped and lnk
  ped <- merge(ped, lnk, by="IID")

  # merge clin df into ped to form seperate df
  ped.clin <- merge(ped, clin, by="ID")
  ped.clin$SEX <- ped.clin$SEX - 1

  # derive subset of dnm calls that are LoF and in an LoF-intolerant gene(LOEUF<10%)
  loeuf1<-loeuf[(loeuf$oe_lof_upper_bin==0) & (is.na(loeuf$oe_lof_upper_bin)==F), "gene"]
  cdnm.lof.loeuf1 <- subset(cdnm, 
                              (GENE %in% loeuf1) 
                              &
                              (EFF %in% c("non","splice","frameshift"))
                             )
  cat("N DNMs with LoF annotation and genes with LOEUF < 10% :", nrow(cdnm.lof.loeuf1),"\n")
  cat("N carriers of at least one LoF DNM with LOEUF < 10% :",
      length(unique(cdnm.lof.loeuf1[,1])), "\n")

  # derive subset of dnm calls with MPC>2
  cdnm.misd.mpcgt2 <- subset(cdnm,
                             (EFF == "misD")
                             &
                             (VARID %in% mpc.vars)
                            )
  cat("N DNMs with missense annotation and MPC >2 :", nrow(cdnm.misd.mpcgt2),"\n")
  cat("N carriers of at least one missense DNM with MPC >2 :",
      length(unique(cdnm.misd.mpcgt2[,1])), "\n")

  # print total DNMs
  cdnm.dmg <- rbind(cdnm.lof.loeuf1, cdnm.misd.mpcgt2)
  cat("N DNMs classified as damaging (LoF/LOEUF1, or misD/MPC>2) :", nrow(cdnm.dmg),"\n")
  cat("N carriers of at least one damaging :",
      length(unique(cdnm.dmg[,1])), "\n")
  
  # get counts of lof/loeuf1, misd/mpcgt2
  ped.clin$lof_LOEUF1 <- ifelse(ped.clin$IID %in% cdnm.lof.loeuf1[,1], 1, 0)
  ped.clin$misD_MPCgt2 <- ifelse(ped.clin$IID %in% cdnm.misd.mpcgt2[,1], 1, 0)
  ped.clin$dmg_dnm <- ifelse(ped.clin$lof_LOEUF1 == 1 |
                             ped.clin$misD_MPCgt2 == 1, 
                             1, 
                             0
                            )

  # add column to ped.clin for MALE status
  ped.clin$MALE <- ifelse(ped.clin$SEX==0, 1, 0) 

  # perform test of burden of vars on the following binary traits:
  # 1. SEX (male)
  ped.clin.sex <- subset(ped.clin, is.na(SEX)==F)
  res_df <- GetLogisticRes("dmg_dnm", "MALE", ped.clin.sex, 
                           covars=NULL, res_df=NULL)
  
  # 2. TICS 
  ped.clin.tics01 <- subset(ped.clin, is.na(TICS)==F)
  ped.clin.tics01$TICS <- ped.clin.tics01$TICS - 1
  print(summary(lm(dmg_dnm ~ SEX + TICS, data=ped.clin.tics01)))
  res_df <- GetLogisticRes("dmg_dnm", "TICS", ped.clin.tics01, 
                           covars=c("MALE"), res_df=res_df)

  # 3. SKIN-PICKING (SKIN3)
  ped.clin.skin01 <- subset(ped.clin, is.na(SKIN3)==F)
  ped.clin.skin01$SKIN3 <- ped.clin.skin01$SKIN3 - 1
  print(summary(lm(dmg_dnm ~ SEX + SKIN3, data=ped.clin.skin01)))
  res_df <- GetLogisticRes("dmg_dnm", "SKIN3", ped.clin.skin01, 
                           covars=NULL, res_df=res_df)
 
  # 4. TRICH
  ped.clin.trich01 <- subset(ped.clin, is.na(TRICH3)==F)
  ped.clin.trich01$TRICH3 <- ped.clin.trich01$TRICH3 - 1
  print(summary(lm(dmg_dnm ~ SEX + TRICH3, data=ped.clin.trich01)))
  res_df <- GetLogisticRes("dmg_dnm", "TRICH3", ped.clin.trich01, 
                           covars=NULL, res_df=res_df)

  print(res_df)

  # write results to file
  write.csv(res_df, file=paste0(OUTROOT,".features.cmp.csv"),
            row.names=F, quote=F)

  # plot results
  res_df$feature <- c("male sex","tics","skin-picking","trichotillomania")
  res_df$feature <- factor(res_df$feature, levels=res_df$feature) 
  res_df$pval <- PvalFix(res_df$p_value_logistic)
  pd <- position_dodge(width = 0.9)
  gg.res <- ggplot(res_df, aes(x=feature, y=odds_ratio_logistic, 
                               label=pval))
  gg.res <- gg.res + geom_hline(yintercept=1,color='black')
  gg.res <- gg.res + geom_linerange(aes(ymin=ci95_low_logistic,                   
                                        ymax=ci95_high_logistic),             
                                    position=pd,
                                    color='black'
                                   )
  gg.res <- gg.res + geom_point(shape = 19, size  = 4, 
                                position=pd, shape=19,
                                color='black'
                               )
  gg.res <- gg.res + theme(text=element_text(size=12))
  gg.res <- gg.res + theme(legend.position = 'none',
                           axis.title.x = element_text(face = "bold", size = 13),
                           axis.title.y = element_text(face = "bold", size = 13),
                           axis.text.x = element_text(face = "bold", size = 13,
                                                      angle = 0),
                           axis.text.y = element_text(hjust=1,
                                                      face = "bold", size = 13, 
                                                      angle = 45)
                          )
  gg.res <- gg.res + geom_text(color="black", angle=90, size=4,
                             position=position_nudge(x=0.2))
  gg.res <- gg.res + ylab("odds ratio")
  gg.res <- gg.res + xlab("clinical feature")
  ggsave(paste0(OUTROOT,".features.cmp.pdf"))

  # followup - focus on sex imbalance. Spread across LoF/LoFintol and MisD/MPC>2?
  res.df <- data.frame(dnm_classif=character(),
                       n_m_1=numeric(), n_f_1=numeric(), 
                       n_m_0=numeric(), n_f_0=numeric(),
                       perc_m_1=numeric(), perc_m_0=numeric(),
                       fet_m_or=numeric(), 
                       fet_m_or_95ci_lower=numeric(), fet_m_or_95ci_upper=numeric(),
                       fet_m_p=numeric())

  # get male/female sex proportions in samples with/without LoF/LoFintol DNM
  ped.0 <- subset(ped, (IID %in% cdnm.lof.loeuf1$IID)==F)
  ped.1 <- subset(ped, (IID %in% cdnm.lof.loeuf1$IID)==T)
  n.m.0 <- table(ped.0$SEX)[[1]]
  n.f.0 <- table(ped.0$SEX)[[2]]
  n.m.1 <- table(ped.1$SEX)[[1]]
  n.f.1 <- table(ped.1$SEX)[[2]]
  conting.tbl <- data.frame(carrier=c(n.m.1, n.f.1), noncarrier=c(n.m.0, n.f.0))
  fet.res <- fisher.test(conting.tbl)
  res.df <- rbind(res.df,
                  data.frame(dnm_classif="LOF_LOEUF1",
                             n_m_1=n.m.1, n_f_1=n.f.1,
                             n_m_0=n.m.0, n_f_0=n.f.0,
                             perc_m_1=n.m.1/(n.m.1+n.f.1), 
                             perc_m_0=n.m.0/(n.m.0+n.f.0),
                             fet_m_or=fet.res$estimate,
                             fet_m_or_95ci_lower=fet.res$conf.int[1], 
                             fet_m_or_95ci_upper=fet.res$conf.int[2],
                             fet_m_p=fet.res$p.value)
                 )

  # get male/female sex proportions in samples with/without MPC missense DNM
  ped.0 <- subset(ped, (IID %in% cdnm.misd.mpcgt2$IID)==F)
  ped.1 <- subset(ped, (IID %in% cdnm.misd.mpcgt2$IID)==T)
  n.m.0 <- table(ped.0$SEX)[[1]]
  n.f.0 <- table(ped.0$SEX)[[2]]
  n.m.1 <- table(ped.1$SEX)[[1]]
  n.f.1 <- table(ped.1$SEX)[[2]]
  conting.tbl <- data.frame(carrier=c(n.m.1, n.f.1), noncarrier=c(n.m.0, n.f.0))
  fet.res <- fisher.test(conting.tbl)
  res.df <- rbind(res.df,
                  data.frame(dnm_classif="mis_MPCgt2",
                             n_m_1=n.m.1, n_f_1=n.f.1,
                             n_m_0=n.m.0, n_f_0=n.f.0,
                             perc_m_1=n.m.1/(n.m.1+n.f.1),
                             perc_m_0=n.m.0/(n.m.0+n.f.0),
                             fet_m_or=fet.res$estimate,
                             fet_m_or_95ci_lower=fet.res$conf.int[1],
                             fet_m_or_95ci_upper=fet.res$conf.int[2],
                             fet_m_p=fet.res$p.value)
                 )

  # get male/female sex proportions in samples with/without 
  # either an LoF DNM in an LoFintol gene, or an MPC missense DNM
  qual.iids <- unique(c(cdnm.lof.loeuf1$IID, cdnm.misd.mpcgt2$IID))
  ped.0 <- subset(ped, (IID %in% qual.iids)==F)
  ped.1 <- subset(ped, (IID %in% qual.iids)==T)
  n.m.0 <- table(ped.0$SEX)[[1]]
  n.f.0 <- table(ped.0$SEX)[[2]]
  n.m.1 <- table(ped.1$SEX)[[1]]
  n.f.1 <- table(ped.1$SEX)[[2]]
  conting.tbl <- data.frame(carrier=c(n.m.1, n.f.1), noncarrier=c(n.m.0, n.f.0))
  fet.res <- fisher.test(conting.tbl)
  res.df <- rbind(res.df,
                  data.frame(dnm_classif="LOF_LOEUF1 OR mis_MPCgt2",
                             n_m_1=n.m.1, n_f_1=n.f.1,
                             n_m_0=n.m.0, n_f_0=n.f.0,
                             perc_m_1=n.m.1/(n.m.1+n.f.1),
                             perc_m_0=n.m.0/(n.m.0+n.f.0),
                             fet_m_or=fet.res$estimate,
                             fet_m_or_95ci_lower=fet.res$conf.int[1],
                             fet_m_or_95ci_upper=fet.res$conf.int[2],
                             fet_m_p=fet.res$p.value)
                 )

  # write subsetted cdnms
  write.table(cdnm.lof.loeuf1, file=paste0(OUTROOT,".lof_lofintol.cdnm"),
              row.names=F, col.names=T, sep="\t", quote=F)
  write.table(cdnm.misd.mpcgt2, file=paste0(OUTROOT,".mis_mpcgt2.cdnm"),
              row.names=F, col.names=T, sep="\t", quote=F)

  # write result dataframe to csv
  write.csv(res.df, file=paste0(OUTROOT,".sexbias.res.csv"),
            row.names=F, quote=F)

} 

GetLogisticRes <- function(outcome, pred, df, 
                           covars=NULL, res_df=NULL) {
  if (is.null(covars)) {
    mdl_str <- paste0(outcome,"~",pred)
  } else {
    mdl_str <- paste0(outcome,"~",paste(covars,collapse="+"),"+",pred)
  }
  mdl.lg <- glm(as.formula(mdl_str),data=df, family=binomial)
  res.lg <- summary(mdl.lg)
  print(res.lg$coefficients)
  res.df.i <- data.frame(predictor=pred,
                         odds_ratio_logistic=exp(coef(mdl.lg)[pred]),
                         ci95_low_logistic=exp(confint(mdl.lg)[pred,1]),
                         ci95_high_logistic=exp(confint(mdl.lg)[pred,2]),
                         p_value_logistic=res.lg$coefficients[pred,4]
                        )
  if (is.null(res_df)==F) {
    res_df <- rbind(res_df, res.df.i)
  } else {
    res_df <- res.df.i
  }
  return(res_df)
}

PvalFix <- function(pvals) {
  pvals.x <- ifelse(pvals <= 1,
                    paste0("p = ",formatC(pvals,format='e',digits=2)),
                    "")
  return(pvals.x)
}

if (interactive() == F) {
  Main()
}
