
Main <- function(){
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 4) {
    cat("qualvar_carrier_analysis.R <caco.sampleped>",
        "<case.ddb.csv> <ctrl.ddb.csv< <caco.qualvar>\n")
    q()
  }
  # get ARGS
  caco.sampleped <- ARGS[1]
  case.ddb.csv<-ARGS[2]
  ctrl.ddb.csv<-ARGS[3]
  caco.qualvar<-ARGS[4]
  # read files
  ped<-read.table(caco.sampleped,stringsAsFactors=F)
  colnames(ped)<-c("FID","IID","PID","MID","SEX","PHE","SEQTYPE",
                   "EXOMEPREPKIT")
  ca.ddb<-read.csv(case.ddb.csv,stringsAsFactors=F)
  co.ddb<-read.csv(ctrl.ddb.csv, stringsAsFactors=F)
  qv<-read.table(caco.qualvar, stringsAsFactors=F)
  colnames(qv) <- c("GENE","VARID","IID","GT")

  # form one ddb csv
  ddb<-rbind(ca.ddb, co.ddb)

  # change PHE and SEX from 1/2 to 0/1
  ped$SEX <- ped$SEX - 1
  ped$PHE <- ped$PHE - 1

  # add column in ped for prob_Caucasian between 0.5 and 0.99
  ddb.lowconf <- subset(ddb, Caucasian_prob <= 0.99)
  ped$P_CAUCASIAN_LT99 <- ifelse(ped$IID %in% ddb.lowconf$CHGVID,
                                 1,
                                 0)
    
  # assemble p_caucasian classifications
  ddb.lowconf <- subset(ddb, Caucasian_prob <= 0.90)
  ped$P_CAUCASIAN_LT90 <- ifelse(ped$IID %in% ddb.lowconf$CHGVID,1,0)
  ddb.lowconf <- subset(ddb, (Caucasian_prob>0.9) & (Caucasian_prob<=0.99))
  ped$P_CAUCASIAN_90TO99 <- ifelse(ped$IID %in% ddb.lowconf$CHGVID,1,0)
  ddb.lowconf <- subset(ddb, (Caucasian_prob<=0.99))
  ped$P_CAUCASIAN_LT99 <- ifelse(ped$IID %in% ddb.lowconf$CHGVID,1,0)
  ped$P_CAUCASIAN_GE99 <- ifelse(ped$IID %in% ddb.lowconf$CHGVID,0,1)  

  # question 1 : any assoc between binary predictors and 
  # case/control status?
  q1_df <- NULL
  for (pred in c("SEX","P_CAUCASIAN_LT90","P_CAUCASIAN_90TO99",
                 "P_CAUCASIAN_GE99")) {
    q1_df <- Regression("PHE",pred, ped, res.df=q1_df)
  }
  print("Q1 : association between binary predictors and case/control status?")
  print(q1_df)

  # collect samples with at least one qualifying variant in SLITRK5
  qv.slitrk5 <- subset(qv, GENE=="SLITRK5")
  ped$SLITRK5_QV_CARRIER <- ifelse(ped$IID %in% qv.slitrk5$IID, 1, 0)

  # question 2 : anny association between binary traits above and SLITRK5 carrier status?
  q2_df <- NULL
  for (pred in c("SEX","P_CAUCASIAN_LT90","P_CAUCASIAN_90TO99",
                 "P_CAUCASIAN_GE99")) {
    q2_df <- Regression("SLITRK5_QV_CARRIER", pred, ped, res.df=q2_df)
  }
  print("Q2 : assocation between binary traits and SLITRK5 carrier status?")
  print(q2_df)

  # question 3: any assoc btwn traits and SLITRK5 carrier status, partitioned into ca and co?
  ped.ca<-subset(ped, PHE==1)
  ped.co<-subset(ped, PHE==0)
  q3_ca_df <- NULL
    print("Q3 : assocation between binary traits and SLITRK5 carrier status, split into cases and controls?")
  for (pred in c("SEX","P_CAUCASIAN_LT90","P_CAUCASIAN_90TO99",
                 "P_CAUCASIAN_GE99")) {
    q3_ca_df <- Regression("SLITRK5_QV_CARRIER", pred, ped.ca, res.df=q3_ca_df)
  }
  print("CASES")
  print(q3_ca_df)
  q3_co_df <- NULL
  for (pred in c("SEX","P_CAUCASIAN_LT90","P_CAUCASIAN_90TO99",
                 "P_CAUCASIAN_GE99")) {
    q3_co_df <- Regression("SLITRK5_QV_CARRIER", pred, ped.co, res.df=q3_co_df)
  }
  print("CONTROLS")
  print(q3_co_df)

  # question 4 : is there a difference between p_caucasian of cases
  # and of controls?
  print("Q4 : difference between case and control p_caucasian?")
  pcauc <- ddb[,c("CHGVID","Caucasian_prob")]
  colnames(pcauc)[1]<-"IID"
  ped<-merge(ped, pcauc, by="IID")
  co.p_cauc <- subset(ped, PHE==0)
  ca.p_cauc <- subset(ped, PHE==1)
  q3_mdl<-glm(SLITRK5_QV_CARRIER ~ Caucasian_prob, family=binomial,data=ped)
  print(summary(q3_mdl))

  # print SLITRK5 case qv carriers
  print(subset(ped.ca, SLITRK5_QV_CARRIER == 1))

  # print SLITRK5 ctrl qv carriers 
  print(subset(ped.co, SLITRK5_QV_CARRIER == 1))

  q()
 
}

Regression <- function(outcome, pred , df, res.df=NULL) {

  # form mdl strings
  eqn.str.lm <- paste0(pred,"~",outcome)
  eqn.str.lg <- paste0(outcome,"~",pred)

  # linear regression
  mdl.lm <- lm(as.formula(eqn.str.lm), data=df)
  res.lm <- summary(mdl.lm)

  # logistic regression
  mdl.lg <- glm(as.formula(eqn.str.lg),data=df, family=binomial)
  res.lg <- summary(mdl.lg)
  res.df.i <- data.frame(predictor=pred,
                         estimate_linear=res.lm$coefficients[outcome,1],
                         ci95_low_linear=confint(mdl.lm)[outcome,1],
                         ci95_high_linear=confint(mdl.lm)[outcome,2], 
                         p_value_linear=res.lm$coefficients[outcome,4],
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

if (interactive() == F) {
  Main()
}
