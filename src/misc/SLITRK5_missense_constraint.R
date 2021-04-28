
library(abind)
library(ggplot2)

## PARAM
NCA <- list(
            "Cluster 0"=599,
            "Cluster 1"=11,
            "Cluster 2"=22,
            "Cluster 3"=205,
            "Cluster 4"=139,
            "Cluster 5"=111,
            "Cluster 6"=51,
            "Cluster 7"=26,
            "Cluster 8"=14,
            "Cluster 9"=8,
            "Cluster 10"=77
           )
NCO <- list(
            "Cluster 0"=2449,
            "Cluster 1"=2139,
            "Cluster 2"=1947,
            "Cluster 3"=1168,
            "Cluster 4"=1113,
            "Cluster 5"=953,
            "Cluster 6"=735,
            "Cluster 7"=313,
            "Cluster 8"=282,
            "Cluster 9"=276,
            "Cluster 10"=205 
           )
            

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 5) {
    cat("SLITRK5_missense_constraint.R ",
        "<in_qv.csv> <lva.csv> <in_qv.ann.csv> ",
        "<results.cmh.csv> <results.nvar.csv>\n")
    q()
  }

  # get ARGS
  in_qv.csv <- ARGS[1]
  lva.csv <- ARGS[2]
  in_qv.ann.csv <- ARGS[3]
  cmh.csv <- ARGS[4]
  nvar.csv <- ARGS[5]

  # read qv csv, lva csv
  qv <- read.csv(in_qv.csv, stringsAsFactors=F)
  ann <- read.csv(lva.csv, stringsAsFactors=F)

  # keep subset of ann that has missense annots of interest
  ann <- ann[,c("Variant.ID","MPC","MTR","CCR.Percentile")]

  # merge tables 
  qv_ann <- merge(qv, ann, by='Variant.ID', all.x=T)

  # write merged table to file
  write.csv(qv_ann,
            file=in_qv.ann.csv,
            row.names=F)

  # subset of table that are missense only
  qv_ann <- subset(qv_ann,
                   (grepl("missense", Effect)) &
                   (Variant.Type=="snv")
                  )
  qv_ann$PHE <- ifelse(qv_ann$Sample.Phenotype=="case", 1, 0)

  # baseline odds ratio (ppn2 >= 0.957)
  cmh_res<-cmh_test(qv_ann, "ppn2_damaging", cmh_res=NULL)
  nvar_ca_baseline <- nrow(subset(qv_ann, PHE==1))
  nvar_co_baseline <- nrow(subset(qv_ann, PHE==0))

  # get thresholds
  mpc_threshold <- summary(qv_ann$MPC)[['3rd Qu.']]
  mtr_threshold <- summary(qv_ann$MTR)[['1st Qu.']]
  ccr_threshold <- summary(qv_ann$CCR.Percentile)[['3rd Qu.']]
  cat("Quartile-based threshold for MPC constraint is >=",mpc_threshold,"\n")
  cat("Quartile-based threshold for MTR constraint is <=",mtr_threshold,"\n") 
  cat("Quartile-based threshold for CCR constraint is >=",ccr_threshold,"\n")

  # constraint comparisons (threshold set as upper quartile derived
  # across these particular variants identified in slitrk5 in 
  # cases and controls)
                        
  # constraint 1 : MPC
  qv_ann.mpc <- subset(qv_ann, MPC>=mpc_threshold)
  cmh_res<-cmh_test(qv_ann.mpc, paste0("ppn2_damaging, MPC>=",mpc_threshold),cmh_res=cmh_res)
  nvar_res<-nvar_fisher(qv_ann.mpc, nvar_ca_baseline, nvar_co_baseline, 
                        paste0("ppn2_damaging, MPC>=",mpc_threshold), nvar_res=NULL)

  # constraint 2 : MTR 
  qv_ann.mtr <- subset(qv_ann, MTR<=mtr_threshold)
  cmh_res<-cmh_test(qv_ann.mtr, paste0("ppn2_damaging, MTR<=",mtr_threshold), cmh_res=cmh_res)
  nvar_res<-nvar_fisher(qv_ann.mtr, nvar_ca_baseline, nvar_co_baseline,
                        paste0("ppn2_damaging, MTR<=",mtr_threshold), nvar_res=nvar_res)

  # constraint 3 : CCR 
  qv_ann.ccr <- subset(qv_ann, CCR.Percentile>=ccr_threshold)
  cmh_res<-cmh_test(qv_ann.ccr, paste0("ppn2_damaging, CCR>=",ccr_threshold), cmh_res=cmh_res)
  nvar_res<-nvar_fisher(qv_ann.ccr, nvar_ca_baseline, nvar_co_baseline,
                        paste0("ppn2_damaging, CCR>=",ccr_threshold), nvar_res=nvar_res)
 
  # write results  to csv
  write.csv(cmh_res, file=cmh.csv, row.names=F)
  write.csv(nvar_res, file=nvar.csv, row.names=F)
}

cmh_test <- function(qv_ann, constraint_str, cmh_res=NULL) {
  tbl_3d <- NULL
  qv_ann_ca <- subset(qv_ann, PHE==1)
  qv_ann_co <- subset(qv_ann, PHE==0)
  nca_tot_1 <- 0
  nca_tot_0 <- 0
  nco_tot_1 <- 0
  nco_tot_0 <- 0
  for (cluster in names(NCA)) {
    nca_1 <- nrow(subset(qv_ann_ca, Cluster==cluster))
    nco_1 <- nrow(subset(qv_ann_co, Cluster==cluster))
    nca_0 <- NCA[[cluster]] - nca_1
    nco_0 <- NCO[[cluster]] - nco_1
    tbl <- data.frame(ca=c(nca_1, nca_0), 
                      co=c(nco_1, nco_0)
                     )
    nca_tot_1 <- nca_tot_1 + nca_1
    nca_tot_0 <- nca_tot_0 + nca_0
    nco_tot_1 <- nco_tot_1 + nco_1
    nco_tot_0 <- nco_tot_0 + nco_0
    if (is.null(tbl_3d)) {
      tbl_3d <- tbl
    } else {
      tbl_3d <- abind(tbl_3d, tbl, along=3)
    }
  }
  res <- mantelhaen.test(tbl_3d, exact=T)
  if (is.null(cmh_res)) {
    cmh_res <- data.frame(constraint=character(),
                          nca_1=numeric(),
                          nca_0=numeric(),
                          nco_1=numeric(),
                          nco_0=numeric(),
                          perc_ca_1=numeric(),
                          perc_co_1=numeric(),
                          cmh_or_95ci_l=numeric(),
                          cmh_or_95ci_u=numeric(),
                          cmh_or=numeric(),
                          cmh_p=numeric()
                         )
  }
  cmh_res <- rbind(cmh_res,
                   data.frame(
                              constraint=constraint_str,
                              nca_1=nca_tot_1,
                              nca_0=nca_tot_0,
                              nco_1=nco_tot_1,
                              nco_0=nco_tot_0,
                              perc_ca_1=nca_tot_1/(nca_tot_1+nca_tot_0),
                              perc_co_1=nco_tot_1/(nco_tot_1+nco_tot_0),
                              cmh_or_95ci_l=res$conf.int[1],
                              cmh_or_95ci_u=res$conf.int[2],
                              cmh_or=res$estimate,
                              cmh_p=res$p.value
                             )
                  ) 
                  
  return(cmh_res)
}

nvar_fisher <- function(qv_ann, nvar_ca_baseline, nvar_co_baseline, 
                        constraint_str, nvar_res=NULL) {
  nvar_ca_1 <- nrow(subset(qv_ann, PHE==1))
  nvar_co_1 <- nrow(subset(qv_ann, PHE==0))
  nvar_ca_0 <- nvar_ca_baseline - nvar_ca_1
  nvar_co_0 <- nvar_co_baseline - nvar_co_1
  df <- data.frame(ca=c(nvar_ca_1, nvar_ca_0),
                   co=c(nvar_co_1, nvar_co_0)
                  )
  res <- fisher.test(df, alternative='greater')
  if (is.null(nvar_res)) {
    nvar_res <- data.frame(constraint=character(),
                           nvar_ca_1=numeric(), nvar_ca_0=numeric(),
                           nvar_ca_1=numeric(), nvar_ca_0=numeric(),
                           perc_ca_1=numeric(), perc_co_1=numeric(),
                           fet_or_95ci_l=numeric(), fet_or_95ci_u=numeric(),
                           fet_or=numeric(), fet_p=numeric())
  }
  nvar_res <- rbind(nvar_res,
                    data.frame(constraint=constraint_str,
                               nvar_ca_1=nvar_ca_1, nvar_ca_0=nvar_ca_0,
                               nvar_co_1=nvar_co_1, nvar_co_0=nvar_co_0,
                               perc_ca_1=nvar_ca_1/(nvar_ca_1+nvar_ca_0),
                               perc_co_1=nvar_co_1/(nvar_co_1+nvar_co_0),
                               fet_or_95ci_l=res$conf.int[1], 
                               fet_or_95ci_u=res$conf.int[2],
                               fet_or=res$estimate, 
                               fet_p=res$p.value)
                   )

  return(nvar_res)
}

if (interactive() == F) {
  Main()
}
