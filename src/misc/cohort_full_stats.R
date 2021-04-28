

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 18) {
    cat("cohort_full_stats.R <cases.ddb.csv> <ctrls.ddb.csv>",
        "<prepids.keep.list.txt>",
        "<clps.c0.sampleped> <clps.c1.sampleped>",
        "<clps.c2.sampleped> <clps.c3.sampleped>",
        "<clps.c4.sampleped> <clps.c5.sampleped>",
        "<clps.c6.sampleped> <clps.c7.sampleped>",
        "<clps.c8.sampleped> <clps.c9.sampleped>",
        "<clps.c10.sampleped>",
        "<lofrate.European.sampleped> <dnm.trios.sampleped>",
        "<dnm.quartets.sampleped> <outroot>")
    q()
  }

  # get ARGS
  cases.ddb.csv <- ARGS[1]
  ctrls.ddb.csv <- ARGS[2]
  prepids.keep.list.txt <- ARGS[3]
  clps.c0.sampleped <- ARGS[4]
  clps.c1.sampleped <- ARGS[5]
  clps.c2.sampleped <- ARGS[6]
  clps.c3.sampleped <- ARGS[7]
  clps.c4.sampleped <- ARGS[8]
  clps.c5.sampleped <- ARGS[9]
  clps.c6.sampleped <- ARGS[10]
  clps.c7.sampleped <- ARGS[11]
  clps.c8.sampleped <- ARGS[12]
  clps.c9.sampleped <- ARGS[13]
  clps.c10.sampleped <- ARGS[14]
  lofrate.european.sampleped <- ARGS[15]
  dnm.trios.sampleped <- ARGS[16]
  dnm.quartets.sampleped <- ARGS[17]
  outroot <- ARGS[18]

  # read dragendb csvs
  ca.ddb <- read.csv(cases.ddb.csv, stringsAsFactors=F)
  co.ddb <- read.csv(ctrls.ddb.csv, stringsAsFactors=F)

  # read list of prep ids to keep, only keep ddb rows if they
  # have a prepid that is within this file
  prepids <- scan(prepids.keep.list.txt, what=character(), quiet=T)
  ca.ddb <- subset(ca.ddb, prepID %in% prepids)
  co.ddb <- subset(co.ddb, prepID %in% prepids)

  # make sure control ddb doesn't have ocd in it
  co.ddb <- subset(co.ddb, 
                   (grepl("ocd",CHGVID)==F) &
                   (BroadPhenotype!="obsessive compulsive disorder")
                  )

  # general filtration on both ddb dataframes
  ca.ddb <- subset(ca.ddb, SeqType %in% c("Exome","Genome"))
  co.ddb <- subset(co.ddb, SeqType %in% c("Exome","Genome"))

  # read samplepeds
  clps.c0 <- ReadSampleped(clps.c0.sampleped)
  clps.c1 <- ReadSampleped(clps.c1.sampleped)
  clps.c2 <- ReadSampleped(clps.c2.sampleped)
  clps.c3 <- ReadSampleped(clps.c3.sampleped)
  clps.c4 <- ReadSampleped(clps.c4.sampleped)
  clps.c5 <- ReadSampleped(clps.c5.sampleped)
  clps.c6 <- ReadSampleped(clps.c6.sampleped)
  clps.c7 <- ReadSampleped(clps.c7.sampleped)
  clps.c8 <- ReadSampleped(clps.c8.sampleped)
  clps.c9 <- ReadSampleped(clps.c9.sampleped)
  clps.c10 <- ReadSampleped(clps.c10.sampleped)
  lofrate.european <- ReadSampleped(lofrate.european.sampleped)
  dnm.trios <- ReadSampleped(dnm.trios.sampleped)
  dnm.quartets <- ReadSampleped(dnm.quartets.sampleped)

  # make a merged sampleped for collapsing
  clps.list <- list("Cluster0"=clps.c0,
                    "Cluster1"=clps.c1,
                    "Cluster2"=clps.c2,
                    "Cluster3"=clps.c3,
                    "Cluster4"=clps.c4,
                    "Cluster5"=clps.c5,
                    "Cluster6"=clps.c6,
                    "Cluster7"=clps.c7,
                    "Cluster8"=clps.c8,
                    "Cluster9"=clps.c9,
                    "Cluster10"=clps.c10)
  clps <- rbind(clps.c0, clps.c1, clps.c2, clps.c3, clps.c4,
                clps.c5, clps.c6, clps.c7, clps.c8, clps.c9,
                clps.c10)
  clps.ca <- subset(clps, clps[,6]==2)
  clps.co <- subset(clps, clps[,6]==1)

  # apply duplicate sample filter to control ddb dataframe, 
  # as some samples were sequenced twice. We used the samples with the 
  # largest (ie. most recent) prepid
  # ca.ddb <- DuplicateSampleFilter(ca.ddb)
  # co.ddb <- DuplicateSampleFilter(co.ddb)

  clps.list[["TOTAL"]] <- clps
 
  # make sure samples are OCD, not nyspi
  ca.ddb <- subset(ca.ddb, BroadPhenotype=="obsessive compulsive disorder")
  ca.ddb <- subset(ca.ddb, SeqType %in% c("Exome","Genome"))
  ca.ddb <- subset(ca.ddb, grepl("nyspiocd", ca.ddb$CHGVID) == F)

  # subset on kits included in collapsing or lof rate (or WGS)
  kits <- c(clps$EXOMEPREPKIT, lofrate.european$EXOMEPREPKIT)
  kits <- unique(sort(kits))
  kits <- c(kits, "N/A")
  seqtypes <- c("Exome","Genome")
  ca.ddb <- subset(ca.ddb, (exomeKit %in% kits) & (SeqType %in% seqtypes))
  co.ddb <- subset(co.ddb, (exomeKit %in% kits) & (SeqType %in% seqtypes))

  # print stats (total number of cases)
  cat("Total number of processed cases :",nrow(ca.ddb), "\n")

  # combine case and control ddb, form a df with iids only in collapsing,
  # another df with iids only in clps rate test
  ddb <- rbind(ca.ddb, co.ddb)
  ddb.clps <- subset(ddb, CHGVID %in% clps$IID)
  ddb.lofrate <- subset(ddb, CHGVID %in% lofrate.european$IID)
  iids <-ddb.lofrate$CHGVID
  iids <- sort(as.character(iids))

  # print total number of samples in collapsing
  cat("Total number of cases in collapsing meta :", 
      nrow(subset(clps, PHE==2)), "\n")
  cat("Total number of controls in collapsing meta :",
      nrow(subset(clps, PHE==1)), "\n")

  # print total number of samples in lof rate comparison
  cat("Total number of cases in LoF rate comparison :",
      nrow( subset( ddb.lofrate, grepl("ocd", CHGVID) ) ), "\n")
  cat("Total number of controls in LoF rate comparison :",
      nrow( subset( ddb.lofrate, grepl("ocd", CHGVID)==F) ), "\n")

  # print total number of trios
  cat("Total number of trios analyzed (DNM analysis) :", 
      length(unique(dnm.trios$FID)), "\n")
  cat("Total number of quartets analyzed (shared DNM analysis) :",
      length(unique(dnm.quartets$FID)), "\n")

  # get all broad phenotype, to use in forming pheno counts df
  broadphenos <- sort(unique(ddb.clps$BroadPhenotype))
  broadphenos <- broadphenos[broadphenos != "obsessive compulsive disorder"]
  broadpheno_counts <- data.frame(broad_phenotype=broadphenos)
  rownames(broadpheno_counts) <- broadpheno_counts$broad_phenotype

  # establish count df for keeping track of nca, nco per group
  cohort_df <- data.frame(analysis=character(), cohort=character(),
                          n_cases=numeric(), n_controls=numeric(),
                          n_cases_trio=numeric(), n_cases_quartet=numeric(),
                          n_cases_singleton=numeric(),
                          n_controls_healthy=numeric(),
                          percent_controls_healthy=numeric(),
                          cases_ccdsbasescov10x_mean=numeric(),
                          controls_ccdsbasescov10x_mean=numeric(),
                          cases_ccdsbasescov10x_median=numeric(),
                          controls_ccdsbasescov10x_median=numeric(),
                          n_cases_exome=numeric(),
                          n_controls_exome=numeric(),
                          percent_cases_exome=numeric(),
                          percent_controls_exome=numeric(),
                          n_cases_genome=numeric(),
                          n_controls_genome=numeric(),
                          percent_cases_genome=numeric(),
                          percent_controls_genome=numeric(),
                          n_cases_exome_Roche=numeric(),
                          n_controls_exome_Roche=numeric(),
                          percent_cases_exome_Roche=numeric(),
                          percent_controls_exome_Roche=numeric(),
                          n_cases_exome_IDTERPv1=numeric(),
                          n_controls_exome_IDTERPv1=numeric(),
                          percent_cases_exome_IDTERPv1=numeric(),
                          percent_controls_exome_IDTERPv1=numeric(),
                          n_cases_exome_65MB=numeric(),
                          n_controls_exome_65MB=numeric(),
                          percent_cases_exome_65MB=numeric(),
                          percent_controls_exome_65MB=numeric(),
                          n_cases_exome_AgilentV5=numeric(),
                          n_controls_exome_AgilentV5=numeric(),
                          percent_cases_exome_AgilentV5=numeric(),
                          percent_controls_exome_AgilentV5=numeric(),
                          n_cases_exome_RocheV2=numeric(),
                          n_controls_exome_RocheV2=numeric(),
                          percent_cases_exome_RocheV2=numeric(),
                          percent_controls_exome_RocheV2=numeric(),
                          n_cases_exome_RocheV1=numeric(),
                          n_controls_exome_RocheV1=numeric(),
                          percent_cases_exome_RocheV1=numeric(),
                          percent_controls_exome_RocheV1=numeric()
                         )

  # get cohort df table values per cohort (gene-based collapsing)
 cohorts <- c("Cluster0",
               "Cluster1",
               "Cluster2",
               "Cluster3",
               "Cluster4",
               "Cluster5",
               "Cluster6",
               "Cluster7",
               "Cluster8",
               "Cluster9",
               "Cluster10",
               "TOTAL")
  healthy_phe <- c("control","healthy family member")
  for (cohort in cohorts) {
    ca.iids.c <- subset(clps.list[[cohort]], PHE==2)$IID
    co.iids.c <- subset(clps.list[[cohort]], PHE==1)$IID
    ddb.ca.c <- subset(ddb.clps, CHGVID %in% ca.iids.c)
    ddb.co.c <- subset(ddb.clps, CHGVID %in% co.iids.c)
    ddb.co.c.h <- subset(ddb.co.c, BroadPhenotype %in% healthy_phe)
    nca.c <- length(ca.iids.c)
    nco.c <- length(co.iids.c)

    # do stats on n and perc of controls by broadphenotype
    n_col <- paste0("n_",cohort, "_collapsing")
    perc_col <- paste0("perc_",cohort, "_collapsing")
    broadpheno_counts[[n_col]] <- rep(0, nrow(broadpheno_counts))
    broadpheno_counts[[perc_col]] <- rep(0, nrow(broadpheno_counts))
    for (broadpheno in rownames(broadpheno_counts)) {
      n.b.i <- nrow(subset(ddb.co.c, BroadPhenotype == broadpheno))
      perc.b.i <- n.b.i/nco.c
      broadpheno_counts[broadpheno, n_col] <- n.b.i
      broadpheno_counts[broadpheno, perc_col] <- perc.b.i 
    }

    # do stats on CCDSBasesCov10X
    ca.mean.c <- mean(ddb.ca.c$CCDSBasesCov10X[is.na(ddb.ca.c$CCDSBasesCov10X)==F])
    co.mean.c <- mean(ddb.co.c$CCDSBasesCov10X[is.na(ddb.co.c$CCDSBasesCov10X)==F])
    ca.median.c <- median(ddb.ca.c$CCDSBasesCov10X[is.na(ddb.ca.c$CCDSBasesCov10X)==F])
    co.median.c <- median(ddb.co.c$CCDSBasesCov10X[is.na(ddb.co.c$CCDSBasesCov10X)==F])

    # get number of cases/controls by kit/WGS
    ca.kitcounts.c <- KitCounts(subset(clps.list[[cohort]], PHE==2))
    co.kitcounts.c <- KitCounts(subset(clps.list[[cohort]], PHE==1))
    
    # get number of healthy/nonhealthy listed controls
    nco.healthy.c <- nrow(subset(clps.list[[cohort]],
                                 (IID %in% ddb.co.c.h$CHGVID) &
                                 (PHE==1)
                                )
                         )
    nco.nonhealthy.c <- nco.c - nco.healthy.c

    # get breakdown of cases for cohort (ie. singleton, trio, quartet)
    fam.iids <- c(dnm.trios$IID, dnm.quartets$IID)
    nca.c.trio <- length(ca.iids.c[ca.iids.c %in% dnm.trios$IID])
    nca.c.quartet <- length(ca.iids.c[ca.iids.c %in% dnm.quartets$IID])
    nca.c.singleton <- length(ca.iids.c[(ca.iids.c %in% fam.iids)==F])

    # append cluster-level counts to dataframe
    cohort_df <- rbind(cohort_df,
                       data.frame(analysis="gene-based collapsing",
                                  cohort=cohort, 
                                  n_cases=nca.c,
                                  n_controls=nco.c, 
                                  n_cases_trio=nca.c.trio,
                                  n_cases_quartet=nca.c.quartet,
                                  n_cases_singleton=nca.c.singleton,
                                  n_controls_healthy=nco.healthy.c,
                                  percent_controls_healthy=nco.healthy.c/nco.c,
                                  cases_ccdsbasescov10x_mean=ca.mean.c,
                                  controls_ccdsbasescov10x_mean=co.mean.c,
                                  cases_ccdsbasescov10x_median=ca.median.c,
                                  controls_ccdsbasescov10x_median=co.median.c,
                                  n_cases_exome=ca.kitcounts.c[["Exome"]],
                                  n_controls_exome=co.kitcounts.c[["Exome"]],
                                  percent_cases_exome=ca.kitcounts.c[["Exome"]]/nca.c,
                                  percent_controls_exome=co.kitcounts.c[["Exome"]]/nco.c,
                                  n_cases_genome=ca.kitcounts.c[["Genome"]],
                                  n_controls_genome=co.kitcounts.c[["Genome"]],
                                  percent_cases_genome=ca.kitcounts.c[["Genome"]]/nca.c,
                                  percent_controls_genome=co.kitcounts.c[["Genome"]]/nco.c,
                                  n_cases_exome_Roche=ca.kitcounts.c[["Roche"]],
                                  n_controls_exome_Roche=co.kitcounts.c[["Roche"]],
                                  percent_cases_exome_Roche=ca.kitcounts.c[["Roche"]]/nca.c,
                                  percent_controls_exome_Roche=co.kitcounts.c[["Roche"]]/nco.c,
                                  n_cases_exome_IDTERPv1=ca.kitcounts.c[["IDTERPv1"]],
                                  n_controls_exome_IDTERPv1=co.kitcounts.c[["IDTERPv1"]],
                                  percent_cases_exome_IDTERPv1=ca.kitcounts.c[["IDTERPv1"]]/nca.c,
                                  percent_controls_exome_IDTERPv1=co.kitcounts.c[["IDTERPv1"]]/nco.c,
                                  n_cases_exome_65MB=ca.kitcounts.c[["65MB"]],
                                  n_controls_exome_65MB=co.kitcounts.c[["65MB"]],
                                  percent_cases_exome_65MB=ca.kitcounts.c[["65MB"]]/nca.c,
                                  percent_controls_exome_65MB=co.kitcounts.c[["65MB"]]/nco.c,
                                  n_cases_exome_AgilentV5=ca.kitcounts.c[["AgilentV5"]],
                                  n_controls_exome_AgilentV5=co.kitcounts.c[["AgilentV5"]],
                                  percent_cases_exome_AgilentV5=ca.kitcounts.c[["AgilentV5"]]/nca.c,
                                  percent_controls_exome_AgilentV5=co.kitcounts.c[["AgilentV5"]]/nco.c,
                                  n_cases_exome_RocheV2=0,
                                  n_controls_exome_RocheV2=0,
                                  percent_cases_exome_RocheV2=0,
                                  percent_controls_exome_RocheV2=0,
                                  n_cases_exome_RocheV1=0,
                                  n_controls_exome_RocheV1=0,
                                  percent_cases_exome_RocheV1=0,
                                  percent_controls_exome_RocheV1=0
                                 )
                              )

  }

  # add stats for lof rate test
  ca.iids.r <- subset(lofrate.european, PHE==2)$IID
  co.iids.r <- subset(lofrate.european, PHE==1)$IID
  nca.r <- length(unique(ca.iids.r))
  nco.r <- length(unique(co.iids.r))
  ddb.ca.r <- subset(ddb.lofrate, CHGVID %in% ca.iids.r)
  ddb.co.r <- subset(ddb.lofrate, CHGVID %in% co.iids.r)
  ddb.co.r.h <- subset(ddb.co.r, BroadPhenotype %in% healthy_phe)

  # do stats on CCDSBasesCov10X
  ca.mean.r <- mean(ddb.ca.r$CCDSBasesCov10X[is.na(ddb.ca.r$CCDSBasesCov10X)==F])
  co.mean.r <- mean(ddb.co.r$CCDSBasesCov10X[is.na(ddb.co.r$CCDSBasesCov10X)==F])
  ca.median.r <- median(ddb.ca.r$CCDSBasesCov10X[is.na(ddb.ca.r$CCDSBasesCov10X)==F])
  co.median.r <- median(ddb.co.r$CCDSBasesCov10X[is.na(ddb.co.r$CCDSBasesCov10X)==F])

  # get number of cases/controls by kit/WGS
  ca.kitcounts.r <- KitCounts(subset(lofrate.european, PHE==2))
  co.kitcounts.r <- KitCounts(subset(lofrate.european, PHE==1))

  # get breakdown of cases for cohort (ie. singleton, trio, quartet)
  fam.iids <- c(dnm.trios$IID, dnm.quartets$IID)
  nca.r.trio <- length(ca.iids.r[ca.iids.r %in% dnm.trios$IID])
  nca.r.quartet <- length(ca.iids.r[ca.iids.r %in% dnm.quartets$IID])
  nca.r.singleton <- length(ca.iids.r[(ca.iids.r %in% fam.iids)==F])

  nco.healthy.r <- nrow(subset(lofrate.european,
                               (IID %in% ddb.co.r.h$CHGVID) &
                               (PHE==1)
                              )
                       )
  nco.nonhealthy.r <- nco.r - nco.healthy.r
  cohort_df <- rbind(cohort_df,
                     data.frame(analysis="LoF rate",
                                cohort="European",
                                n_cases=nca.r,
                                n_controls=nco.r,
                                n_cases_trio=nca.r.trio,
                                n_cases_quartet=nca.r.quartet,
                                n_cases_singleton=nca.r.singleton,
                                n_controls_healthy=nco.healthy.r,
                                percent_controls_healthy=nco.healthy.r/nco.r,
                                cases_ccdsbasescov10x_mean=ca.mean.r,
                                controls_ccdsbasescov10x_mean=co.mean.r,
                                cases_ccdsbasescov10x_median=ca.median.r,
                                controls_ccdsbasescov10x_median=co.median.r,
                                n_cases_exome=ca.kitcounts.r[["Exome"]],
                                n_controls_exome=co.kitcounts.r[["Exome"]],
                                percent_cases_exome=ca.kitcounts.r[["Exome"]]/nca.r,
                                percent_controls_exome=co.kitcounts.r[["Exome"]]/nco.r,
                                n_cases_genome=ca.kitcounts.r[["Genome"]],
                                n_controls_genome=co.kitcounts.r[["Genome"]],
                                percent_cases_genome=ca.kitcounts.r[["Genome"]]/nca.r,
                                percent_controls_genome=co.kitcounts.r[["Genome"]]/nco.r,
                                n_cases_exome_Roche=ca.kitcounts.r[["Roche"]],
                                n_controls_exome_Roche=co.kitcounts.r[["Roche"]],
                                percent_cases_exome_Roche=ca.kitcounts.r[["Roche"]]/nca.r,
                                percent_controls_exome_Roche=co.kitcounts.r[["Roche"]]/nco.r,
                                n_cases_exome_IDTERPv1=ca.kitcounts.r[["IDTERPv1"]],
                                n_controls_exome_IDTERPv1=co.kitcounts.r[["IDTERPv1"]],
                                percent_cases_exome_IDTERPv1=ca.kitcounts.r[["IDTERPv1"]]/nca.r,
                                percent_controls_exome_IDTERPv1=co.kitcounts.r[["IDTERPv1"]]/nco.r,
                                n_cases_exome_65MB=ca.kitcounts.r[["65MB"]],
                                n_controls_exome_65MB=co.kitcounts.r[["65MB"]],
                                percent_cases_exome_65MB=ca.kitcounts.r[["65MB"]]/nca.r,
                                percent_controls_exome_65MB=co.kitcounts.r[["65MB"]]/nco.r,
                                n_cases_exome_AgilentV5=ca.kitcounts.r[["AgilentV5"]],
                                n_controls_exome_AgilentV5=co.kitcounts.r[["AgilentV5"]],
                                percent_cases_exome_AgilentV5=ca.kitcounts.r[["AgilentV5"]]/nca.r,
                                percent_controls_exome_AgilentV5=co.kitcounts.r[["AgilentV5"]]/nco.r,
                                n_cases_exome_RocheV2=0,
                                n_controls_exome_RocheV2=0,
                                percent_cases_exome_RocheV2=0,
                                percent_controls_exome_RocheV2=0,
                                n_cases_exome_RocheV1=0,
                                n_controls_exome_RocheV1=0,
                                percent_cases_exome_RocheV1=0,
                                percent_controls_exome_RocheV1=0
                             
                            )
  
                    )

  # add stats for control count per broad pheno in LoF rate cohort
  n_col <- "n_European_LoFrate"
  perc_col <- "perc_European_LoFrate"
  broadpheno_counts[[n_col]] <- rep(0, nrow(broadpheno_counts))
  broadpheno_counts[[perc_col]] <- rep(0, nrow(broadpheno_counts))
  for (broadpheno in rownames(broadpheno_counts)) {
    n_i <- nrow(subset(ddb.lofrate, BroadPhenotype==broadpheno))
    perc_i <- n_i / nco.r
    broadpheno_counts[broadpheno, n_col] <- n_i
    broadpheno_counts[broadpheno, perc_col] <- perc_i
  }

  # add row for trio analysis
  ca.iids.t <- subset(dnm.trios, PHE==2)$IID
  nca.t <- length(unique(ca.iids.t))
  nco.t <- 1911
  ddb.ca.t <- subset(ca.ddb, CHGVID %in% ca.iids.t)
  ca.mean.t <- mean(ddb.ca.t$CCDSBasesCov10X)
  ca.median.t <- median(ddb.ca.t$CCDSBasesCov10X)

  # get number of cases/controls by kit/WGS
  ca.kitcounts.t <- KitCounts(subset(dnm.trios, PHE==2))
  cohort_df <- rbind(cohort_df,
                     data.frame(analysis="de novo mutations",
                                cohort="Trios",
                                n_cases=nca.t,
                                n_controls=nco.t,
                                n_cases_trio=nca.t,
                                n_cases_quartet=0,
                                n_cases_singleton=0,
                                n_controls_healthy=nco.t,
                                percent_controls_healthy=1,
                                cases_ccdsbasescov10x_mean=ca.mean.t,
                                controls_ccdsbasescov10x_mean=NA,
                                cases_ccdsbasescov10x_median=ca.median.t,
                                controls_ccdsbasescov10x_median=NA,
                                n_cases_exome=587,
                                n_controls_exome=nco.t,
                                percent_cases_exome=1,
                                percent_controls_exome=1,
                                n_cases_genome=0,
                                n_controls_genome=0,
                                percent_cases_genome=0,
                                percent_controls_genome=0,
                                n_cases_exome_Roche=ca.kitcounts.t[["Roche"]],
                                n_controls_exome_Roche=0,
                                percent_cases_exome_Roche=ca.kitcounts.t[["Roche"]]/nca.t,
                                percent_controls_exome_Roche=0,
                                n_cases_exome_IDTERPv1=ca.kitcounts.t[["IDTERPv1"]],
                                n_controls_exome_IDTERPv1=0,
                                percent_cases_exome_IDTERPv1=ca.kitcounts.t[["IDTERPv1"]]/nca.t,
                                percent_controls_exome_IDTERPv1=0,
                                n_cases_exome_65MB=0,
                                n_controls_exome_65MB=0,
                                percent_cases_exome_65MB=0,
                                percent_controls_exome_65MB=0,
                                n_cases_exome_AgilentV5=0,
                                n_controls_exome_AgilentV5=0,
                                percent_cases_exome_AgilentV5=0,
                                percent_controls_exome_AgilentV5=0,
                                n_cases_exome_RocheV2=0,
                                n_controls_exome_RocheV2=nco.t-19,
                                percent_cases_exome_RocheV2=0,
                                percent_controls_exome_RocheV2=(nco.t-19)/nco.t,
                                n_cases_exome_RocheV1=0,
                                n_controls_exome_RocheV1=19,
                                percent_cases_exome_RocheV1=0,
                                percent_controls_exome_RocheV1=19/nco.t
 
                               )
                             )

  # add pheno stats for trio controls
  n_col <- "n_trioprobands_DNM"
  perc_col <- "perc_trioprobands_DNM"
  broadpheno_counts[[n_col]] <- rep(0, nrow(broadpheno_counts))
  broadpheno_counts[[perc_col]] <- rep(0, nrow(broadpheno_counts))
  broadpheno_counts["healthy family member", n_col] <- 1911
  broadpheno_counts["healthy family member", perc_col] <- 1

  # add row for quartet analysis
  ca.iids.q <- subset(dnm.quartets, PHE==2)$IID
  nca.q <- length(unique(ca.iids.q))
  nco.q <- 1911

  # get number of cases/controls by kit/WGS
  ca.kitcounts.q <- KitCounts(subset(dnm.quartets, PHE==2))
  ddb.ca.q <- subset(ca.ddb, CHGVID %in% ca.iids.q)
  ca.mean.q <- mean(ddb.ca.q$CCDSBasesCov10X)
  ca.median.q <- median(ddb.ca.q$CCDSBasesCov10X)
  cohort_df <- rbind(cohort_df,
                     data.frame(analysis="de novo mutations",
                                cohort="Quartets",
                                n_cases=nca.q,
                                n_controls=nco.q,
                                n_cases_trio=nca.q,
                                n_cases_quartet=0,
                                n_cases_singleton=0,
                                n_controls_healthy=nco.q,
                                percent_controls_healthy=1,
                                cases_ccdsbasescov10x_mean=ca.mean.q,
                                controls_ccdsbasescov10x_mean=NA,
                                cases_ccdsbasescov10x_median=ca.median.q,
                                controls_ccdsbasescov10x_median=NA,
                                n_cases_exome=nca.q,
                                n_controls_exome=nco.q,
                                percent_cases_exome=1,
                                percent_controls_exome=1,
                                n_cases_genome=0,
                                n_controls_genome=0,
                                percent_cases_genome=0,
                                percent_controls_genome=0,
                                n_cases_exome_Roche=ca.kitcounts.q[["Roche"]],
                                n_controls_exome_Roche=0,
                                percent_cases_exome_Roche=ca.kitcounts.q[["Roche"]]/nca.q,
                                percent_controls_exome_Roche=0,
                                n_cases_exome_IDTERPv1=ca.kitcounts.q[["IDTERPv1"]],
                                n_controls_exome_IDTERPv1=0,
                                percent_cases_exome_IDTERPv1=ca.kitcounts.q[["IDTERPv1"]]/nca.q,
                                percent_controls_exome_IDTERPv1=0,
                                n_cases_exome_65MB=0,
                                n_controls_exome_65MB=0,
                                percent_cases_exome_65MB=0,
                                percent_controls_exome_65MB=0,
                                n_cases_exome_AgilentV5=0,
                                n_controls_exome_AgilentV5=0,
                                percent_cases_exome_AgilentV5=0,
                                percent_controls_exome_AgilentV5=0,
                                n_cases_exome_RocheV2=0,
                                n_controls_exome_RocheV2=nco.q-19,
                                percent_cases_exome_RocheV2=0,
                                percent_controls_exome_RocheV2=(nco.q-19)/nco.q,
                                n_cases_exome_RocheV1=0,
                                n_controls_exome_RocheV1=19,
                                percent_cases_exome_RocheV1=0,
                                percent_controls_exome_RocheV1=19/nco.q

                               )
                             )

  # adjust column names, write to csv
  colnames(cohort_df) <- c("Analysis", "Cohort", "N_Cases", "N_Controls",
                           "N_Cases_Trio","N_Cases_Quartet", "N_Cases_Singleton",
                           "N_Controls_Healthy","Percent_Controls_Healthy",
                           "Cases_CCDSBasesCovered10X_Mean","Controls_CCDSBasesCovered10X_Mean",
                           "Cases_CCDSBasesCovered10X_Median","Controls_CCDSBasesCovered10X_Median",
                           "N_Cases_Exome","N_Controls_Exome", "Percent_Cases_Exome",
                           "Percent_Controls_Exome", "N_Cases_Genome",
                           "N_Controls_Genome", "Percent_Cases_Genome",
                           "Percent_Controls_Genome",
                           "N_Cases_Exome_RocheV3","N_Controls_Exome_RocheV3",
                           "Percent_Cases_Exome_RocheV3","Percent_Controls_Exome_RocheV3",
                           "N_Cases_Exome_IDTERPv1","N_Controls_Exome_IDTERPv1",
                           "Percent_Cases_Exome_IDTERPv1","Percent_Controls_Exome_IDTERPv1",
                           "N_Cases_Exome_65MB","N_Controls_Exome_65MB",
                           "Percent_Cases_Exome_65MB","Percent_Controls_Exome_65MB",
                           "N_Cases_Exome_AgilentV5","N_Controls_Exome_AgilentV5",
                           "Percent_Cases_Exome_AgilentV5","Percent_Controls_Exome_AgilentV5",
                           "N_Cases_Exome_RocheV2","N_Controls_Exome_RocheV2",
                           "Percent_Cases_Exome_RocheV2","Percent_Controls_Exome_RocheV2",
                           "N_Cases_Exome_RocheV1","N_Controls_Exome_RocheV1",
                           "Percent_Cases_Exome_RocheV1","Percent_Controls_Exome_RocheV1")
  
  write.csv(cohort_df, row.names=F, quote=F,
            file=paste0(outroot, ".info.csv"))

  # before writing broad phenotypes to file, indicate blank broad_phenotype
  # as 'undefined'
  x<-as.character(broadpheno_counts$broad_phenotype)
  broadpheno_counts$broad_phenotype <- x
  broadpheno_counts["", "broad_phenotype"] <- "undefined"

  # round broad phenotype percentages to fourth decimal place
  for (col in colnames(broadpheno_counts)) {
    if (grepl("perc_", col)) {
      broadpheno_counts[[col]] <- round(broadpheno_counts[[col]], 4)
    }
  }

  # misc adjustments to columns 
  colnames(broadpheno_counts) <- gsub("_TOTAL_", "_Meta_",
                                      colnames(broadpheno_counts)
                                     )

  # write broad phenotypes to csv as well
  write.csv(broadpheno_counts, row.names=F, quote=F,
            file=paste0(outroot, ".broad_phenotype_counts.csv")
           )

}

ReadSampleped <- function(sampleped.file) {
  ped <- read.table(sampleped.file, stringsAsFactors=F)
  colnames(ped)<- c("FID","IID","PID",
                    "MID","SEX","PHE",
                    "SEQTYPE","EXOMEPREPKIT")
  return(ped)
}

KitCounts <- function(ped.x, 
                      kits=c('65MB','AgilentV5','IDTERPv1','Roche')) {
  kit_counts <- list()
  kit_counts[["Exome"]] <- nrow(subset(ped.x,
                                       SEQTYPE!="Genome_as_fake_exome"))
  kit_counts[["Genome"]] <- nrow(subset(ped.x,
                                        SEQTYPE=="Genome_as_fake_exome"))
  
  ped.x <- subset(ped.x, SEQTYPE!="Genome_as_fake_exome")
  for (kit.i in kits) {
    ped.x.i <- subset(ped.x, EXOMEPREPKIT==kit.i)
    kit_counts[[kit.i]] <- nrow(ped.x.i)
  }
  return(kit_counts)
}

DuplicateSampleFilter <- function(ddb) {
  # given a dragendb csv, filter samples with multiple CHGVIDs by only keeping
  # the subset with the maximum prepid (since these likely represent the newest
  # sample
  chgvids <- sort(ddb$CHGVID)
  chgvids_dup <- unique(chgvids[duplicated(chgvids)])
  prepids_rm <- c()
  for (chgvid_i in chgvids_dup) {
    ddb.i <- subset(ddb, CHGVID == chgvid_i)
    prepids <- ddb.i$prepID
    ddb.i.j <- subset(ddb.i, sample_status == "In DragenDB")
    prepids <- ddb.i.j$prepID
    prepid_max <- max(prepids)
    prepids_rm <- c(prepids_rm, prepids[prepids!=prepid_max])
  }
  ddb<-subset(ddb, (prepID %in% prepids_rm)==F)
  return(ddb)
}

if (interactive() == F) {
  Main()
}
