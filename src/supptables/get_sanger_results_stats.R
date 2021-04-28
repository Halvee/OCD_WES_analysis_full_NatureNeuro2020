
Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 3) {
    cat("sanger_results_stats.R <sanger_results.csv>",
        "<sanger_results_pruned.csv> <stats.csv\n")
    q()
  }

  # read files
  sanger_res <- read.csv(ARGS[1], stringsAsFactors=F)
  out_sanger_csv <- ARGS[2]
  out_sanger_stats_csv <- ARGS[3]

  # remove suboptimal DNMs from sanger_res
  sanger_res <- subset(sanger_res, notSuboptimal == 1)
  sanger_res$notSuboptimal <- NULL

  # place quotes back into columns where they're needed
  sanger_res$RESULTS_PROBAND <- paste0("\"", sanger_res$RESULTS_PROBAND, "\"")
  sanger_res$RESULTS_PARENTS <- paste0("\"", sanger_res$RESULTS_PARENTS, "\"")

  # write pruned sanger results to csv
  write.csv(sanger_res, row.names=F, quote=F,
            file=out_sanger_csv)

  # form output df
  out_df <- data.frame(SUMS=c("N_TOTAL","N_AVAILABLE","N_ASSESSABLE",
                              "N_ASSESSABLE_PASS",
                              "PERC_AVAILABLE_PASS","PERC_ASSESSABLE_PASS"),
                       SNV=rep(NA,6), INDEL=rep(NA,6), SNVINDEL=rep(NA,6))
  rownames(out_df) <- out_df$SUMS

  # get SNV counts
  res_snv <- subset(sanger_res, isSNV==1)
  out_df["N_TOTAL","SNV"] <- nrow(res_snv)
  res_snv <- subset(res_snv, proband_data_present==1 & n_parent_data_present==2)
  out_df["N_AVAILABLE","SNV"] <- nrow(res_snv)
  res_snv <- subset(res_snv, 
                    proband_data_conclusive==1 & 
                    n_parent_data_conclusive==2)
  out_df["N_ASSESSABLE","SNV"] <- nrow(res_snv)
  res_snv <- subset(res_snv,
                    proband_confirmed_carrier==1 &
                    n_parent_confirmed_noncarrier==2)
  out_df["N_ASSESSABLE_PASS","SNV"] <- nrow(res_snv)

  # get indel counts
  res_indel <- subset(sanger_res, isIndel==1)
  out_df["N_TOTAL","INDEL"] <- nrow(res_indel)
  res_indel <- subset(res_indel, proband_data_present==1 & n_parent_data_present==2)
  out_df["N_AVAILABLE","INDEL"] <- nrow(res_indel)
  res_indel <- subset(res_indel, 
                    proband_data_conclusive==1 & 
                    n_parent_data_conclusive==2)
  out_df["N_ASSESSABLE","INDEL"] <- nrow(res_indel)
  res_indel <- subset(res_indel,
                    proband_confirmed_carrier==1 &
                    n_parent_confirmed_noncarrier==2)
  out_df["N_ASSESSABLE_PASS","INDEL"] <- nrow(res_indel)

  # get total snv/indel counts
  for (row.i in c("N_TOTAL","N_AVAILABLE","N_ASSESSABLE","N_ASSESSABLE_PASS")){
    out_df[row.i, "SNVINDEL"] <- out_df[row.i,"SNV"] + out_df[row.i,"INDEL"]
  }

  # get percentage of available mutations that are passing
  cols <- c("SNV","INDEL","SNVINDEL")
  out_df["PERC_AVAILABLE_PASS",cols] <- out_df["N_ASSESSABLE_PASS",cols] /
                                        out_df["N_AVAILABLE",cols]

  # get percentage of assessable mutations that are passing
  cols <- c("SNV","INDEL","SNVINDEL")
  out_df["PERC_ASSESSABLE_PASS",cols] <- out_df["N_ASSESSABLE_PASS",cols] /
                                     out_df["N_ASSESSABLE",cols]

  # round percentages
  cols <- c("SNV","INDEL","SNVINDEL")
  for (row.i in c("PERC_AVAILABLE_PASS","PERC_ASSESSABLE_PASS")) {
    out_df[row.i,cols] <- round(out_df[row.i,cols],3)
  }

  # write to output csv
  write.csv(out_df, file=out_sanger_stats_csv, row.names=F, quote=F)
}

if (interactive() == F) {
  Main()
}
