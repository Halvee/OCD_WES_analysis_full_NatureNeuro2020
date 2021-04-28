

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 1 & length(ARGS) != 2) {
    cat("cohort_general_stats.R <in.ddb.csv> [in.sampleped]\n")
    q()
  }

  # read dragendb csv
  in.ddb.csv <- ARGS[1]
  ddb <- read.csv(in.ddb.csv, stringsAsFactors=F)

  # if defined, read sampleped, subset ddb on it
  if (length(ARGS) == 2) {
    in.sampleped <- ARGS[2]
    ped <- read.table(in.sampleped, stringsAsFactors=F)
    ddb <- subset(ddb, CHGVID %in% ped[,2])
  }

  # subset on samples in dragendb
  ddb <- subset(ddb, sample_status=="In DragenDB")

  # make sure samples are OCD, not nyspi
  ddb <- subset(ddb, BroadPhenotype=="obsessive compulsive disorder")
  ddb <- subset(ddb, SeqType %in% c("Exome","Genome"))
  ddb <- subset(ddb, grepl("nyspiocd", ddb$CHGVID) == F)

  # print stats (total number of cases)
  cat("Total number of cases :",nrow(ddb), "\n")
  
}

if (interactive() == F) {
  Main()
}
