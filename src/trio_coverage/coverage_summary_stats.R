
Main <- function() {
  
  # get ARGS
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 2) {
    cat("coverage_summary_stats.R <in.samplecov.tsv> <in.genecov.tsv>\n")
    q()
  }
  samplecov.tsv <- ARGS[1]
  genecov.tsv <- ARGS[2]

  # coverage data
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

}

if (interactive() == F) {
  Main()
}
