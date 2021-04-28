

## PARAM
DOT.RM <- T

# read args
ARGS <- commandArgs(trailingOnly=T)
if (length(ARGS) != 2) {
  cat("tsv_to_csv.R <in.tsv> <out.csv>\n")
  q()
}

# read tsv
df <-read.table(ARGS[1], header=T, 
                check.names=F,
                sep="\t", stringsAsFactors=F)

# replace dots with underscores?
if (DOT.RM == T) {
  colnames(df) <- gsub("\\.","_",colnames(df))
}

# write csv
write.csv(df, file=ARGS[2], row.names=F, quote=F)
