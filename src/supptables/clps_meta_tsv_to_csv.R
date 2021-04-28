

## PARAM
DOT.RM <- T

# read args
ARGS <- commandArgs(trailingOnly=T)
if (length(ARGS) != 3) {
  cat("clps_meta_tsv_to_csv.R <in.tsv> <symbol_ensg_entrez.tsv> <out.csv>\n")
  q()
}
in.tsv <- ARGS[1]
symbol_ensg_entrez.tsv <- ARGS[2]
out.csv <- ARGS[3]

# read tsvs
df <-read.table(in.tsv, header=T, 
                check.names=F,
                sep="\t", stringsAsFactors=F)
symbol_ensg_entrez <- read.table(symbol_ensg_entrez.tsv, header=T, sep="\t",
                                 stringsAsFactors=F)

# replace dots with underscores?
if (DOT.RM == T) {
  colnames(df) <- gsub("\\.","_",colnames(df))
}

# merge results with ensg / enstrez id table, resort by cmh p
colnames(symbol_ensg_entrez) <- c("gene","ensg","entrez_id")
df <- merge(symbol_ensg_entrez, df, by="gene")
df <- df[order(df$cmh_p), ]

# write csv
write.csv(df, file=out.csv, row.names=F, quote=F)
