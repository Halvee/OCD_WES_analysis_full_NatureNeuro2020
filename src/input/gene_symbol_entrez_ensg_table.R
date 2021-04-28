#!/usr/bin/env Rscript

## LIBRARY
library(org.Hs.eg.db)

## PARAM
ENSG_PRUNE_MANUAL <- c(
                       "ENSG00000218497", # ZNF84
                       "ENSG00000232748", #ZNF668
                       "ENSG00000236575", #ZNF26
                       "ENSG00000270386", #UGT2A1
                       "ENSG00000184040", #TMEM236
                       "ENSG00000215699", #SRSF10
                       "ENSG00000270466", #SEPT1
                       "ENSG00000269846", #RBL1
                       "ENSG00000215700", #PNRC2
                       "ENSG00000233024", #NPIPA7
                       "ENSG00000183889", #NPIPA7
                       "ENSG00000262621", #NAA60
                       "ENSG00000272781", #MDGA2
                       "ENSG00000272658", #LTB4R2
                       "ENSG00000269099", #LSP1
                       "ENSG00000257207", #LIMS3
                       "ENSG00000217792", #KIR3DL2
                       "ENSG00000231880", #KBTBD4
                       "ENSG00000269783", #GOLGA7B
                       "ENSG00000220903", #FRG2C
                       "ENSG00000233050", #DEFB130
                       "ENSG00000268942", #CKS1B
                       "ENSG00000181464", #CDRT1
                       "ENSG00000255994" #CCDC177
                      )
ENTREZID_PRUNE_MANUAL <- c("654780", "105372382","100124696",
                           "100505381","51072","105376839",
                           "100187828")

Main <- function(){
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 3) {
    cat("ensg_symbol_entrez_table.R <ensg_symbol.tsv> <CCDSr20.geneset> <out.tsv>\n")
    q()
  }
  ensg_symbol.tsv<-ARGS[1]
  ccdsr20.geneset <- ARGS[2]
  out.tsv <- ARGS[3]

  # read input ccds file
  ccds <- scan(ccdsr20.geneset, what=character())

  # read input annotation table
  ensg_symbol <- read.table(ensg_symbol.tsv, header=F, stringsAsFactors=F)
  colnames(ensg_symbol) <- c("ENSEMBL","SYMBOL")
  ensg_symbol <- unique(ensg_symbol)
  print(dim(ensg_symbol))

  # only keep rows where the symbol is in ccds geneset
  ensg_symbol <- subset(ensg_symbol, SYMBOL %in% ccds)

  # use ensg to extract corresponding entrez IDs
  x<-select(org.Hs.eg.db, keys=ensg_symbol$SYMBOL,
            columns=c("ENTREZID"), keytype="SYMBOL")
  print(dim(x))
  ensg_symbol_entrez <- merge(ensg_symbol, x, by="SYMBOL")

  # only retain rows where ensg was in original table
  ensg_symbol_entrez <- subset(ensg_symbol_entrez,
                               ENSEMBL %in% ensg_symbol$ENSEMBL)
  
  # only keep entries with an entrez ID
  ensg_symbol_entrez <- subset(ensg_symbol_entrez,
                               is.na(ENTREZID)==F)
  print(dim(ensg_symbol_entrez))
  print(head(ensg_symbol_entrez))

  # remove any entrez IDs that are in manual prune list
  ensg_symbol_entrez <- subset(ensg_symbol_entrez,
                               (ENTREZID %in% ENTREZID_PRUNE_MANUAL)==F)
  
  # remove any ensg IDs that are in manual prune list
  ensg_symbol_entrez <- subset(ensg_symbol_entrez,
                               (ENSEMBL %in% ENSG_PRUNE_MANUAL)==F)

  # make sure that only unique entry rows are in df
  ensg_symbol_entrez <- unique(ensg_symbol_entrez)
  colnames(ensg_symbol_entrez) <- c("GeneSymbol","ENSG","EntrezID")

  # write to tsv
  write.table(ensg_symbol_entrez, row.names=F, col.names=T, sep="\t", quote=F,
              file=out.tsv)
}

if (interactive() == F) {
  Main()
}
