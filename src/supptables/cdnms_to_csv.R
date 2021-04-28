

## LIBRARY
library(org.Hs.eg.db)


Main <- function() {

  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 6) {
    cat("cdnms_to_csv.R",
        "<trios.cdnm>>",
        "<quartets.cdnm>",
        "<controls.cdnm>",
        "<cappi_trios.cdnm>",
        "<symbol_ensg_entrezid.tsv>",
        "<out.csv>\n")
    q()
  }

  # get ARGS
  trios.cdnm.file <- ARGS[1]
  quartets.cdnm.file <- ARGS[2]
  controls.cdnm.file <- ARGS[3]
  cappi_trios.cdnm.file <- ARGS[4]
  symbol_ensg_entrezid.tsv <- ARGS[5]
  out.csv <- ARGS[6]

  # read files
  trios.cdnm <-read.table(trios.cdnm.file, stringsAsFactors=F)
  quartets.cdnm <- read.table(quartets.cdnm.file, stringsAsFactors=F)
  controls.cdnm <- read.table(controls.cdnm.file, stringsAsFactors=F)
  cappi_trios.cdnm <- read.table(cappi_trios.cdnm.file, stringsAsFactors=F)
  symbol_ensg_entrezid <- read.table(symbol_ensg_entrezid.tsv,
                                     header=T, sep="\t", stringsAsFactors=F)

  # combine cdnm dfs, add column in front that labels source of dnm call
  quartets.cdnm<-unique(quartets.cdnm)
  cdnm <- rbind(trios.cdnm, quartets.cdnm, controls.cdnm, cappi_trios.cdnm)
  colnames(cdnm)<-c("Sample_Name","Chrom","Pos","Variant_ID",
                    "Ref_Allele","Alt_Allele","Gene","EffectShort",
                    "In_Jointly_Covered_Loci")
  cdnm.x <- data.frame(Cohort=c(rep("OCD_JHU_trios", nrow(trios.cdnm)),
                                rep("OCD_JHU_quartets",nrow(quartets.cdnm)),
                                rep("control_Iossifov2014_trios",
                                    nrow(controls.cdnm)),
                                rep("OCD_Cappi2019_trios",
                                    nrow(cappi_trios.cdnm))))
  cdnm <- cbind(cdnm.x, cdnm)

  # add Entrez_ID and Ensembl_Gene_ID columns to cdnm df
  colnames(symbol_ensg_entrezid)[1] <- c("Gene")
  cdnm <- merge(cdnm, symbol_ensg_entrezid, by = "Gene")
  cdnm <- cdnm[,c("Cohort","Sample_Name","Chrom","Pos","Variant_ID",
                  "Ref_Allele","Alt_Allele",
                  "EffectShort","Gene","ENSG","EntrezID",
                  "In_Jointly_Covered_Loci")]

  # write to csvs
  write.csv(cdnm, file=out.csv,
            row.names=F, quote=T)

}

if (interactive() == F) {
  Main()
}
