#!/usr/bin/env Rscript

## LIBRARY
library(data.table)

## DATA
CLPS.CSV<-"data/prep_extTADA/OCD_collapsing_singletonquartetcases_genotypes.LoF.csv"
OUTROOT <- "results/extTADA/OCD.476ca_1761co.lof"
OUT_TSV <- "results/extTADA/OCD.771trios_476ca_1761co.extTADA.input.tsv"
MU_RATE_TSV <- "results/trio_coverage/OCDfams_2019.trios.n1.exp_dnm_rates.tsv"
TRIO_DATASET_NAMES <- c("ocdjhu2020",
                        "ocdcappi2019")
TRIO_DATASET_CDNMS <- c("results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.inkits.cdnm",
                        "results/dnm_casectrl/Cappi_etal_2019_OCD.DNMs.ppn2hdiv.inkits.cdnm")
TRIO_DATASET_NTRIOS <- c(587,
                         184)
TRIO_DATASETS_OVERLAPSCC <- c(TRUE,
                              FALSE)
CACO_DATASET_NAMES <- c("ocdjhu2020caucasian")
CACO_DATASET_NCA <- c(476)
CACO_DATASET_NCO <- c(1761)
CACO_DATASET_INROOTS <- c("OCD.476ca_1761co.")

## PARAM
CC.ANNOTS <- c("lof")

Main <- function() {

  # read csv
  df <- fread(CLPS.CSV, stringsAsFactors=F, data.table=F)

  # apply same filters used in LOF rate tests
  colnames(df) <- gsub(" ",".", colnames(df))
  df$Gene.Name <- gsub("'","",df$Gene.Name)

  # loss of function variants only
  lof.annots <-  c("frameshift_variant",
                   "splice_donor_variant",
                   "splice_acceptor_variant",
                   "stop_gained")
  df <- subset(df, Effect %in% lof.annots)
  
  # subset on het variants
  df <- subset(df, GT == "het")

  # loss of function variants only
  lof.annots <-  c("frameshift_variant",
                   "splice_donor_variant",
                   "splice_acceptor_variant",
                   "stop_gained")
  df <- subset(df, Effect %in% lof.annots)
 
  # derive subset of variants that are high quality, the same way it was
  # done in LoF rate comparisons
  for (g.filter in c("gnomAD.Exome.FILTER","gnomAD.Genome.FILTER")) {
    df[[g.filter]] <- ifelse(is.na(df[[g.filter]]), "PASS", df[[g.filter]])
    df <- df[df[[g.filter]] == "PASS", , drop=F]
  }
  df <- subset(df, FILTER %in% c("PASS"))  
  df <- subset(df, QC.Fail.Case == 0)
  df <- subset(df, QC.Fail.Ctrl == 0)

  # subset on rare variants, the same way it was done in LoF rate comparisons
  df <- subset(df, LOO.AF < 0.00001)
  for (af.col in c("gnomAD.Exome.global_AF","ExAC.global.af")) {
    df[[af.col]] <- ifelse(is.na(df[[af.col]]), 0, df[[af.col]])
    df <- df[df[[af.col]] < 0.00001, , drop=F]
  }

  # loss of function variants only
  lof.annots <-  c("frameshift_variant",
                   "splice_donor_variant",
                   "splice_acceptor_variant",
                   "stop_gained")
  df.lof <- subset(df, Effect %in% lof.annots)
  
  # make qualvar df, write to file
  qv <- df.lof[,c("Gene.Name", "Variant.ID", "Sample.Name")]
  qv$NALLELE <- rep(1, nrow(qv))
  write.table(qv, row.names=F, col.names=F, sep="\t", quote=F,
              file=paste0(OUTROOT, ".qualvar"))

  # read mutation rate tsv
  mu <- read.table(MU_RATE_TSV, header=T, stringsAsFactors=F)
  genes <- mu$Gene

  # collect counts of synon, damaging missense and LoF DNMs for each
  # gene symbol, for each independent DNM cohort
  dnm.counts <- list()
  dnm.varids <- c()
  for (i in 1:length(TRIO_DATASET_NAMES)) {
    name.i <- TRIO_DATASET_NAMES[i]
    cdnm.file.i <- TRIO_DATASET_CDNMS[i]
    cdnm.i <- read.table(cdnm.file.i, stringsAsFactors=F)
    cdnm.i[,8] <- gsub("splice|frameshift|non","lof",cdnm.i[,8])
    if (TRIO_DATASETS_OVERLAPSCC[i] == TRUE) {
      dnm.varids <- c(dnm.varids, cdnm.i[,4])
    }
    lof.counts <- table(cdnm.i[cdnm.i[,8] == "lof", 7])
    misnd.counts <- table(cdnm.i[cdnm.i[,8] %in%
                                 c("mis","misB","misP","misU"), 7])
    misd.counts <- table(cdnm.i[cdnm.i[,8] == "misD", 7])
    syn.counts <- table(cdnm.i[cdnm.i[,8] == "syn", 7])
    dnm.counts[[name.i]] <- list("syn"=syn.counts,
                                 "misND"=misnd.counts,
                                 "misD"=misd.counts,
                                 "lof"=lof.counts)
  }
  dnm.varids <- unique(dnm.varids)

  # for each row in caco dataset, read corresponding qualvar table,
  # filter out variant IDs overlapping any DNM calls,
  # get case/control proportions for each gene symbol
  caco.dataset.names <- c()
  ca.counts <- list()
  co.counts <- list()
  for (i in 1:length(CACO_DATASET_NAMES)) {
    name.i <- CACO_DATASET_NAMES[i]
    nca.i <- CACO_DATASET_NCA[i]
    nco.i <- CACO_DATASET_NCO[i]
    qv.j <- qv

    for (ann.j in CC.ANNOTS) {

      # split into case/control sets of calls
      qv.ca.j <- qv.j[grepl("ocd",qv.j[,3])==T, ]
      qv.co.j <- qv.j[grepl("ocd",qv.j[,3])==F, ]

      # strip away known dnm calls from cases
      qv.ca.j <-  qv.ca.j[(qv.ca.j[,2] %in% dnm.varids) == F, ] 

      # get counts of cases/controls with >= 1 var per gene
      qv.ca.j <- unique(qv.ca.j[,c(1,3)])
      qv.co.j <- unique(qv.co.j[,c(1,3)])
      ca.counts.j <- table(qv.ca.j[,1])
      co.counts.j <- table(qv.co.j[,1])
      if ((name.i %in% names(ca.counts))==F) {
        ca.counts[[name.i]] <- list()
        co.counts[[name.i]] <- list()
      }
      ca.counts[[name.i]][[ann.j]]<- ca.counts.j
      co.counts[[name.i]][[ann.j]]<- co.counts.j
  
    }

  }
  caco.dataset.names <- unique(CACO_DATASET_NAMES)

  # init output dataframe for reading into extTADA, store non-zero values
  # to dataframe from DNM and ca/co counts
  out.df <- data.frame(gene=mu$Gene, 
                       mut_syn=mu$syn, mut_misND=mu$misND, mut_misD=mu$misD, mut_lof=mu$lof)
  rownames(out.df) <- out.df$gene
  for (name.i in TRIO_DATASET_NAMES) {
    for (annot.j in c("syn","misND","misD","lof")) {
      col <- paste("dn",name.i, annot.j, sep="_")
      out.df[[col]] <- rep(0,nrow(out.df))
      gene.dnm.counts <- dnm.counts[[name.i]][[annot.j]]
      gene.dnm.1 <- names(gene.dnm.counts)
      gene.dnm.1 <- intersect(gene.dnm.1, out.df$gene)
      for (gene in gene.dnm.1) {
        out.df[gene, col] <- gene.dnm.counts[[gene]]
      }
    }
  }
  for (name.i in caco.dataset.names) {
    
    # annotations to be added to table
    annots <- CC.ANNOTS
  
    # add case columns
    for (annot.j in annots) {
      
      # add case column
      ca.col <- paste("cc",name.i, annot.j,"ca", sep="_")
      out.df[[ca.col]] <- rep(0,nrow(out.df))
      gene.ca.counts <- ca.counts[[name.i]][[annot.j]] 
      gene.ca.1 <- names(gene.ca.counts) 
      for (gene in gene.ca.1) {
        out.df[gene, ca.col] <- as.numeric(gene.ca.counts[[gene]])
      }

      # add control column
      co.col <- paste("cc",name.i, annot.j,"co", sep="_")
      out.df[[co.col]] <- rep(0, nrow(out.df))
      gene.co.counts <- co.counts[[name.i]][[annot.j]]
      gene.co.1 <- names(gene.co.counts)
      for (gene in gene.co.1) {
        out.df[gene, co.col] <- as.numeric(gene.co.counts[[gene]])
      }

    }
  }

  # make sure any 'NA' rows are pruned out
  out.df <- subset(out.df, is.na(gene)==F)

  # write complete data frame to tsv
  write.table(out.df, file=OUT_TSV,
              row.names=F, col.names=T,
              quote=F, sep="\t")

}

if (interactive() == F) {
  Main()
}
