
## PARAM
EXTTADA.CACO.COLS <- c("cc_ocdjhu2020caucasian_lof_ca",
                       "cc_ocdjhu2020caucasian_lof_co")

Main <- function(){
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 6) {
    cat("merge_dnm_test_results.R <ensg_symbol_entrezid.tsv>",
        "<dnm.misd.res.tsv>",
        "<dnm.misdlof.res.tsv>< <extTADA.in.tsv> <extTADA.out.tsv> <out.csv>\n")
    q()
  }
  ensg_symbol_entrezid.tsv<-ARGS[1]
  dnm.lof.res.tsv<-ARGS[2]
  dnm.misdlof.res.tsv<-ARGS[3]
  exttada.in.tsv<-ARGS[4]
  exttada.res.tsv <- ARGS[5]
  out.csv<-ARGS[6]
 
  # read input tables
  exttada.in <-read.table(exttada.in.tsv, header=T,sep="\t",
                          stringsAsFactors=F)
  ensg_symbol_entrezid <- read.table(ensg_symbol_entrezid.tsv,
                                     header=T,sep="\t",stringsAsFactors=F)
  colnames(ensg_symbol_entrezid) <- c("gene","ensg","entrez_id")
  dnm.lof.res <-read.table(dnm.lof.res.tsv,header=T,stringsAsFactors=F)
  dnm.lof.res<-dnm.lof.res[,c("gene","p.value")]
  colnames(dnm.lof.res)[2] <- "dn_lof_poisson_p"
  dnm.misdlof.res <-read.table(dnm.misdlof.res.tsv,header=T,stringsAsFactors=F)
  dnm.misdlof.res<-dnm.misdlof.res[,c("gene","p.value")]
  colnames(dnm.misdlof.res)[2] <- "dn_lofmisD_poisson_p"
  exttada.res <- read.table(exttada.res.tsv, header=T, 
                            stringsAsFactors=F)

  # keep DNM counts per annot, merge with ensg/entrezid
  exttada.in <- exttada.in[,c("gene","mut_misD","mut_lof",
                              "dn_ocdjhu2020_misD","dn_ocdjhu2020_lof",
                              "dn_ocdcappi2019_misD","dn_ocdcappi2019_lof")]
  df <- merge(ensg_symbol_entrezid, exttada.in, by='gene')
 
  # add total de novo counts
  df$dn_misD <- df$dn_ocdjhu2020_misD + df$dn_ocdcappi2019_misD
  df$dn_lof <- df$dn_ocdjhu2020_lof + df$dn_ocdcappi2019_lof

  # merge DNM-only rate tests into df
  df <- merge(df, dnm.lof.res, by="gene")
  df <- merge(df, dnm.misdlof.res, by="gene")
  print(head(df))
  print(head(exttada.res))
  
  # retrieve extTADA results and merge into dataframe
  exttada.res <- exttada.res[,c("Gene","BF","PP","qvalue",
                                EXTTADA.CACO.COLS)]
  colnames(exttada.res) <- c("gene","extTADA_BF","extTADA_PP","extTADA_qvalue",
                             EXTTADA.CACO.COLS)
  df <- merge(df, exttada.res, by='gene')
  df <- df[order(df$extTADA_qvalue), ]
  print(head(df))

  #gene  mut_syn mut_misND mut_misD  mut_lof dn_ocdjhu2020_syn
  #dn_ocdjhu2020_misND dn_ocdjhu2020_misD  dn_ocdjhu2020_lof dn_ocdcappi2019_syn
  #dn_ocdcappi2019_misND dn_ocdcappi2019_misD  dn_ocdcappi2019_lof
  
  # write to csv
  write.csv(df, row.names=F, quote=F,
            file=out.csv) 

}

if (interactive() == F) {
  Main()
}
