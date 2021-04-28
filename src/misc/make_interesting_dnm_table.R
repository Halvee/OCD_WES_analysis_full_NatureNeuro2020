

Main <- function(){
  ARGS <-commandArgs(trailingOnly=T)
  if (length(ARGS)!=8) {
    cat("make_interesting_dnm_table.R <loeuf_bins.tsv> <mpc_gt_2.varlist>",
        "<NDD.geneset> TS.geneset> <trios.cdnm> <quartets.cdnm>",
        "<cappi_trios.cdnm> <out.csv>\n")
    q()
  }

  # get ARGS
  loeuf_bins.tsv <- ARGS[1]
  mpc_gt_2.varlist <- ARGS[2]
  ndd.geneset.file <- ARGS[3]
  ts.geneset.file <- ARGS[4]
  trios.cdnm.file <- ARGS[5]
  quartets.cdnm.file <- ARGS[6]
  cappi_trios.cdnm.file <- ARGS[7]
  out.csv <- ARGS[8]

  # read input files
  loeuf_bins <- read.table(loeuf_bins.tsv,header=T,stringsAsFactors=F)
  mpc_gt_2 <- scan(mpc_gt_2.varlist,what=character())
  ndd.geneset <- scan(ndd.geneset.file,what=character())
  ts.geneset <- scan(ts.geneset.file,what=character())
  trios.cdnm <- read.table(trios.cdnm.file, stringsAsFactors=F)
  quartets.cdnm <- read.table(quartets.cdnm.file, stringsAsFactors=F)
  cappi_trios.cdnm <- read.table(cappi_trios.cdnm.file, stringsAsFactors=F)

  # assign cdnm colnames
  cols <- c("IID","CHR","POS","VARID","REF","ALT","GENE","EFF","INKITS")
  colnames(trios.cdnm) <- cols
  colnames(quartets.cdnm) <- cols
  colnames(cappi_trios.cdnm) <- cols

  # unique rows only
  trios.cdnm <- unique(trios.cdnm)
  print(nrow(quartets.cdnm))
  quartets.cdnm <- unique(quartets.cdnm)

  # get subset of cdnms where variant is lof in a gene with loeuf<10%
  lof.effs <- c("non","splice","frameshift")
  loeuf1.genes <- subset(loeuf_bins, oe_lof_upper_bin==0)$gene
  t.cdnm.lof_loeuf1 <- subset(trios.cdnm,
                              (EFF %in% lof.effs) & (GENE %in% loeuf1.genes))
  q.cdnm.lof_loeuf1 <- subset(quartets.cdnm,
                              (EFF %in% lof.effs) & (GENE %in% loeuf1.genes))
  print(dim(t.cdnm.lof_loeuf1))
  print(dim(q.cdnm.lof_loeuf1)) 

  # get subset of cdnms with a misD variant with MPC > 2
  t.cdnm.misd_mpcgt2 <- subset(trios.cdnm,
                               (EFF=="misD") & (VARID %in% mpc_gt_2))
  q.cdnm.misd_mpcgt2 <- subset(quartets.cdnm,
                               (EFF=="misD") & (VARID %in% mpc_gt_2))

  # get subset of cdnms with a misD or LOF variant in NDD genes (n=187 genes)
  misdlof.effs <- c("misD","non","splice","frameshift")
  t.cdnm.misdlof_ndd <- subset(trios.cdnm,
                               (EFF %in% misdlof.effs) & 
                               (GENE %in% ndd.geneset)
                              )
  q.cdnm.misdlof_ndd <- subset(quartets.cdnm,
                               (EFF %in% misdlof.effs) & 
                               (GENE %in% ndd.geneset)
                              )
  
  # get subset of cdnms with a misD or LOF variant in TS genes (n=6 genes)
  misdlof.effs <- c("misD","non","splice","frameshift")
  t.cdnm.misdlof_ts <- subset(trios.cdnm,
                              (EFF %in% misdlof.effs) & 
                              (GENE %in% ts.geneset)
                             )
  q.cdnm.misdlof_ts <- subset(quartets.cdnm,
                              (EFF %in% misdlof.effs) & 
                              (GENE %in% ts.geneset)
                             )
 
  # get subset of cdnms where gene impacted is hit >=2x with a misD/LOF DNM
  misdlof.effs <- c("misD","non","splice","frameshift")
  genes.exclude <- c("TTN")
  trios_all.cdnm <- rbind(trios.cdnm, cappi_trios.cdnm)
  misdlof.genes <- subset(trios_all.cdnm, (EFF %in% misdlof.effs))$GENE
  misdlof.genes <- unique(misdlof.genes[duplicated(misdlof.genes)==T])
  misdlof.genes <- misdlof.genes[(misdlof.genes %in% genes.exclude)==F]
  t.cdnm.misdlof_recurr <- subset(trios.cdnm,
                                  (EFF %in% misdlof.effs) &
                                  (GENE %in% misdlof.genes)
                                 )
  q.cdnm.misdlof_recurr <- subset(quartets.cdnm,
                                  (EFF %in% misdlof.effs) &
                                  (GENE %in% misdlof.genes)
                                 )
  print(misdlof.genes)
  print(t.cdnm.misdlof_recurr)
  print(q.cdnm.misdlof_recurr)

  # get variant ids for all of the above
  qualvars <- c(t.cdnm.lof_loeuf1$VARID, q.cdnm.lof_loeuf1$VARID,
                t.cdnm.misd_mpcgt2$VARID, q.cdnm.misd_mpcgt2$VARID,
                t.cdnm.misdlof_ndd$VARID, q.cdnm.misdlof_ndd$VARID,
                t.cdnm.misdlof_ts$VARID, q.cdnm.misdlof_ts$VARID,
                t.cdnm.misdlof_recurr$VARID, q.cdnm.misdlof_recurr$VARID)
  qualvars <- unique(sort(qualvars))
  t.qv <- subset(trios.cdnm, VARID %in% qualvars)
  q.qv <- subset(quartets.cdnm, VARID %in% qualvars)

  # subset columns, add extra that tells which cohort the var is from
  cols.keep <- c("IID","COHORT","CHR","POS","VARID","REF","ALT","GENE","EFF")
  t.qv$COHORT <- rep("Trios", nrow(t.qv))
  q.qv$COHORT <- rep("Quartets",nrow(q.qv))
  t.qv <- t.qv[,cols.keep]
  q.qv <- q.qv[,cols.keep]
  print(dim(q.qv)) 
  print(dim(t.qv)) 
  qv <- rbind(t.qv, q.qv)

  # add columns for variant membership in groups of varids
  qv[["LOF/LOEUF<10%"]] <- ifelse(qv$VARID %in% 
                                  c(t.cdnm.lof_loeuf1$VARID,
                                    q.cdnm.lof_loeuf1$VARID),
                                  "YES",
                                  "NO")
  qv[["misD/MPC>2"]] <- ifelse(qv$VARID %in% 
                                c(t.cdnm.misd_mpcgt2$VARID,
                                  q.cdnm.misd_mpcgt2$VARID),
                                "YES",
                                "NO")
  qv[["misDLOF/NDD genes"]] <- ifelse(qv$VARID %in% 
                                      c(t.cdnm.misdlof_ndd$VARID,
                                        q.cdnm.misdlof_ndd$VARID),
                                      "YES",
                                      "NO")
  qv[["misDLOF/TS genes"]] <- ifelse(qv$VARID %in% 
                                      c(t.cdnm.misdlof_ts$VARID,
                                        q.cdnm.misdlof_ts$VARID),
                                      "YES",
                                      "NO")
  qv[["misDLOF/recurrently hit genes"]]<-ifelse(qv$VARID %in%
                                                c(t.cdnm.misdlof_recurr$VARID,
                                                  q.cdnm.misdlof_recurr$VARID),
                                                "YES",
                                                "NO")
  qv[["AffectedSiblingShared"]] <- ifelse(qv$VARID == "2-240036868-G-A", "YES", "NO")

  # write to csv
  write.csv(qv, row.names=F, quote=F,
            file=out.csv)




  q()

}

if (interactive() == F) {
 Main()
}
