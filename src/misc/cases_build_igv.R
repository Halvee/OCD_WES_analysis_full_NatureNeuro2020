
# GLOBAL VARS
FLANK <- 60
SAMPLEID.PREFIX <- "ocd"
ORIG.ID.COL <- "OrigID"
IGM.ID.COL <- "CHGVID"
PREP.ID.COL <- "pseudo_prepid"
ALNSEQLOC.COL <- "AlignSeqFileLoc"
SEQDB.COLS.KEEP <- c(ORIG.ID.COL, IGM.ID.COL)

main <- function() {
  ARGS <- commandArgs(trailingOnly = T)
  if (length(ARGS) != 6) {
    cat("DNM_setup_output.R <proj_name> <cohort.sampleped> <cohort.seqdb.csv>",
        "<samplevar.tsv> <samtools.path>",
        "<outdir>\n")
    q()
  }

  proj.name <- ARGS[1]
  ped.file <- ARGS[2]
  seqdb.csv <- ARGS[3]
  samplevar.tsv <- ARGS[4]
  samtools.path <- ARGS[5]
  outdir <- ARGS[6]

  # read input (ped, seqdb and samplevar files)
  ped <- read.table(ped.file, stringsAsFactors=F)
  seqdb <- read.csv(seqdb.csv, stringsAsFactors=F) 
  samplevar <- read.table(samplevar.tsv, stringsAsFactors=F)
  colnames(samplevar)<- c("Sample.Name","Variant.ID")

  # only keep samplevar, seqdb entries where proband is found in ped file
  samplevar<-subset(samplevar, Sample.Name %in% ped[,2])
  seqdb <- subset(seqdb, CHGVID %in% ped[,2])

  # read seqdb csv file and get bam locations
  seqdb <- read.csv(seqdb.csv, stringsAsFactors=F)
  seqdb <- subset(seqdb, grepl(SAMPLEID.PREFIX, seqdb[[IGM.ID.COL]]) == T)
  seqdb$bampath <- paste0(seqdb[[ALNSEQLOC.COL]], "/", 
                          seqdb[[IGM.ID.COL]], ".", seqdb[[PREP.ID.COL]], "/",  
                          seqdb[[IGM.ID.COL]], ".", seqdb[[PREP.ID.COL]], ".realn.recal.bam")
  seqdb <- seqdb[, c(SEQDB.COLS.KEEP, "bampath")]

  # from seqdb df, store bampaths
  bampaths <- list()
  for (i in 1:nrow(seqdb)) {
    bampaths[[ seqdb[i, IGM.ID.COL] ]] <- seqdb[i, "bampath"]
  }

  # form bam outroot as outdir + proj.name
  bam.outroot <- paste0(outdir, "/", proj.name)
  
  # write bash script for getting bam slices as input to IGV
  bam.paths <- c()
  bam.paths.child <- c()
  interval.strs <- c()
  fileConn<-file(paste0(bam.outroot, ".samtools.sh"), "w")
  writeLines("#!/bin/bash", fileConn)
  writeLines("#$ -cwd", fileConn) 
  writeLines(paste0("#$ -o logs/",proj.name,".bamslices.samtools.out"), fileConn)
  writeLines(paste0("#$ -e logs/",proj.name,".bamslices.samtools.err"), fileConn)
  writeLines("#$ -V", fileConn)
  writeLines("#$ -pe threaded 1", fileConn)
  writeLines(paste0("#$ -N ",proj.name,"_bamslices_samtools"), fileConn) 
  for (i in 1:nrow(samplevar)) {
    sample.name <- samplevar[i, "Sample.Name"]
    var <- as.character(samplevar[i, "Variant.ID"])
    var.info <- strsplit(var, "-")[[1]]
    chrom <- var.info[1]
    pos <- as.integer(var.info[2])
    pos.left <- pos - FLANK
    pos.right <- pos + FLANK
    interval.str <- paste0(chrom, ":", pos.left, "-", pos.right)
    bam.out <- paste0(bam.outroot, ".", sample.name, "_", var, ".aln.bam")
    writeLines(paste(samtools.path, "view -b", bampaths[[sample.name]], "-o", bam.out, interval.str), fileConn)
    writeLines(paste(samtools.path, "index", bam.out), fileConn)

    # store info for creation of igv batch file
    bam.out.str <- paste0("BAMDIRECTORY/",basename(bam.out))
    bam.paths <- c(bam.paths, bam.out.str)
    interval.strs <- c(interval.strs, interval.str)
  }
  close(fileConn)

  # write IGV script for visualizing single vars
  fileConn<-file(paste0(bam.outroot, ".igv.batch"), "w")
  writeLines(c("genome hg19", "snapshotDirectory SNAPSHOTDIRECTORY"), fileConn)
  for ( i in 1:length(bam.paths)) {
    bam.path <- bam.paths[i]
    interval.str <- interval.strs[i]
    aln.png <- gsub(".bam", ".aln.png", bam.path) 
    aln.png <- gsub("BAMDIRECTORY/","", aln.png)
    writeLines("new", fileConn)
    writeLines(paste0("load ", bam.path), fileConn)
    writeLines(paste("goto", interval.str), fileConn)
    writeLines(c("sort position", "collapse"), fileConn)
    writeLines(paste("snapshot", aln.png), fileConn)
  }
  writeLines("exit", fileConn)
  close(fileConn)

}

if (interactive() == F) {
  main()
}
