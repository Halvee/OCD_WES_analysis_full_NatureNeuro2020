

Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 4) {
    cat("gnomad_maf_prune.R <gnomad_mafs.tsv> <maf.max> <in.cdnm> <out.cdnm\n")
    q()
  }

  gnomad_mafs.tsv <- ARGS[1]
  maf.max <- as.numeric(ARGS[2])
  in.cdnm <- ARGS[3]
  out.cdnm <- ARGS[4]

  # read files
  gnomad <- read.table(gnomad_mafs.tsv, header=T, stringsAsFactors=F)
  cdnm <- read.table(in.cdnm, stringsAsFactors=F, sep="\t")

  # get all gnomad vars with popmax MAF > MAF_MAX
  vars.rm <- c()
  for (col in c("gnomAD_AF",
                "gnomAD_AFR_AF",
                "gnomAD_AMR_AF",
                "gnomAD_ASJ_AF",
                "gnomAD_EAS_AF",
                "gnomAD_FIN_AF",
                "gnomAD_NFE_AF",
                "gnomAD_OTH_AF",
                "gnomAD_SAS_AF")) {
    gnomad[[col]] <- gsub("^-$", "0", gnomad[[col]])  
    gnomad[[col]] <- as.numeric(gnomad[[col]])
    gnomad.i <- gnomad[gnomad[[col]] >= maf.max, , drop=F]
    vars.rm <- c(vars.rm, gnomad.i$Uploaded_variation)
  }
  vars.all <- unique(sort(gnomad$Uploaded_variation))
  vars.rm <- unique(sort(vars.rm))
  cdnm <- cdnm[(cdnm[,4] %in% vars.rm) == F, , drop=F]

  # write pruned cdnm to file
  write.table(cdnm,
              file = out.cdnm,
              row.names=F, col.names=F,
              sep="\t", quote=F)
}

if (interactive() == F) {
  Main()
}
