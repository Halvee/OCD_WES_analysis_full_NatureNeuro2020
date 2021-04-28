
Main <- function() {
  ARGS <- commandArgs(trailingOnly = T)
  if (length(ARGS) != 4) {
    cat("cdnm_trios_quartets_split.R <in.cdnm>",
        "<in.sampleped> <out.trios.cdnm> <out.quartets.cdnm>\n")
    q()
  }

  in.cdnm <- ARGS[1]
  in.sampleped <- ARGS[2]
  out.trios.cdnm <- ARGS[3]
  out.quartets.cdnm <-ARGS[4]

  cdnm <- read.table(in.cdnm, stringsAsFactors=F)
  ped <- read.table(in.sampleped, stringsAsFactors=F)

  # get trio, quartet iids
  fids <- ped[,1]
  fid.counts <- table(fids)
  fids.t <- names(fid.counts[fid.counts==3])
  fids.q <- names(fid.counts[fid.counts>3])
  iids.t <- ped[ped[,1] %in% fids.t, 2]
  iids.q <- ped[ped[,1] %in% fids.q, 2]

  # split cdnm into trio-only, quartet-only
  cdnm.t <- cdnm[cdnm[,1] %in% iids.t, , drop=F]
  cdnm.q <- cdnm[cdnm[,1] %in% iids.q, , drop=F] 

  # write to output file
  write.table(cdnm.t, file=out.trios.cdnm,
              row.names=F, col.names=F, sep="\t", quote=F)
  write.table(cdnm.q, file=out.quartets.cdnm,
              row.names=F, col.names=F, sep="\t", quote=F)

}

if (interactive() == F) {
  Main()
}
