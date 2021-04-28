
Main <- function() {
  ARGS <- commandArgs(trailingOnly=T)
  if (length(ARGS) != 3) {
    cat("make_extTADA_csv.R <extTADA.input.tsv> <extTADA.output.tsv>",
        "<extTADA.output_formatted.csv>\n")
    q()

    # read extTADA input, output tsvs
    
    # rearrange so that sumstats are in the first few columns,
    # de novo counts are listed per cohort (JHU, Cappi),
    # lof counts are listed per each of the 6 case/ctrl cohorts

    # note : need to include ENSG or entrez id alongside gene symbols

    # write to csv
  }
}

if (interactive() == F) {
  Main()
}
