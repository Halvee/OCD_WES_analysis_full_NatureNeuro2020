#!/bin/bash

# make sure output dir exists
mkdir -p results/extTADA/

# get to the extra R packages you'll need
export R_LIBS=$PWD/lib/R/

# run gene-based dnm rate test Rscript
bin/Rscript 07.prep_extTADA.R \
> results/extTADA/prep_extTADA.out \
2> results/extTADA/prep_extTADA.err

exit
