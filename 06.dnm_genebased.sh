#!/bin/bash

# get to the extra R packages you'll need
export R_LIBS=$PWD/lib/R/

# run gene-based dnm rate test Rscript
bin/Rscript 06.dnm_genebased.R \
> results/dnm_genebased/dnm_genebased.out \
2> results/dnm_genebased/dnm_genebased.err

exit
