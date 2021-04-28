#!/bin/bash

# get to the extra R packages you'll need
export R_LIBS=$PWD/lib/R/

# run dnm rate caco analysis Rscript
bin/Rscript 05.dnm_rate_caco.R \
> results/dnm_rate_caco/dnm_rate_caco.out \
2> results/dnm_rate_caco/dnm_rate_caco.err

exit
