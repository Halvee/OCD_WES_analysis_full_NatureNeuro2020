#!/bin/bash

# needed for some R scripts - access to installed packages
export R_LIBS=$PWD/lib/R/

# make sure output dir exists
mkdir -p results/clps_rate/

# run case/control rate analysis Rscript
# (full case/control cohort)
bin/Rscript 01.run_clps_rate_analysis.R \
data/clps_rate/OCD_clps_rate.allcases.sampleped \
data/clps_rate/OCD_clps_rate.allcases.evec \
data/clps_rate/OCD_clps_rate.allcases.genotypes.syn.csv \
data/clps_rate/OCD_clps_rate.allcases.genotypes.LoF.csv \
results/clps_rate/OCD_clps_rate.allcases \
> results/clps_rate/clps_rate_analysis.allcases.log \
2> results/clps_rate/clps_rate_analysis.allcases.err

exit

