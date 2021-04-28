#!/bin/bash

# make sure output dir exists
mkdir -p results/clinical/

# run clinical analysis script
bin/Rscript 09.clinical.R \
> results/clinical/clinical.analysis.log
2> results/clinical/clinical.analysis.err

# supplemental analysis :
# LoF rate in male versus female OCD cases
# (non-trio samples, representative of a comparison
# independent of the DNM burden tests above)
bin/Rscript 01.run_clps_rate_analysis.R \
data/clinical/ocd.nontriocases.malevsfemale.sampleped \
data/clinical/ocd.nontriocases.malevsfemale.evec \
data/clinical/ocd.nontriocases.malevsfemale.syn.csv \
data/clinical/ocd.nontriocases.malevsfemale.lof.csv \
results/clinical/ocd.nontriocases.malevsfemale.lofrate \
> results/clinical/ocd.nontriocases.malevsfemale.lofrate.log \
2>results/clinical/ocd.nontriocases.malevsfemale.lofrate.err

