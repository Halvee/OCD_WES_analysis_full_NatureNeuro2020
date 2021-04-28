#!/bin/bash
#$ -S /bin/bash
#$ -o logs/run_extTADA/run_extTADA.202006.out
#$ -e logs/run_extTADA/run_extTADA.202006.err
#$ -V
#$ -cwd
#$ -l h_rt=72:00:00
#$ -l mem_free=32G
#$ -N OCD_extTADA

# how fast is CPU clock speed?
cat /proc/cpuinfo | grep MHz

# system-specific full paths that server the following purposes :
# 1. additions to PATH and LD_LIBRARY_PATH so that extTADA runs properly
# 2. full path to R/Rscript v3.4.3, needed for extTADA to work
source cfg/exttada.cfg

# add path to R_PKGS
export R_LIBS=$PWD/lib/R/:$R_LIBS

# create Makevars file to use for building rstan
echo "CXX14 = g++" > results/extTADA/Makevars
echo "CXX14FLAGS = -std=c++14 -fPIC -O3" >> results/extTADA/Makevars

# define path to Makevars file required for using rstan
export R_MAKEVARS_USER=$PWD/results/extTADA/Makevars

# multi-caco, columns left separate
$RSCRIPT343 \
08.run_extTADA.R \
results/extTADA/OCD.771trios_476ca_1761co.final.extTADA.20000perm 72963 20000 \
> results/extTADA/OCD.771trios_476ca_1761co.final.extTADA.20000perm.log \
2> results/extTADA/OCD.771trios_476ca_1761co.final.extTADA.20000perm.err

exit
