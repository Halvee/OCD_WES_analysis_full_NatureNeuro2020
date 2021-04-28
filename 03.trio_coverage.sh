#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o logs/trio_coverage/trio_coverage.out
#$ -e logs/trio_coverage/trio_coverage.err
#$ -V
#$ -l mem_free=32G
#$ -l h_rt=24:00:00
#$ -N OCD.trio_coverage

## FILES
TRIOS_SAMPLEPED=results/input/OCDfams_2019.cohort.trios.sampleped
PARENT_GENESET=results/input/CCDSr20_autosomal.geneset
MU_LIS=results/mut_rate_tbls/hs37d5.CCDSr20.mu.lis
SCORE_LIS=results/mut_rate_tbls/hs37d5.CCDSr20.dbNSFP2.9.scores_ppn2hdiv.gnomADexomeMAFlt0.0005.lis
COVDIR=../OCDfams_2018_trio_cov_10x20x/results/jointcov_rgn/

## PARAM
COV="10x"

# needed for some python scripts - access to installed packages
export PYTHONPATH="$PWD/lib/python/:$PYTHONPATH"

# needed for some R scripts - access to installed packages
export R_LIBS=$PWD/lib/R/

# create subsets of sampleped that are roche kit only, IDTERPv1 only
grep "Roche" results/input/OCDfams_2019.cohort.trios.sampleped \
> results/trio_coverage/OCDfams_2019.trios.RocheOnly.sampleped
grep "IDTERPv1" results/input/OCDfams_2019.cohort.trios.sampleped \
> results/trio_coverage/OCDfams_2019.trios.IDTERPv1Only.sampleped

# assemble text files where each line is a full path to a .rgn where intervals
# are joint coverage >= 10 or 20x in the referenced trio
bin/python src/trio_coverage/assemble_covpaths.py \
$COVDIR \
$TRIOS_SAMPLEPED \
$COV \
results/trio_coverage/OCDfams_2019.trios.jointcov_covpaths.$COV.txt

# assemble covpaths specific to Roche, IDTERPv1 kits
bin/python src/trio_coverage/assemble_covpaths.py \
$COVDIR \
results/trio_coverage/OCDfams_2019.trios.RocheOnly.sampleped \
$COV \
results/trio_coverage/OCDfams_2019.trios.jointcov_covpaths.$COV.RocheOnly.txt
bin/python src/trio_coverage/assemble_covpaths.py \
$COVDIR \
results/trio_coverage/OCDfams_2019.trios.IDTERPv1Only.sampleped \
$COV \
results/trio_coverage/OCDfams_2019.trios.jointcov_covpaths.$COV.IDTERPv1Only.txt

# compute joint coverage metrics (full cohort)
bin/python src/trio_coverage/collect_joint_cov_metrics.py \
results/trio_coverage/OCDfams_2019.trios.jointcov_covpaths.$COV.txt \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV

# compute joint coverage metrics (Roche, IDTERPv1 only)
bin/python src/trio_coverage/collect_joint_cov_metrics.py \
results/trio_coverage/OCDfams_2019.trios.jointcov_covpaths.$COV.RocheOnly.txt \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.RocheOnly
bin/python src/trio_coverage/collect_joint_cov_metrics.py \
results/trio_coverage/OCDfams_2019.trios.jointcov_covpaths.$COV.IDTERPv1Only.txt \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.IDTERPv1Only

# collect summary statistics for coverage (full cohort)
bin/Rscript src/trio_coverage/coverage_summary_stats.R \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.sampleProfileSummary.tsv \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.genicProfileSummary.tsv \
> results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.coverage_summary_stats.log

# collect summary statistics for coverage (Roche, IDTERPv1 only)
bin/Rscript src/trio_coverage/coverage_summary_stats.R \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.RocheOnly.sampleProfileSummary.tsv \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.RocheOnly.genicProfileSummary.tsv \
> results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.RocheOnly.coverage_summary_stats.log
bin/Rscript src/trio_coverage/coverage_summary_stats.R \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.IDTERPv1Only.sampleProfileSummary.tsv \
results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.IDTERPv1Only.genicProfileSummary.tsv \
> results/trio_coverage/OCDfams_2019.trios.cvg_stats.$COV.IDTERPv1Only.coverage_summary_stats.log

# for each CCDS base, compute the number of input trios that meet joint cov specifications
bin/python src/trio_coverage/cov_samples_per_pos.py \
$MU_LIS \
results/trio_coverage/OCDfams_2019.trios.jointcov_covpaths.$COV.txt \
$TRIOS_SAMPLEPED \
results/trio_coverage/OCDfams_2019.trios.$COV.ncov

# from mutation rate table and joint cov metrics compute mutation rate for each gene,
# taking into account the total num trios covered >= 10x at each site
bin/python src/mut_rate_tbls/compute_exp_dnm_rates.py \
results/input/CCDSr20_autosomal.geneset \
$MU_LIS \
$SCORE_LIS \
results/trio_coverage/OCDfams_2019.trios.$COV.ncov \
> results/trio_coverage/OCDfams_2019.trios.$COV.exp_dnm_rates.tsv

# from mutation rate table and joint cov metrics compute mutation rate for each gene,
# taking into account the total num trios covered >= 10x at each site
bin/python src/mut_rate_tbls/compute_exp_dnm_rates.py \
results/input/CCDSr20.geneset \
$MU_LIS \
$SCORE_LIS \
1 \
> results/trio_coverage/OCDfams_2019.trios.n1.exp_dnm_rates.tsv

exit
