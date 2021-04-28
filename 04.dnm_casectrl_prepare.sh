#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o logs/dnm_casectrl_prepare/dnm_casectrl_prepare.out
#$ -e logs/dnm_casectrl_prepare/dnm_casectrl_prepare.err
#$ -V
#$ -N OCD.dnm_casectrl_prepare

## BIN

## DATA
CASE_DNMS=(
           "results/input/OCDfams_2019.cohort.triosquartets.10x.dnm"
           "results/input/Cappi_etal_2019_OCD.DNMs.tsv"
          )
CTRL_DNMS=(
           "results/input/Iossifov_etal_2014_healthysib_DNMs.tsv"
           "results/input/Iossifov_etal_2014_sibling_sharedDNMs.tsv"
          )
KIT_INTERSECT_BED=../RocheV2_coverage/OCD_RocheV2_addjCCDSr20_highcov_intersect.0.9.bed
SCORE_LIS=results/mut_rate_tbls/hs37d5.CCDSr20.dbNSFP2.9.scores_ppn2hdiv.lis
RGN_FILE=data/input/addjusted.CCDS.genes.index.r20.hg19.txt

# needed for some python scripts - access to installed packages
export PYTHONPATH="$PWD/lib/python/:$PYTHONPATH"

# needed for some R scripts - access to installed packages
export R_LIBS=$PWD/lib/R/

# store clineff license path to variable
source data/clineff_snpsift/clineff_license_path

## make sure output dir exists
mkdir -p results/dnm_casectrl/

DNMS="${CASE_DNMS[@]} ${CTRL_DNMS[@]}"

for DNM in ${DNMS[@]}
do
  
  # init outroot for all files
  OUTROOT=`basename $DNM | sed 's/\.dnm$//g; s/\.tsv$//g'`

  # make VCF file from DNM file
  awk '{
    OFS="\t"; 
    if ((length($5)==1) || (length($5)!=length($6))) {
      print $2,$3,$4,$5,$6,".",".","."
    }
  }' $DNM |
    sort -k1,1 -k2,2n |
    bin/python src/dnm_casectrl_prepare/vcf_rgn_intersect.py \
    stdin \
    $RGN_FILE \
    > results/dnm_casectrl/$OUTROOT.ccds.vcf

  # annotate VCF using clineff
  ANNOT_VCF=results/dnm_casectrl/$OUTROOT.ccds.annotated.vcf
  ANNOT_LOG=results/dnm_casectrl/$OUTROOT.ccds.AnnotatedVCF.log
  echo "annotating $OUTROOT.ccds.vcf."
  echo "output VCF : $OUTROOT.ccds.annotated.vcf"
  echo "log file : $OUTROOT.ccds.AnnotatedVCF.log"

  # go to dir with clineff and workflow to execute command
  cd data/clineff_snpsift/

  # execute clineff command
  cmd="../../bin/java \
       -Xmx6g \
       -jar ../../bin/ClinEff.jar \
       -v \
       -license $CLINEFF_LICENSE \
       GRCh37.87 \
       ../../results/dnm_casectrl/$OUTROOT.ccds.vcf \
       > ../../results/dnm_casectrl/$OUTROOT.ccds.annotated.vcf \
       2> ../../results/dnm_casectrl/$OUTROOT.ccds.AnnotatedVCF.log"
  echo $cmd
  eval $cmd

  # go back to parent dir
  cd ../../

  # convert to cdnm
  cmd="bin/python src/dnm_casectrl_prepare/dnm_vcfann_to_cdnm.py \
       results/mut_rate_tbls/GRCh37.87.ccds.enst.txt \
       $DNM \
       results/dnm_casectrl/$OUTROOT.ccds.annotated.vcf \
       results/dnm_casectrl/$OUTROOT.cdnm"
  echo $cmd
  eval $cmd

done

## annotate each cdnm file using precomputed matrix of scores (Polyphen2 HDIV)
## alongside removing variants that have a POPMAX MAF >= 0.0005 in 
## gnomAD exome dataset
cmd="bin/Rscript src/dnm_casectrl_prepare/get_scores_cdnm_mafprune.R \
     $SCORE_LIS \
     data/input/gnomAD_exome_popmaxMAFge0.0005.varlist \
     results/dnm_casectrl/OCDfams_2019.cohort.triosquartets.10x.cdnm results/dnm_casectrl/OCDfams_2019.cohort.triosquartets.10x.ppn2hdiv.cdnm \
     results/dnm_casectrl/Cappi_etal_2019_OCD.DNMs.cdnm results/dnm_casectrl/Cappi_etal_2019_OCD.DNMs.ppn2hdiv.cdnm \
     results/dnm_casectrl/Iossifov_etal_2014_healthysib_DNMs.cdnm results/dnm_casectrl/Iossifov_etal_2014_healthysib_DNMs.ppn2hdiv.cdnm \
     results/dnm_casectrl/Iossifov_etal_2014_sibling_sharedDNMs.cdnm results/dnm_casectrl/Iossifov_etal_2014_sibling_sharedDNMs.ppn2hdiv.cdnm"
echo $cmd
eval $cmd

# add column indicating whether variants falling SeqCap v2 / SeqCap v3 /
# IDTERPv1 intersect regions
cmd="bin/python src/dnm_casectrl_prepare/bed_cdnm_intersect.py \
     $KIT_INTERSECT_BED \
     results/dnm_casectrl/OCDfams_2019.cohort.triosquartets.10x.ppn2hdiv.cdnm results/dnm_casectrl/OCDfams_2019.cohort.triosquartets.10x.ppn2hdiv.inkits.cdnm \
     results/dnm_casectrl/Cappi_etal_2019_OCD.DNMs.ppn2hdiv.cdnm results/dnm_casectrl/Cappi_etal_2019_OCD.DNMs.ppn2hdiv.inkits.cdnm \
     results/dnm_casectrl/Iossifov_etal_2014_healthysib_DNMs.ppn2hdiv.cdnm results/dnm_casectrl/Iossifov_etal_2014_healthysib_DNMs.ppn2hdiv.inkits.cdnm \
     results/dnm_casectrl/Iossifov_etal_2014_sibling_sharedDNMs.ppn2hdiv.cdnm results/dnm_casectrl/Iossifov_etal_2014_sibling_sharedDNMs.ppn2hdiv.inkits.cdnm"
echo $cmd
eval $cmd

## split OCD cdnm into trios, quartets
cmd="bin/Rscript src/dnm_casectrl_prepare/cdnm_trios_quartets_split.R \
     results/dnm_casectrl/OCDfams_2019.cohort.triosquartets.10x.ppn2hdiv.cdnm \
     results/input/OCDfams_2019.cohort.triosquartets.sampleped \
     results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.cdnm \
     results/dnm_casectrl/OCDfams_2019.cohort.quartets.10x.ppn2hdiv.cdnm"
echo $cmd
eval $cmd
cmd="bin/Rscript src/dnm_casectrl_prepare/cdnm_trios_quartets_split.R \
     results/dnm_casectrl/OCDfams_2019.cohort.triosquartets.10x.ppn2hdiv.inkits.cdnm \
     results/input/OCDfams_2019.cohort.triosquartets.sampleped \
     results/dnm_casectrl/OCDfams_2019.cohort.trios.10x.ppn2hdiv.inkits.cdnm \
     results/dnm_casectrl/OCDfams_2019.cohort.quartets.10x.ppn2hdiv.inkits.cdnm"
echo $cmd
eval $cmd

exit
