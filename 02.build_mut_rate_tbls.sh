#!/bin/bash
#$ -S /bin/bash -cwd
#$ -o logs/build_mut_rate_tbls/build_mut_rate_tbls.out
#$ -e logs/build_mut_rate_tbls/build_mut_rate_tbls.err
#$ -l mem_free=32G
#$ -l h_rt=24:00:00
#$ -V
#$ -N OCD.build_Hg19_CCDSr20_mut_rates

## PARAM
N_NT_FLANK=1

# make sure output dir exists
mkdir -p results/mut_rate_tbls/

# user needs to make sure pythonpath contains dirs with required module 'twobitreader'
export PYTHONPATH="$PWD/lib/python/:$PYTHONPATH"

# get ENST IDs that are a part of CCDS
zcat data/input/Homo_sapiens.GRCh37.87.gtf.gz |
  awk -F"\t" '{if ($3=="transcript") {print $9}}' |
  grep 'tag "CCDS"' |
  sed 's/\; /;/g' |
  tr ";" "\n" |
  grep transcript_id |
  cut -d" " -f2 |
  tr -d '"' \
  > results/mut_rate_tbls/GRCh37.87.ccds.enst.txt

# make .2bit file for hs37d5 genome reference
bin/faToTwoBit \
data/input/hs37d5.fa \
results/mut_rate_tbls/hs37d5.2bit

# make mutation rate table across CCDSr20
bin/python \
src/mut_rate_tbls/build_mu_tbl.py \
data/input/mutation_rate_by_trinucleotide_matrix.txt \
$N_NT_FLANK \
results/mut_rate_tbls/hs37d5.2bit \
data/input/addjusted.CCDS.genes.index.r20.hg19.txt \
> results/mut_rate_tbls/hs37d5.CCDSr20.mu.lis

# form single unsorted vcf from mu tbl
tail -n+2 results/mut_rate_tbls/hs37d5.CCDSr20.mu.lis |
  awk '{OFS="\t";
        print $2,$3,$4,"A";
        print $2,$3,$4,"C";
        print $2,$3,$4,"G";
        print $2,$3,$4,"T"}' |
  awk '{OFS="\t"; if ($3!=$4) {
    print $1,$2,$1"-"$2"-"$3"-"$4,$3,$4,".",".","."
  }}' > results/mut_rate_tbls/hs37d5.CCDSr20.unsorted.vcf

# init sorted vcf, sort each chrom and append to vcf
> results/mut_rate_tbls/hs37d5.CCDSr20.vcf
for chrom in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
  echo $chrom
  awk -v CHR=$chrom '{if ($1==CHR) {print $0}}' results/mut_rate_tbls/hs37d5.CCDSr20.unsorted.vcf |
    sort -k1,1 -k2,2n |
    uniq \
    >> results/mut_rate_tbls/hs37d5.CCDSr20.vcf

done

# remove intermediate files
rm results/mut_rate_tbls/hs37d5.CCDSr20.unsorted.vcf

# run clineff on all possible SNVs in CCDSr20 loci (Ensembl 87 CCDS transcripts)
cd data/input/

../../bin/java \
  -Xm20g \
  -jar ../../bin/ClinEff.jar \
  -license clinEff.license \
  GRCh37.87 \
  ../../results/mut_rate_tbls/hs37d5.CCDSr20.vcf \
  2> ../../results/mut_rate_tbls/hs37d5.CCDSr20.AnnotateVCF.log \
| ../../bin/bgzip -c \
> ../../results/mut_rate_tbls/hs37d5.CCDSr20.annotated.vcf.gz

cd ../../

# annotate with dbnsfp
bin/java \
-Xmx20g \
-jar bin/SnpSift.jar dbnsfp \
-db data/input/dbNSFP2.9.txt.gz \
-f Polyphen2_HDIV_score,Polyphen2_HDIV_pred \
-v results/mut_rate_tbls/hs37d5.CCDSr20.annotated.vcf.gz \
2> results/mut_rate_tbls/hs37d5.CCDSr20.annotated.dbNSFP2.9.log \
| bin/bgzip -c \
> results/mut_rate_tbls/hs37d5.CCDSr20.annotated.dbNSFP2.9.vcf.gz

# for each input SNV in VCF, derive maximum annotation :
# 3 == synon, 2 == LoF (stop-gain OR splice donor/acceptor),
# 0 otherwise
# write as output table in same format as mu table
zcat results/mut_rate_tbls/hs37d5.CCDSr20.annotated.vcf.gz \
| bin/python src/mut_rate_tbls/init_score_table.py \
  results/input/CCDSr20.geneset \
  results/mut_rate_tbls/GRCh37.87.ccds.enst.txt \
  stdin \
  results/mut_rate_tbls/hs37d5.CCDSr20.mu.lis \
> results/mut_rate_tbls/hs37d5.CCDSr20.annot_synlof.lis

echo "syn / mis / lof annotation done."

# recode missense scores in table based on PPN2 HDIV scores from dbnsfp
zcat results/mut_rate_tbls/hs37d5.CCDSr20.annotated.dbNSFP2.9.vcf.gz \
| bin/python src/mut_rate_tbls/score_tbl_missense_dbnsfp.py \
 stdin \
 results/mut_rate_tbls/hs37d5.CCDSr20.annot_synlof.lis \
> results/mut_rate_tbls/hs37d5.CCDSr20.dbNSFP2.9.scores_ppn2hdiv.lis

echo "missense annotation done."

# recode sites to '4' if found at MAF >= 0.0005 in at least one gnomAD 
# subpopulation
bin/python src/mut_rate_tbls/score_tbl_exclude.py \
results/mut_rate_tbls/hs37d5.CCDSr20.dbNSFP2.9.scores_ppn2hdiv.lis \
data/input/gnomAD_exome_popmaxMAFge0.0005.varlist \
> results/mut_rate_tbls/hs37d5.CCDSr20.dbNSFP2.9.scores_ppn2hdiv.gnomADexomeMAFlt0.0005.lis

# produce table of gene-based de novo mutation rates per annotation
bin/python src/mut_rate_tbls/compute_exp_dnm_rates.py \
results/input/CCDSr20.geneset \
results/mut_rate_tbls/hs37d5.CCDSr20.mu.lis \
results/mut_rate_tbls/hs37d5.CCDSr20.dbNSFP2.9.scores_ppn2hdiv.gnomADexomeMAFlt0.0005.lis \
1 \
> results/mut_rate_tbls/hs37d5.CCDSr20.dbNSFP2.9.n1.exp_dnm_rates.tsv

exit
