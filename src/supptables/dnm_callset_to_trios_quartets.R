

## LIBRARY
library(org.Hs.eg.db)

# excess DNM count (> 5)
FIDS.PRUNE <- c("ocdl")
# treated like a trio but actually are quartets. were excluded from analysis.
FIDS.PRUNE <- c(FIDS.PRUNE, c("ocda","ocdg","ocdk","ocdv"))

ARGS <- commandArgs(trailingOnly=T)
if (length(ARGS) != 9) {
  cat("callset_to_trios_quartets.R",
      "<trios.sampleped>",
      "<quartets.sampleped>",
      "<callset.csv>",
      "<sanger.csv>",
      "<trios.cdnm>>",
      "<quartets.cdnm>",
      "<symbol_ensg_entrezid.tsv>",
      "<callset.trios.csv>",
      "<callset.quartets.csv>\n")
  q()
}

# get ARGS
trios.sampleped <- ARGS[1]
quartets.sampleped <- ARGS[2]
callset.csv <- ARGS[3]
sanger.csv <- ARGS[4]
trios.cdnm.file <- ARGS[5]
quartets.cdnm.file <- ARGS[6]
symbol_ensg_entrezid.tsv <- ARGS[7]
callset.out.trios.csv <- ARGS[8]
callset.out.quartets.csv <- ARGS[9]

# read files
t.ped <- read.table(trios.sampleped, stringsAsFactors=F)
q.ped <- read.table(quartets.sampleped, stringsAsFactors=F)
t.cdnm <-read.table(trios.cdnm.file, stringsAsFactors=F)
q.cdnm <- read.table(quartets.cdnm.file, stringsAsFactors=F)
callset <- read.csv(callset.csv, stringsAsFactors=F)
sanger <- read.csv(sanger.csv, stringsAsFactors=F)
symbol_ensg_entrezid <- read.table(symbol_ensg_entrezid.tsv,
                                   header=T, sep="\t", stringsAsFactors=F)

# replace dots in callset columns with underscores
colnames(callset) <- gsub("\\.","_",colnames(callset))

# exclude families from callset that were excluded from trio/quartet analyses
callset <- subset(callset, (Family_ID %in% FIDS.PRUNE)==F)

# remove extraneous columns
cols.rm <-  c(
              "HGMDm2site","HGMDm1site","HGMD_site","HGMD_Disease","HGMD_PMID",
              "HGMD_Class","HGMDp1site","HGMDp2site","HGMD_indel_9bpflanks",
              "ClinVar","ClinVar_RS_Number","ClinVar_Disease",
              "ClinVar_Clinical_Significance","ClinVar_PMID","ClinVar_pathogenic_indels","ClinVar_all_indels","ClinVar_Pathogenic_Indel_Count","Clinvar_Pathogenic_CNV_Count","ClinVar_Pathogenic_SNV_Splice_Count","ClinVar_Pathogenic_SNV_Nonsense_Count","ClinVar_Pathogenic_SNV_Missense_Count","ClinVar_Pathogenic_Last_Pathogenic_Location","ClinGen","ClinGen_HaploinsufficiencyDesc","ClinGen_TriplosensitivityDesc","OMIM_Disease","RecessiveCarrier","ACMG","dbDSM_Disease","dbDSM_Classification","dbDSM_PubmedID","DenovoDB_Phenotype","DenovoDB_PubmedID",
              "Evs_All_Maf","Evs_All_Genotype_Count","Evs_Filter_Status",
              "All_Effect_Gene_Transcript_HGVS_p_Polyphen_Humdiv_Polyphen_Humvar",
              "X","Chrom","Pos","Ref","Alt","Ref_Allele","Alt_Allele",
              "Polyphen_Humdiv_Score", "Polyphen_Humdiv_Prediction",
              "Polyphen_Humdiv_Prediction__CCDS_","Polyphen_Humvar_Score",
              "Polyphen_Humvar_Prediction","Polyphen_Humvar_Score__CCDS_",
              "Polyphen_Humvar_Prediction__CCDS_")
for (col in cols.rm) {
  callset[[col]] <- NULL
}

# add column names to cdnm files
q.cdnm<-unique(q.cdnm)
cdnm <- rbind(t.cdnm, q.cdnm)
colnames(cdnm)<-c("Sample_Name","Chrom","Pos","Variant_ID",
                  "Ref_Allele","Alt_Allele","Gene","EffectShort",
                  "In_Jointly_Covered_Loci")

# add Entrez_ID and Ensembl_Gene_ID columns to cdnm df
colnames(symbol_ensg_entrezid)[1] <- c("Gene")
cdnm <- merge(cdnm, symbol_ensg_entrezid, by = "Gene")
cdnm <- cdnm[,c("Sample_Name","Chrom","Pos","Variant_ID",
                "Ref_Allele","Alt_Allele",
                "EffectShort","Gene","ENSG","EntrezID",
                "In_Jointly_Covered_Loci")]

# form samplevar columns
callset$SampleName_VariantID <- paste(callset$Sample_Name,
                                      callset$Variant_ID,
                                      sep=":")
cdnm$SampleName_VariantID <- paste(cdnm$Sample_Name, 
                                   cdnm$Variant_ID,
                                   sep=":")
sanger$VariantID <- paste(sanger$Chr, sanger$Pos, 
                          sanger$Ref, sanger$Alt,
                          sep="-")
sanger$SampleName_VariantID <- paste0(sanger$Sample.Name,
                                      ":",
                                      sanger$VariantID)

# add a column to cdnm df for status of sanger validation.
# status will be one of the following :
# not_attempted : no sanger validation attempted for the DNM.
# unassessable : sanger data unable to be completely assessed
#                due to one of or some combination of the following:
#                1. absence of data in proband or parents
#                2. sanger data that is too messy to be assessed
#                   due to a failed run, sample mixing, etc
# assessed_pass : de novo mutation assessed as present in at least
#                 one of two strands in a proband, and absent
#                 in both parents in at least one strand.

# derive add sanger status column to cdnm
cdnm$Sanger_Validation <- rep("not_attempted", nrow(cdnm))
sanger.pass <- subset(sanger, 
                      (assessable_snv_pass==1) |
                      (assessable_indel_pass==1)
                     )
cdnm$Sanger_Validation <- ifelse(cdnm$SampleName_VariantID %in% 
                                 sanger.pass$SampleName_VariantID,
                                 "assessable_pass",
                                 cdnm$Sanger_Validation)
sanger.assess <- subset(sanger,
                        (total_snv==1 & assessable_snv==0) |
                        (total_indel==1 & assessable_indel==0)
                       )
cdnm$Sanger_Validation <- ifelse(cdnm$SampleName_VariantID %in% 
                                 sanger.assess$SampleName_VariantID,
                                 "unassessable",
                                 cdnm$Sanger_Validation)

# merge cdnm df with callset df
callset$SampleName_VariantID <- paste(callset$Sample_Name,
                                      callset$Variant_ID,
                                      sep=":")
cdnm$SampleName_VariantID <- paste(cdnm$Sample_Name, 
                                   cdnm$Variant_ID,
                                   sep=":")
callset$Sample_Name <- NULL
callset$Variant_ID <- NULL
callset <- merge(cdnm, callset, by="SampleName_VariantID")
callset$SampleName_VariantID <- NULL

# get trio and quartet fids
t.fids <- unique(t.ped[,1])
q.fids <- unique(q.ped[,1])

# split callset into trios, quartets
t.callset <- subset(callset, Family_ID %in% t.fids)
q.callset <- subset(callset, Family_ID %in% q.fids)

# write to csvs
write.csv(t.callset, file=callset.out.trios.csv,
          row.names=F, quote=T)
write.csv(q.callset, file=callset.out.quartets.csv,
          row.names=F, quote=T)

