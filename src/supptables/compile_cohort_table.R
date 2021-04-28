
## PARAM
DDB.COLS.RM <- c(
                 "individual_external_name","sample_aka",
                 "experimentID","prepID",
                 "subproject_id", "SubProject","CurrProjLeader","FundCode",
                 "Protocol","GAFbin","DetailedPhenotype","AgeAtCollection",
                 "SelfDeclGender","AvaiContUsed","ExternalData","end_point",
                 "Submission.Date","Sequencing.Complete.Date","release_status",
                 "sample_status","status_time","pseudo_prepID","seqGender",
                 "raw_bc_date","raw_bc_status","raw_bc_file",
                 "raw_bc_rginfo","raw_fc_rg_group","raw_fc_date",
                 "Relatedness_Summary","Relatedness_report","sample_type",
                 "WGS_Dragen_gVCF","MT_Mean_Cov","WGS_Mean_Cov",
                 "WGS_Median_Cov","ConcordanceWithProduction",
                 "QCMessage","pseudo_prepid","FamilyRelationProband",
                 "AlignSeqFileLoc"
                )

ARGS <- commandArgs(trailingOnly=T) 
if (length(ARGS) != 6) {
  cat("compile_cohort_table.R",
      "<combined.caco.sampleped>",
      "<rate.caco.sampleped>",
      "<triosquartets.sampleped>",
      "<fam.iids.keep.txt>",
      "<dragendb.csv>",
      "<out.csv>")

  q()
}

# get ARGS
combined.caco.sampleped <- ARGS[1]
rate.caco.sampleped <- ARGS[2]
triosquartets.sampleped <- ARGS[3]
fam.iids.keep.txt <- ARGS[4]
dragendb.csv <- ARGS[5]
out.csv <- ARGS[6]

# read samplepeds for cacos, trioquartet studies 
caco.ped <- read.table(combined.caco.sampleped, stringsAsFactors=F)
rate.ped <- read.table(rate.caco.sampleped, stringsAsFactors=F)
tq.ped <- read.table(triosquartets.sampleped, stringsAsFactors=F)
fam.iids.keep <- scan(fam.iids.keep.txt,what=character(),quiet=T)

# only trio/quartet samples in iids keep list. 
tq.ped <- subset(tq.ped, tq.ped[,2] %in% fam.iids.keep)

# assign colnames to peds
cols <-c("FamilyID","sample_internal_name","paternal_internal_name",
         "maternal_internal_name","PedGender","PedPheno","SeqType","exomeKit")
colnames(tq.ped)<-cols
colnames(caco.ped)<-cols
colnames(rate.ped)<-cols

# trim down peds : 1) take ctrls out of caco ped
caco.ped <- subset(caco.ped, PedPheno == 2)

# split tq ped into trio, quartet
fid.counts <- table(tq.ped$FamilyID)
t.fids <- names(fid.counts[fid.counts==3])
q.fids <- names(fid.counts[fid.counts>3])
t.ped <- subset(tq.ped, FamilyID %in% t.fids)
q.ped <- subset(tq.ped, FamilyID %in% q.fids)

# store which iids were in tq, caco, rate
t.iids <- t.ped$sample_internal_name
q.iids <- q.ped$sample_internal_name 
caco.iids <- caco.ped$sample_internal_name
rate.iids <- rate.ped$sample_internal_name

# read dragendb csv
ddb <- read.csv(dragendb.csv, stringsAsFactors=F)

# prunings from ddb : 1) remove RNAseq samples from table
ddb <- subset(ddb, sample_type %in% c("Exome","Genome"))

# modifications to ddb : 1) replace empty FamilyRelationProband with NA
#                        2) replace empty DNAsrc with unknown
ddb$FamilyRelationProband <- ifelse(ddb$FamilyRelationProband=="", 
                                    NA, ddb$FamilyRelationProband)
ddb$DNAsrc <- ifelse(ddb$DNAsrc == "", "unknown", ddb$DNAsrc)

# remove select columns
cols.rm <- c("FamilyID","SeqGender","SeqType","exomeKit")
for (col in cols.rm) {ddb[[col]] <- NULL}

# add list for parent internal_id to external_id 
internal.iids <- ddb$sample_internal_name
external.iids <- ddb$sample_external_name
iids.trans <- list()
for (i in 1:length(internal.iids)) {
  iids.trans[[internal.iids[i]]] <- external.iids[i]
}
tq.ped$paternal_external_name <- rep("0",nrow(tq.ped))
tq.ped$maternal_external_name <- rep("0",nrow(tq.ped))
for (i in 1:nrow(tq.ped)) {
  pid <- as.character(tq.ped[i,"paternal_internal_name"])
  if (pid != "0") {
    iid <- iids.trans[[pid]]
    tq.ped[i,"paternal_external_name"] <- iid
  }
  mid <- as.character(tq.ped[i,"maternal_internal_name"])
  if (mid != "0") {
    iid <- iids.trans[[mid]]
    tq.ped[i,"maternal_external_name"] <- iid
  }

}
caco.ped$paternal_external_name <- rep("0",nrow(caco.ped)) 
caco.ped$maternal_external_name <- rep("0",nrow(caco.ped))

# take iids out of caco ped that are already in trios quartets
caco.ped <- subset(caco.ped, 
                   (sample_internal_name %in%
                    tq.ped$sample_internal_name)==F)
full.ped <- rbind(tq.ped, caco.ped)
sample_external_names <- c()
for (iid in full.ped$sample_internal_name) {
  iid.new <- iids.trans[[iid]]
  sample_external_names <- c(sample_external_names,
                             iid.new)
}
full.ped$sample_external_name <- sample_external_names

# add columns for whether sample was included in trios, quartets,
# rate analyses
iids <- full.ped$sample_internal_name
full.ped$in_trio_dnm_analysis <- ifelse(iids %in% t.iids, "Y", "N")
full.ped$in_quartet_dnm_analysis <- ifelse(iids %in% q.iids, "Y", "N")
full.ped$in_casectrl_collapsing <- ifelse(iids %in% caco.iids, "Y", "N") 
full.ped$in_casectrl_rate <- ifelse(iids %in% rate.iids, "Y", "N")

# edit SeqType so that WGS samples listed as 'Genome' and not
# 'Genome_as_fake_exome'
full.ped$SeqType <- ifelse(full.ped$SeqType=="Genome_as_fake_exome",
                           "Genome", full.ped$SeqType)

# reorder columns
cols.ordering <- c("FamilyID",
                        "sample_internal_name", "sample_external_name",
                        "paternal_internal_name","paternal_external_name",
                        "maternal_internal_name","maternal_external_name",
                        "PedGender","PedPheno",
                        "SeqType","exomeKit","in_trio_dnm_analysis",
                        "in_quartet_dnm_analysis","in_casectrl_collapsing",
                        "in_casectrl_rate")
full.ped <- full.ped[,cols.ordering]

# add column indicating fam history
full.ped$fam_history <- ifelse(full.ped$SeqType=="Genome", "familial",
                               "unknown")
full.ped$fam_history <- ifelse((full.ped$SeqType!="Genome") &
                               ((full.ped$paternal_internal_name!="0") |
                                (full.ped$PedPheno == 1)),
                               "sporadic",full.ped$fam_history)

# prune some ddb columns before merge
for (col in DDB.COLS.RM) {ddb[[col]] <- NULL} 

# prune ddb sample_external_name
cols.rm <- c("sample_external_name")
for (col in cols.rm) {ddb[[col]] <- NULL}

# merge with dragendb table
full.db <- merge(full.ped, ddb, by='sample_internal_name')

# write to output csv
write.csv(full.db, file=out.csv,
          row.names=F, quote=T)


