#!/usr/bin/env Rscript

## PARAM
GENE_COUNTS_TSV="results/extTADA/OCD.771trios_476ca_1761co.extTADA.input.tsv"
NTRIOS <- 587+184
NCAS <- c(476
         )
NCOS <- c(1761
         )
NCA<-sum(NCAS)
NCO<-sum(NCOS)
TRIO_DATASETS <- c("ocdjhu2020","ocdcappi2019")
CACO_DATASET <- "ocdjhu2020caucasian"
CACO_DATASETS <- c("ocdjhu2020caucasian")
EXTTADA_SCRIPT_DIR="../src/extTADA/script/"
OUTROOT <- "results/extTADA/OCD.771trios_476ca_1761co.final.extTADA"
NO_FRAMESHIFT_DNMS <- T
SWITCH_LOF_RATES_DENOVOLYZER <- F

# incorporate multi-core processing
options(mc.cores = parallel::detectCores())

ARGS <- commandArgs(trailingOnly=T)
if (length(ARGS) >= 1) {
  OUTROOT <- ARGS[1]
}
if (length(ARGS) >= 2) {
  set.seed(ARGS[2])
} 
if (length(ARGS) >= 3) {
  NITER <- as.numeric(ARGS[3])
}

## switchout input gene counts tsv if provided as arg
fileR <- dir(EXTTADA_SCRIPT_DIR, ".R$")
for (ii in fileR) {
  source(paste0(EXTTADA_SCRIPT_DIR, ii))
}

## make rstan run faster?
set_cppo(mode="fast")

N = list(dn = NTRIOS, ca=NCA, cn=NCO)
inputData <- read.table(GENE_COUNTS_TSV, 
                        header = TRUE, as.is = TRUE)
inputData <- inputData[is.na(inputData[,1]) == F, ]

## if desired by user switch out lof mutation rate with denovolyzeR vals
if (SWITCH_LOF_RATES_DENOVOLYZER == T) {
  library(denovolyzeR)
  mu.tbl <- viewProbabilityTable()
  mu.tbl <- mu.tbl[, c("hgncSymbol","lof")]
  colnames(mu.tbl) <- c("gene", "mut_lof")
  inputData$mut_lof <- NULL
  inputData <- merge(inputData, mu.tbl)

}

## combine de novo mutation calls across JHU and Cappi datasets
inputData$dn_misD <- rowSums(inputData[,paste0("dn_",TRIO_DATASETS,"_misD")])
inputData$dn_lof <- rowSums(inputData[,paste0("dn_",TRIO_DATASETS,"_lof")])

## get rid of rows where for whatever reason an expected mutation rate is zero
inputData <- subset(inputData, mut_misD > 0)
inputData <- subset(inputData, mut_lof > 0)

inputData$mut_misD <- inputData$mut_misD 
inputData$mut_lof <- inputData$mut_lof 

allDNData <- inputData[, paste0("dn_", c("misD", "lof"))]
allMutData <- inputData[,paste0("mut_", c("misD", "lof"))]

head(data.frame(allMutData, allDNData))
print(table(inputData[["dn_misD"]]))
print(table(inputData[["dn_lof"]]))

## make sure input cols are numeric
for (i in colnames(allDNData)) {allDNData[[i]] <- as.numeric(allDNData[[i]])}
for (i in colnames(allMutData)) {
  allMutData[[i]] <- as.numeric(allMutData[[i]])
  print(min(allMutData[[i]]))
  print(max(allMutData[[i]]))
}

## form merged cc counts
inputData$cc_OCDjhu_lof_ca <- rowSums(inputData[,paste0("cc_",CACO_DATASETS,"_lof_", "ca"),drop=F])
inputData$cc_OCDjhu_lof_co <- rowSums(inputData[,paste0("cc_",CACO_DATASETS,"_lof_", "co"),drop=F])

## collect ca/co data
dataCCcase <- inputData[,paste0("cc_",CACO_DATASETS,"_lof_ca"),drop=F]
dataCCcontrol <- inputData[,paste0("cc_",CACO_DATASETS,"_lof_co"),drop=F]

mcmcDD <- extTADAmcmc(modelName = DNandCCextTADA, #extTADA for only de novo data
                      dataDN = allDNData, mutRate = allMutData,
                      Ndn = rep(N$dn, 2), ## There are two de novo categories
            		      dataCCcase = dataCCcase, dataCCcontrol = dataCCcontrol,
		                  Ncase = NCAS, Ncontrol = NCOS,
                      nCore = 3, nChain = 3,
                      nIteration = NITER)

options(repr.plot.width=5, repr.plot.height=5)
stan_trace(mcmcDD)

mcmcDD 

gamma.cc<-c('hyperGammaMeanCC[1]')
beta.cc<-c('hyperBetaCC[1]')
pars0 <- estimatePars(pars = c('pi0',
                               'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]',
                               gamma.cc, beta.cc),
                      mcmcResult = mcmcDD)

pars0 

options(repr.plot.width=4, repr.plot.height=3)
par(mfrow = c(1, 2))
plotParHeatmap(pars = c("pi0", "hyperGammaMeanDN[1]"), mcmcResult = mcmcDD)
plotParHeatmap(pars = c("pi0", "hyperGammaMeanDN[2]"), mcmcResult = mcmcDD)

##Get gene list
geneName <- inputData[, 1]
##Set parameters: use pars0 above
parsFDR <- list(gammaMeanDN = pars0[, 1][c('hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]')],
	        betaDN = pars0[, 1][c('hyperBetaDN[1]', 'hyperBetaDN[2]')],
		gammaMeanCC = pars0[,1][gamma.cc],
		betaCC = pars0[, 1][beta.cc],
		pi0 = pars0[, 1][1],
	        nfamily = rep(N$dn, 2),
		ncase = NCAS,
		ncontrol = NCOS)
print(parsFDR)
dataFDR <- calculateFDR(pars = parsFDR,
			dnData = allDNData, mutData = allMutData,
			caseData=dataCCcase, controlData=dataCCcontrol,
			geneName = geneName)

## order by increasing q-value
dataFDR <- dataFDR[order(dataFDR$qvalue), ]

head(dataFDR)

dim(dataFDR[dataFDR$qvalue < 0.1, ])
dim(dataFDR[dataFDR$qvalue < 0.05, ])

## write table to file
write.table(dataFDR, file=paste0(OUTROOT, ".gene_res.tsv"),
	    row.names=F, col.names=T,
	    sep="\t", quote=F)
