library(ggplot2)

## PARAM
COLS.TOLINTOL <- c(rgb(0,0.45,0.70,1),rgb(0.80,0.40,0,1))

Main <- function() {

  ARGS <- commandArgs(trailingOnly = T)
  if (length(ARGS) != 4) {
    cat("plot_dnm_rate_caco_analysis.R <rate.csv> <caco.csv> <hc.csv> <outroot>\n")
    q()
  }
  
  rate.csv <- ARGS[1]
  caco.csv <- ARGS[2]
  hc.csv <- ARGS[3]
  outroot <- ARGS[4]

  rate.tbl.all<-read.csv(rate.csv,stringsAsFactors=F)
  caco.tbl.all<-read.csv(caco.csv,stringsAsFactors=F)
  hc.tbl <- read.csv(hc.csv,stringsAsFactors=F)

  # subset on tol/intol
  rate.tbl.exome <- subset(rate.tbl.all, geneset == "exome")
  caco.tbl.exome <- subset(caco.tbl.all, geneset == "exome")
  rate.tbl <- subset(rate.tbl.all, geneset %in% c("tolerant","intolerant") )
  caco.tbl <- subset(caco.tbl.all, geneset %in% c("tolerant","intolerant") )

  # fix pval cols
  PvalFix <- function(pvals) {
    pvals.x <- ifelse(pvals <= 1,
                    paste0("p = ",formatC(pvals,format='e',digits=2)),
                    "")
    print(pvals.x)
    return(pvals.x)
  }
  rate.tbl.exome$ptest_p <- PvalFix(rate.tbl.exome$ptest_p)
  rate.tbl$ptest_p <- PvalFix(rate.tbl$ptest_p) 
  caco.tbl.exome$pois_p <-PvalFix(caco.tbl.exome$pois_p)
  caco.tbl$pois_p <-PvalFix(caco.tbl$pois_p)
  hc.tbl$fet_p <-PvalFix(hc.tbl$fet_p)

  # Figure S1A (exome, rate vs exp(sequence))
  effs <- c("Synon","Mis","MisND","MisD","LoF","Mis/LoF")
  rate.tbl.exome <- subset(rate.tbl.exome, eff %in% effs)
  pd <- position_dodge(width = 0.4)
  rate.tbl.exome$eff <- factor(rate.tbl.exome$eff, levels=effs)
  rate.tbl.exome$geneset <- factor(rate.tbl.exome$geneset,levels=c("exome"))
  rate.tbl.exome$rate_obs <- rate.tbl.exome$rate_obs - rate.tbl.exome$rate_exp
  rate.tbl.exome$rate_95ci_l <- rate.tbl.exome$rate_95ci_l - rate.tbl.exome$rate_exp
  rate.tbl.exome$rate_95ci_u <- rate.tbl.exome$rate_95ci_u - rate.tbl.exome$rate_exp
  gg.rate.exome <- ggplot(rate.tbl.exome, aes(y=rate_obs, x=eff, label=ptest_p))
  gg.rate.exome <- gg.rate.exome + geom_linerange(aes(ymin=rate_95ci_l,
                                               ymax=rate_95ci_u),
                                      position=pd )
  gg.rate.exome <- gg.rate.exome + geom_point(shape = 19, size  = 2, position=pd, shape=19)
  gg.rate.exome <- gg.rate.exome + geom_text(position=position_nudge(x=0.2),
                                             angle=90,size=4)
  gg.rate.exome <- gg.rate.exome + ylab("rate (obs - exp, mutability)")
  gg.rate.exome <- gg.rate.exome + xlab("Annotation")
  gg.rate.exome <- gg.rate.exome + geom_hline(yintercept=0)
  gg.rate.exome <- gg.rate.exome + theme(legend.position = 'none')
  ggsave(paste0(outroot, ".exome.rate.pdf"))

  # Figure S1B (exome, rate vs exp(controls))
  effs <- c("Synon","Mis","MisND","MisD","LoF","Mis/LoF")
  caco.tbl.exome <- subset(caco.tbl.exome, eff %in% effs)
  caco.tbl.exome$eff <- factor(caco.tbl.exome$eff, levels=effs)
  caco.tbl.exome$geneset <- factor(caco.tbl.exome$geneset,
                                   levels=c("tolerant","intolerant"))
  caco.tbl.exome$caco_or_est <- as.numeric(caco.tbl.exome$caco_or_est)
  caco.tbl.exome$ca_rate <- caco.tbl.exome$ca_rate - caco.tbl.exome$co_rate
  caco.tbl.exome$ca_pois_95ci_l <- caco.tbl.exome$ca_pois_95ci_l - caco.tbl.exome$co_rate
  caco.tbl.exome$ca_pois_95ci_u <- caco.tbl.exome$ca_pois_95ci_u - caco.tbl.exome$co_rate
  gg.caco.exome <- ggplot(caco.tbl.exome, aes(x=eff, y=ca_rate,
                                   label=pois_p))
  gg.caco.exome <- gg.caco.exome + geom_hline(yintercept=0,color='black')
  gg.caco.exome <- gg.caco.exome + geom_linerange(aes(ymin=ca_pois_95ci_l,
                                          ymax=ca_pois_95ci_u),
                                      position=pd )
  gg.caco.exome <- gg.caco.exome + geom_point(shape = 19, size  = 2, position=pd)
  gg.caco.exome <- gg.caco.exome + ylab("rate (obs - exp, control trios)")
  gg.caco.exome <- gg.caco.exome + xlab("Annotation")
  gg.caco.exome <- gg.caco.exome + theme(legend.position = 'none')
  gg.caco.exome <- gg.caco.exome + geom_text(angle=90,size=4,
                                             position=position_nudge(x=0.2))
  ggsave(paste0(outroot, ".exome.caco.pdf"))

  # for rest of figures only keep certain annots
  effs <- c("Synon","Mis","LoF","Mis/LoF")
  rate.tbl <- subset(rate.tbl.all, eff %in% effs)
  caco.tbl <- subset(caco.tbl.all, eff %in% effs)
  print(caco.tbl$geneset)
  rate.tbl <- subset(rate.tbl, geneset %in% c("tolerant","intolerant"))
  caco.tbl <- subset(caco.tbl, geneset %in% c("tolerant","intolerant"))

  # Figure 2A
  print(rate.tbl)
  pd <- position_dodge(width = 0.4)
  rate.tbl$eff <- factor(rate.tbl$eff, levels=effs)
  rate.tbl$geneset <- factor(rate.tbl$geneset,
                             levels=c("tolerant","intolerant") )
  rate.tbl$rate_obs <- rate.tbl$rate_obs - rate.tbl$rate_exp
  rate.tbl$rate_95ci_l <- rate.tbl$rate_95ci_l - rate.tbl$rate_exp
  rate.tbl$rate_95ci_u <- rate.tbl$rate_95ci_u - rate.tbl$rate_exp
  rate.tbl$ptest_p <- PvalFix(rate.tbl$ptest_p) 
  gg.rate <- ggplot(rate.tbl, aes(y=rate_obs, x=eff, label=ptest_p, color=geneset)) 
  gg.rate <- gg.rate + geom_linerange(aes(ymin=rate_95ci_l,
                                               ymax=rate_95ci_u), 
                                      position=pd ) 
  gg.rate <- gg.rate + geom_point(shape = 19, size  = 1, position=pd)
  gg.rate <- gg.rate + geom_text(position=position_dodge(0.6), 
                                 angle=90,size=3)
  gg.rate <- gg.rate + ylab("rate (obs - exp, mutability)")
  gg.rate <- gg.rate + xlab("Annotation")
  gg.rate <- gg.rate + scale_color_manual(values=COLS.TOLINTOL)
  gg.rate <- gg.rate + geom_hline(yintercept=0)
  gg.rate <- gg.rate + theme(legend.position = 'none')  
  ggsave(paste0(outroot, ".tolintol.rate.pdf"))

  # Figure 2B
  caco.tbl$eff <- factor(caco.tbl$eff, levels=effs)
  caco.tbl$geneset <- factor(caco.tbl$geneset,
                             levels=c("tolerant","intolerant") )
  caco.tbl$caco_or_est <- as.numeric(caco.tbl$caco_or_est)
  caco.tbl$ca_rate <- caco.tbl$ca_rate - caco.tbl$co_rate
  caco.tbl$ca_pois_95ci_l <- caco.tbl$ca_pois_95ci_l - caco.tbl$co_rate
  caco.tbl$ca_pois_95ci_u <- caco.tbl$ca_pois_95ci_u - caco.tbl$co_rate
  caco.tbl$pois_p <-PvalFix(caco.tbl$pois_p)
  gg.caco <- ggplot(caco.tbl, aes(x=eff, y=ca_rate, 
                                   label=pois_p, fill=geneset, col=geneset))
  gg.caco <- gg.caco + geom_hline(yintercept=0,color='black')
  gg.caco <- gg.caco + geom_linerange(aes(ymin=ca_pois_95ci_l,                   
                                          ymax=ca_pois_95ci_u),             
                                      position=pd )
  gg.caco <- gg.caco + geom_point(shape = 19, size  = 1, position=pd)
  gg.caco <- gg.caco + scale_color_manual(values=COLS.TOLINTOL)
  gg.caco <- gg.caco + ylab("rate (obs - exp, control trios)")
  gg.caco <- gg.caco + xlab("Annotation")
  gg.caco <- gg.caco + theme(legend.position = 'none')  
  # gg.caco <- gg.caco + ylim(-0.051,0.121)
  gg.caco <- gg.caco + geom_text(position=position_dodge(width=0.6),
                                   angle=90,size=3)
  # position_dodge(width=0.25)
  ggsave(paste0(outroot, ".tolintol.caco.pdf"))

  # Figure 2C
  pd <- position_dodge(width = 0.9)
  # hc.tbl <- subset(hc.tbl, Annotation != "Synon")
  hc.tbl$or_est <- as.numeric(hc.tbl$or_est)
  hc.tbl$Annotation <- factor(hc.tbl$Annotation, 
                              levels=rev(hc.tbl$Annotation))
  gg.hc <- ggplot(hc.tbl, aes(x=Annotation, y=or_est, label=fet_p))
  gg.hc <- gg.hc + geom_hline(yintercept=1,color='black')
  gg.hc <- gg.hc + geom_linerange(aes(ymin=or_95ci_l,                   
                                      ymax=or_95ci_u),             
                                      position=pd,
                                      color=COLS.TOLINTOL[2])
  gg.hc <- gg.hc + geom_point(shape = 19, size  = 1, position=pd, 
                              color=COLS.TOLINTOL[2])
  gg.hc <- gg.hc + theme(text=element_text(size=12))
  gg.hc <- gg.hc + theme(axis.text.y=element_text(angle=0,hjust=1)) 
  gg.hc <- gg.hc +  theme(legend.position = 'none')
  gg.hc <- gg.hc + geom_text(color=COLS.TOLINTOL[2], 
                             position=position_nudge(x=0.3))
  gg.hc <- gg.hc + coord_flip()
  gg.hc <- gg.hc + ylab("case/control odds ratio")
  gg.hc <- gg.hc + xlab("Annotation")
  ggsave(paste0(outroot, ".hc.caco.pdf"))

  pdf(paste0(outroot, ".rate_caco_hc.pdf"))
  multiplot(gg.rate, gg.caco, gg.hc, cols=1)
  dev.off()

  pdf(paste0(outroot, ".exome.rate_caco.pdf"))                                  
  multiplot(gg.rate.exome, gg.caco.exome, cols=1)                               
  dev.off() 

  # make a table focused on polyphen partitioned missense across exome:
  # grab subset of vars with polyphen annots
  caco.tbl <- subset(caco.tbl.all, eff %in% effs)
  effs.pphn <- c("MisB","MisP","MisD")
  rate.tbl.pphn <- subset(rate.tbl.all, eff %in% effs.pphn)
  caco.tbl.pphn <- subset(caco.tbl.all, eff %in% effs.pphn)

  # subset on exome
  rate.tbl.pphn <- subset(rate.tbl.pphn, geneset=="exome")
  caco.tbl.pphn <- subset(caco.tbl.pphn, geneset=="exome")

  # pval fix
  rate.tbl.pphn$ptest_p <- PvalFix(rate.tbl.pphn$ptest_p) 
  caco.tbl.pphn$pois_p <-PvalFix(caco.tbl.pphn$pois_p)
   
  # make rate section of figure
  pd <- position_dodge(width = 0.4)
  rate.tbl.pphn$eff <- factor(rate.tbl.pphn$eff, levels=effs.pphn)
  rate.tbl.pphn$rate_obs <- rate.tbl.pphn$rate_obs - rate.tbl.pphn$rate_exp
  rate.tbl.pphn$rate_95ci_l <- rate.tbl.pphn$rate_95ci_l - rate.tbl.pphn$rate_exp
  rate.tbl.pphn$rate_95ci_u <- rate.tbl.pphn$rate_95ci_u - rate.tbl.pphn$rate_exp
  gg.rate <- ggplot(rate.tbl.pphn, aes(y=rate_obs, x=eff, label=ptest_p)) 
  gg.rate <- gg.rate + geom_linerange(aes(ymin=rate_95ci_l,
                                          ymax=rate_95ci_u), 
                                      position=pd ) 
  gg.rate <- gg.rate + geom_point(shape = 19, size  = 1, position=pd)
  gg.rate <- gg.rate + geom_text(nudge_x=0.1, 
                                 angle=90,size=3)
  gg.rate <- gg.rate + ylab("rate (obs - exp, mutability)")
  gg.rate <- gg.rate + xlab("Annotation")
  gg.rate <- gg.rate + geom_hline(yintercept=0)
  gg.rate <- gg.rate + theme(legend.position = 'none')  
  ggsave(paste0(outroot, ".exome_mis_polyphen2.rate.pdf"))

  # Figure 2B
  caco.tbl.pphn$eff <- factor(caco.tbl.pphn$eff, levels=effs.pphn)
  caco.tbl.pphn$caco_or_est <- as.numeric(caco.tbl.pphn$caco_or_est)
  caco.tbl.pphn$ca_rate <- caco.tbl.pphn$ca_rate - caco.tbl.pphn$co_rate
  caco.tbl.pphn$ca_pois_95ci_l <- caco.tbl.pphn$ca_pois_95ci_l - caco.tbl.pphn$co_rate
  caco.tbl.pphn$ca_pois_95ci_u <- caco.tbl.pphn$ca_pois_95ci_u - caco.tbl.pphn$co_rate
  gg.caco <- ggplot(caco.tbl.pphn, aes(x=eff, y=ca_rate, 
                                       label=pois_p))
  gg.caco <- gg.caco + geom_hline(yintercept=0,color='black')
  gg.caco <- gg.caco + geom_linerange(aes(ymin=ca_pois_95ci_l,                   
                                          ymax=ca_pois_95ci_u),             
                                      position=pd )
  gg.caco <- gg.caco + geom_point(shape = 19, size  = 1, position=pd)
  gg.caco <- gg.caco + ylab("rate (obs - exp, control trios)")
  gg.caco <- gg.caco + xlab("Annotation")
  gg.caco <- gg.caco + theme(legend.position = 'none')  
  gg.caco <- gg.caco + ylim(-0.051,0.121)
  gg.caco <- gg.caco + geom_text(nudge_x=0.1,
                                 angle=90,size=3)
  ggsave(paste0(outroot, ".exome_mis_polyphen2.caco.pdf"))

  # multiplot for polyphen classifs
  pdf(paste0(outroot, ".exome.mis_polyphen2.rate_caco.pdf"))                                  
  multiplot(gg.rate, gg.caco, cols=1)                               
  dev.off() 

  # make a table focused on polyphen partitioned missense across exome:
  # grab subset of vars with polyphen annots
  # make a table focused on LoF variation in LOEUF deciles
  caco.tbl <- subset(caco.tbl.all, eff %in% ("LoF"))
  rate.tbl.lof <- subset(rate.tbl.all, eff %in% ("LoF"))
  rate.tbl.lof <- subset(rate.tbl.lof, grepl("LOEUF_",geneset))
  caco.tbl.lof <- subset(caco.tbl.all, eff %in% ("LoF"))
  caco.tbl.lof <- subset(caco.tbl.lof, grepl("LOEUF_",geneset))

  # fix geneset column
  rate.tbl.lof$geneset <- gsub("LOEUF_","",rate.tbl.lof$geneset)
  caco.tbl.lof$geneset <- gsub("LOEUF_","",caco.tbl.lof$geneset) 

  # pval fix
  rate.tbl.lof$ptest_p <- PvalFix(rate.tbl.lof$ptest_p) 
  caco.tbl.lof$pois_p <-PvalFix(caco.tbl.lof$pois_p)

  # add column for LOEUF decile
  rate.tbl.lof$LOEUF_decile <- factor(rate.tbl.lof$geneset,
                                      levels=1:10)
  caco.tbl.lof$LOEUF_decile <- factor(caco.tbl.lof$geneset,
                                      levels=1:10)

  # make rate section of figure
  pd <- position_dodge(width = 0.4)
  rate.tbl.lof$rate_obs <- rate.tbl.lof$rate_obs - rate.tbl.lof$rate_exp
  rate.tbl.lof$rate_95ci_l <- rate.tbl.lof$rate_95ci_l - rate.tbl.lof$rate_exp
  rate.tbl.lof$rate_95ci_u <- rate.tbl.lof$rate_95ci_u - rate.tbl.lof$rate_exp
  gg.rate <- ggplot(rate.tbl.lof, aes(y=rate_obs, x=LOEUF_decile, label=ptest_p)) 
  gg.rate <- gg.rate + geom_linerange(aes(ymin=rate_95ci_l,
                                          ymax=rate_95ci_u), 
                                      position=pd ) 
  gg.rate <- gg.rate + geom_point(shape = 19, size  = 1, position=pd)
  gg.rate <- gg.rate + geom_text(nudge_x=0.2, 
                                 angle=90,size=3)
  gg.rate <- gg.rate + ylab("rate (obs - exp, mutability)")
  gg.rate <- gg.rate + xlab("LOEUF decile")
  gg.rate <- gg.rate + geom_hline(yintercept=0)
  gg.rate <- gg.rate + theme(legend.position = 'none')  
  ggsave(paste0(outroot, ".exome_lof_LOEUF_deciles.rate.pdf"))

  # Figure 2B
  caco.tbl.lof$caco_or_est <- as.numeric(caco.tbl.lof$caco_or_est)
  caco.tbl.lof$ca_rate <- caco.tbl.lof$ca_rate - caco.tbl.lof$co_rate
  caco.tbl.lof$ca_pois_95ci_l <- caco.tbl.lof$ca_pois_95ci_l - caco.tbl.lof$co_rate
  caco.tbl.lof$ca_pois_95ci_u <- caco.tbl.lof$ca_pois_95ci_u - caco.tbl.lof$co_rate
  gg.caco <- ggplot(caco.tbl.lof, aes(x=LOEUF_decile, y=ca_rate, 
                                       label=pois_p))
  gg.caco <- gg.caco + geom_hline(yintercept=0,color='black')
  gg.caco <- gg.caco + geom_linerange(aes(ymin=ca_pois_95ci_l,                   
                                          ymax=ca_pois_95ci_u),             
                                      position=pd )
  gg.caco <- gg.caco + geom_point(shape = 19, size  = 1, position=pd)
  gg.caco <- gg.caco + ylab("rate (obs - exp, control trios)")
  gg.caco <- gg.caco + xlab("LOEUF decile")
  gg.caco <- gg.caco + theme(legend.position = 'none')  
  gg.caco <- gg.caco + geom_text(nudge_x=0.2,
                                 angle=90,size=3)
  ggsave(paste0(outroot, ".exome_lof_LOEUF_deciles.caco.pdf"))

  # multiplot for polyphen classifs
  pdf(paste0(outroot, ".exome.lof_LOEUF_deciles.rate_caco.pdf"))                                  
  multiplot(gg.rate, gg.caco, cols=1)                               
  dev.off() 


}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
  
PvalFix <- function(pvals) {
  pvals.x <- ifelse(pvals <= 1,
                    paste0("p = ",formatC(pvals,format='e',digits=2)),
                    "")
  return(pvals.x)
}

if (interactive() == F) {
  Main()
}
