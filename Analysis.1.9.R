##########################################################
### This script is used to analyze and produce figures ###
### in the yeast mutagenesis project                   ###
###                                                    ###
### Date: 2025-02-27                                   ###
### Version: 1.9                                       ###
### Author: BY and EWM                                 ###
##########################################################

#######################################
#####*****STEP 0: Preparation*****#####
#######################################
###***0.1. Set working directory and clean working space***###
rm(list=ls())
options(warn=-1)


parent.dir <- "x"
output.dir <- "x"

###***0.2. Load libraries***###
library(ggplot2)
library(ggjoy)
library(RColorBrewer)
library(robustbase)
library(rstatix)
library(dplyr)

###***0.3. Necessary Functions***###

## permutation test for mutational variance ##
permutation.test <- function(pop, draw=300,iter=10000) {
  r <- c()
  for(i in 1:iter) {
    s.pop <- sample(pop,draw,replace=TRUE)
    r <- c(r, var(s.pop))
  }
  return(r)
}

## permutation test for medcouple ##
permutation.test.mc <- function(pop, draw=300,iter=10000) {
  s <- c()
  for(i in 1:iter) {
    s.pop <- sample(pop,draw,replace=TRUE)
    s <- c(s, mc(s.pop))
  }
  return(s)
}

#########################################
#####*****STEP 1: Load data sets*****#####
#########################################
load.types <- c(rep("factor",4), rep("numeric",18))
setwd(parent.dir)

DATA.1.LOAD <- read.table("STRAIN.ESTIMATES.1.txt",sep="\t",header=TRUE,colClasses=load.types)
DATA.2.LOAD <- read.table("STRAIN.ESTIMATES.2.txt",sep="\t",header=TRUE,colClasses=load.types)
DATA.3.LOAD <- read.table("STRAIN.ESTIMATES.3.txt",sep="\t",header=TRUE,colClasses=load.types)
DATA.4.LOAD <- read.table("STRAIN.ESTIMATES.4.txt",sep="\t",header=TRUE,colClasses=load.types)
DATA.5.LOAD <- read.table("STRAIN.ESTIMATES.5.txt",sep="\t",header=TRUE,colClasses=load.types)
DATA.6.LOAD <- read.table("STRAIN.ESTIMATES.6.NEW.txt",sep="\t",header=TRUE,colClasses=load.types)

DATA.0.LOAD <- read.table("SUMMARY.TRANS.txt",sep="\t",header=TRUE)

###***1.1. Rename ALL data columns: DATA.1 to DATA.0***###
DATA.0 <- subset(DATA.0.LOAD, CLASS!="NULL", select=c("CLASS",
	"YFP.MEDIAN.RELATIVE.MEAN","YFP.MEDIAN.RELATIVE.SD",
	"YFP.CV.RELATIVE.MEAN","YFP.CV.RELATIVE.SD",
	"YFP.SD.RELATIVE.MEAN","YFP.SD.RELATIVE.SD",
	"P.VAL.MEDIAN","P.VAL.NOISE"))

DATA.0.ALTER <- subset(DATA.0.LOAD, CLASS!="NULL", select=c("CLASS",
	"YFP.MEDIAN.RELATIVE.MEAN","YFP.MEDIAN.RELATIVE.SD",
	"YFP.CV.RELATIVE.MEAN","YFP.CV.RELATIVE.SD"))

selected.cols <- c("STRAIN",
	"YFP.MEDIAN.R2OWN.MEAN","YFP.MEDIAN.R2OWN.SD",
	"YFP.CV.R2OWN.MEAN","YFP.CV.R2OWN.SD",
	"YFP.SD.R2OWN.MEAN","YFP.SD.R2OWN.SD",
	"pvalues.median.r2own","pvalues.cv.r2own")
selected.cols.alter <- c("STRAIN",
	"YFP.MEDIAN.R2WT.MEAN","YFP.MEDIAN.R2WT.SD",
	"YFP.CV.R2WT.MEAN","YFP.CV.R2WT.SD")

new.names <- c("strain",
	"median.mean","median.sd",
	"cv.mean","cv.sd",
	"sd.mean","sd.sd",
	"pvalues.median","pvalues.noise")
new.names.alter <- c("strain",
	"median.mean","median.sd",
	"cv.mean","cv.sd")

colnames(DATA.0) <- new.names
DATA.1 <- subset(DATA.1.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.1) <- new.names
DATA.2 <- subset(DATA.2.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.2) <- new.names
DATA.3 <- subset(DATA.3.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.3) <- new.names
DATA.4 <- subset(DATA.4.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.4) <- new.names
DATA.5 <- subset(DATA.5.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.5) <- new.names
DATA.6 <- subset(DATA.6.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.6) <- new.names

colnames(DATA.0.ALTER) <- new.names.alter
DATA.1.ALTER <- subset(DATA.1.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols.alter); colnames(DATA.1.ALTER) <- new.names.alter
DATA.2.ALTER <- subset(DATA.2.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols.alter); colnames(DATA.2.ALTER) <- new.names.alter
DATA.3.ALTER <- subset(DATA.3.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols.alter); colnames(DATA.3.ALTER) <- new.names.alter
DATA.4.ALTER <- subset(DATA.4.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols.alter); colnames(DATA.4.ALTER) <- new.names.alter
DATA.5.ALTER <- subset(DATA.5.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols.alter); colnames(DATA.5.ALTER) <- new.names.alter
DATA.6.ALTER <- subset(DATA.6.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols.alter); colnames(DATA.6.ALTER) <- new.names.alter

###***1.2. Label different genotypes***###
DATA.0$genotype <- "wt.1"
DATA.0.ALTER$genotype <- "wt.1"

generate.genotypes <- function(pairs, strains) {
	### pairs: names of the genotypes
	### strains: the labels of the strains in each experiment
	result <- c()
	for(s in strains) {
		if(s %in% c("ems.m76","sham.m76")) { result <- c(result,pairs[1]) }
		else if (s %in% c("ems.tata","sham.tata")) { result <- c(result,pairs[2]) }
		else { result <- c(result,"BY") }
	}
	return(result)
}

DATA.1$genotype <- generate.genotypes(c("m76","tata1"), DATA.1$strain)
DATA.2$genotype <- generate.genotypes(c("m66","tata2"), DATA.2$strain)
DATA.3$genotype <- generate.genotypes(c("tye7","ade6"), DATA.3$strain)
DATA.4$genotype <- generate.genotypes(c("nam7","rap1"), DATA.4$strain)
DATA.5$genotype <- generate.genotypes(c("m22","yps1000"), DATA.5$strain)
DATA.6$genotype <- generate.genotypes(c("wt.2","wt.3"), DATA.6$strain)

DATA.1.ALTER$genotype <- generate.genotypes(c("m76","tata1"), DATA.1.ALTER$strain)
DATA.2.ALTER$genotype <- generate.genotypes(c("m66","tata2"), DATA.2.ALTER$strain)
DATA.3.ALTER$genotype <- generate.genotypes(c("tye7","ade6"), DATA.3.ALTER$strain)
DATA.4.ALTER$genotype <- generate.genotypes(c("nam7","rap1"), DATA.4.ALTER$strain)
DATA.5.ALTER$genotype <- generate.genotypes(c("m22","yps1000"), DATA.5.ALTER$strain)
DATA.6.ALTER$genotype <- generate.genotypes(c("wt.2","wt.3"), DATA.6.ALTER$strain)

###***1.3. Assemble the datasets***###
expids <- c(
	rep(0, nrow(DATA.0)),
	rep(1, nrow(DATA.1)),
	rep(2, nrow(DATA.2)),
	rep(3, nrow(DATA.3)),
	rep(4, nrow(DATA.4)),
	rep(5, nrow(DATA.5)),
	rep(6, nrow(DATA.6)))

DATA.ALL <- rbind(DATA.0, DATA.1, DATA.2, DATA.3, DATA.4, DATA.5, DATA.6)
DATA.ALL$expids <- expids
DATA.ALL.ALTER <- rbind(DATA.0.ALTER, DATA.1.ALTER, DATA.2.ALTER, DATA.3.ALTER, DATA.4.ALTER, DATA.5.ALTER, DATA.6.ALTER)
DATA.ALL.ALTER$expids <- expids

treatment <- c()
for(t in DATA.ALL$strain) {
	if(t %in% c("EMS","ems.m76","ems.tata")) { treatment <- c(treatment,"ems") }
	else if(t %in% c("SHAM","sham.wt","sham.m76","sham.tata")) { treatment <- c(treatment,"sham") }
}

DATA.ALL$treatment <- treatment
DATA.ALL.ALTER$treatment <- treatment

DATA.ALL$genotype <- as.factor(DATA.ALL$genotype)
DATA.ALL$expids <- as.factor(DATA.ALL$expids)
DATA.ALL$treatment <- as.factor(DATA.ALL$treatment)

DATA.ALL.ALTER$genotype <- as.factor(DATA.ALL.ALTER$genotype)
DATA.ALL.ALTER$expids <- as.factor(DATA.ALL.ALTER$expids)
DATA.ALL.ALTER$treatment <- as.factor(DATA.ALL.ALTER$treatment)

DATA.ALL <- droplevels(DATA.ALL)
DATA.ALL.ALTER <- droplevels(DATA.ALL.ALTER)

### make table of isolate counts for each genotype ###

df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("genotype", "sham", "ems"))

gen<-as.vector(unique(DATA.ALL.ALTER$genotype))

for (i in 1:length(unique(DATA.ALL.ALTER$genotype))){
  df[i,1]<-gen[i]
  df[i,2]<-sum(DATA.ALL.ALTER$genotype == gen[i] & DATA.ALL.ALTER$treatment == "sham", na.rm=TRUE)
  df[i,3]<-sum(DATA.ALL.ALTER$genotype == gen[i] & DATA.ALL.ALTER$treatment == "ems", na.rm=TRUE)    
}

treatment_subset<-df[c(2,4:10),]

sham_mean<-mean(treatment_subset$sham)
ems_mean<-mean(treatment_subset$ems)

treatment_subset[9,1]<-"avg"
treatment_subset[9,2]<-sham_mean
treatment_subset[9,3]<-ems_mean

###***1.4. Create z-scores statistics***###
A <- aggregate(median.mean~genotype, data=subset(DATA.ALL.ALTER,
  treatment=="sham"), FUN=mean)
B <- aggregate(median.mean~genotype, data=subset(DATA.ALL.ALTER,
  treatment=="sham"), FUN=sd)
standard.table <- merge(A, B, by="genotype")
colnames(standard.table) <- c("genotype","mu","sigma")

DATA.ALL.ALTER$zscore <- 0
for(i in 1:nrow(DATA.ALL.ALTER)) {
  g <- DATA.ALL.ALTER[i,"genotype"]
  mu <- subset(standard.table,genotype==g)[,"mu"]
  sigma <- subset(standard.table,genotype==g)[,"sigma"]
  DATA.ALL.ALTER[i,"zscore"] <- (DATA.ALL.ALTER[i,"median.mean"]-mu)/sigma
}


#############################################
#####*****Step 2: Reference Strains*****#####
#############################################

# wt.1 is equivalent to ref.0, and refers to the Metzger et al. 2016 data set 
# wt.2 is equivalent to ref.1
# wt.3 is equivalent to ref.2
#yps100 is an alternate wt genotype that was not included in the manuscript

###***2.1. Build the ggplot data structure***###
DATA.FIGURE.2 <- subset(DATA.ALL.ALTER, genotype %in% c("wt.1","wt.2","wt.3","yps1000"))
DATA.FIGURE.2 <- droplevels(DATA.FIGURE.2)
plot.colors <- brewer.pal(4, 'Set3')

### Compare SHAM distributions
fig <- ggplot(subset(DATA.FIGURE.2,treatment=="sham"), aes(x=median.mean,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

A <- subset(DATA.FIGURE.2,treatment=="sham")
alter.genotypes <- c("wt.1","wt.2","wt.3","yps1000")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"median.mean"]
  print(mean(result))
}

str(DATA.FIGURE.2)

sham <-subset(DATA.FIGURE.2,treatment=="sham")
str(sham)

boxplot(median.mean ~ genotype,
        data = sham
)

pw <- sham %>% 
  pairwise_t_test(
    median.mean ~ genotype, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

boxplot(median.sd ~ genotype,
        data = sham
)

pw <- sham %>% 
  pairwise_t_test(
    median.sd ~ genotype, pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

### Compare EMS distributions
fig <- ggplot(subset(DATA.FIGURE.2,treatment=="ems"), aes(x=median.mean,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig

B <- subset(DATA.FIGURE.2,treatment=="ems")
alter.genotypes <- c("wt.1","wt.2","wt.3","yps1000")

for(i in alter.genotypes) {
  result <- subset(B, genotype==i)[,"median.mean"]
  print(mean(result))
}

# t-tests of sham to ems for each genotype

for (i in 1:length(alter.genotypes)){
    
   y <- subset(DATA.FIGURE.2,genotype==alter.genotypes[i])
   yt <- t_test(y, median.mean~treatment, p.adjust.method = "bonferroni") 
   print(alter.genotypes[i])
   print(yt)
  }


for (i in 1:length(alter.genotypes)){
  
  y <- subset(DATA.FIGURE.2,genotype==alter.genotypes[i])
  yt <- t_test(y, median.sd~treatment, p.adjust.method = "bonferroni") 
  print(alter.genotypes[i])
  print(yt)
}

ems <-subset(DATA.FIGURE.2,treatment=="ems")

pw <- ems %>% 
  pairwise_t_test(
    median.mean ~ genotype,
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

pw <- ems %>% 
  pairwise_t_test(
    median.sd ~ genotype,
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

### Compare zscore distributions
fig <- ggplot(subset(DATA.FIGURE.2,treatment=="ems"), aes(x=zscore,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig

A <- subset(DATA.FIGURE.2,treatment=="ems")
alter.genotypes <- c("wt.1","wt.2","wt.3","yps1000")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(mean(result))
}

ems <-subset(DATA.FIGURE.2,treatment=="ems")

boxplot(zscore ~ genotype,
        data = ems
)


pw <- ems %>% 
  pairwise_t_test(
    zscore ~ genotype, 
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

boxplot(median.sd ~ genotype,
        data = sham
)


### Compare YPS1000 and wt replicates mv to wt.1 (ref.0)

raw.dist.mv <- permutation.test(subset(DATA.ALL.ALTER,genotype=="wt.1"&treatment=="ems")[,"zscore"])

raw.dist.mv.sorted <-sort(raw.dist.mv) #must be sorted to calculate intervals

#write.csv(raw.dist.mv.sorted, file = "x.csv") # (save in whatever format you like)

# If you plan to use one reference mv distribution for all genotypes (as opposed to pulling 
# new distributions each time), run the above once and then do not overwrite these 
# objects for subsequent statistics and plotting.

raw.dist.mv.sorted <-read.csv("raw.dist.mv.sorted.csv", header = FALSE) # to pull from the included file
raw.dist.mv.sorted <-raw.dist.mv.sorted[,1]

int.lower<-raw.dist.mv.sorted[250] #for 2 tailed
int.upper<-raw.dist.mv.sorted[9750] #for 2 tailed

int.lower<-raw.dist.mv.sorted[500] #for 1 tailed
int.upper<-raw.dist.mv.sorted[9500] #for 1 tailed

# For figures, if all background plots should be identical (using the same 
# distribution), then the following plot is the base, with ablines to  
# be added for each section:

# base wt.1 distribution and intervals:

fig <- ggplot(data.frame(mv=raw.dist.mv.sorted),aes(mv))+
  geom_histogram(binwidth=0.1, alpha = 0.2)+
  theme_bw() 
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + geom_vline(xintercept=min(raw.dist.mv.sorted), color="black")
fig <- fig + geom_vline(xintercept=max(raw.dist.mv.sorted), color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))
fig <- fig + 
  xlim(0, 18)

fig

# add wt strain and reference replicate comparisons:

alter.genotypes <- c("yps1000","wt.2","wt.3")

A <- aggregate(zscore~genotype, data=subset(DATA.FIGURE.2,treatment=="ems"),
  FUN=var)
for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(length(raw.dist.mv.sorted[which(raw.dist.mv.sorted<result)])/length(raw.dist.mv.sorted))
}

fig <- ggplot(data.frame(mv=raw.dist.mv.sorted),aes(mv))
fig <- fig + geom_histogram(binwidth=0.1, alpha=0.2)
fig <- fig + geom_vline(xintercept=A[2,2], color=plot.colors[4])
fig <- fig + geom_vline(xintercept=A[3,2], color=plot.colors[1])
fig <- fig + geom_vline(xintercept=A[4,2], color=plot.colors[3])
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + annotate("text",x=3.36,y=200,label="YPS1000, 3.36")
fig <- fig + annotate("text",x=5.61,y=250,label="wt.2, 5.61")
fig <- fig + annotate("text",x=4.28,y=300,label="wt.3, 4.28")
fig <- fig + xlim(0, 18)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig
A

### Compare YPS1000 and wt replicates mc to wt.1 (ref.0)

raw.dist.mc <- permutation.test.mc(subset(DATA.ALL.ALTER,genotype=="wt.1"&treatment=="ems")[,"zscore"])

raw.dist.mc.sorted <-sort(raw.dist.mc) #for intervals

#write.csv(raw.dist.mc.sorted, file = "x.csv") # (save in whatever format you like)

# If you plan to use one reference mc distribution for all genotypes (as opposed to pulling 
# new distributions each time), run the above once and then do not overwrite these 
# objects for subsequent statistics and plotting.

raw.dist.mc.sorted <-read.csv("raw.dist.mc.sorted.csv", header = FALSE) # to pull from the included file
raw.dist.mc.sorted <-raw.dist.mc.sorted[,1]

int.lower<-raw.dist.mc.sorted[250] #for 2 tailed
int.upper<-raw.dist.mc.sorted[9750] #for 2 tailed

int.lower<-raw.dist.mc.sorted[500] #for 1 tailed
int.upper<-raw.dist.mc.sorted[9500] #for 1 tailed

# For figures, if all background plots should be identical (using the same 
# distribution), then the following plot is the base, with ablines to  
# be added for each section:

# base wt.1 distribution and intervals:

fig <- ggplot(data.frame(mc=raw.dist.mc.sorted),aes(mc))+
  geom_histogram(binwidth=0.002, alpha=0.2)+
  theme_bw() 
fig <- fig + geom_vline(xintercept=int.lower.mc, color="black")
fig <- fig + geom_vline(xintercept=int.upper.mc, color="black")
fig <- fig + geom_vline(xintercept=min(raw.dist.mc.sorted), color="black")
fig <- fig + geom_vline(xintercept=max(raw.dist.mc.sorted), color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))
fig <- fig + 
  xlim(-0.4, 0.3)

fig

# add wt strain and reference replicate comparisons:

alter.genotypes <- c("yps1000","wt.2","wt.3")

A <- aggregate(zscore~genotype, data=subset(DATA.FIGURE.2,treatment=="ems"),
               FUN=mc)
for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(length(raw.dist.mc.sorted[which(raw.dist.mc.sorted<result)])/length(raw.dist.mc.sorted))
}

fig <- ggplot(data.frame(mc=raw.dist.mc.sorted),aes(mc))
fig <- fig + geom_histogram(binwidth=0.002, alpha=0.2)
fig <- fig + geom_vline(xintercept=A[1,2], color="red")
fig <- fig + geom_vline(xintercept=A[2,2], color="blue")
fig <- fig + geom_vline(xintercept=A[3,2], color="purple")
fig <- fig + geom_vline(xintercept=int.lower.mc, color="black")
fig <- fig + geom_vline(xintercept=int.upper.mc, color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))

fig
A

#######################################
#####*****Step 3: Cis mutants*****#####
#######################################

#m66 = rap1 binding site
#m76 = gcr1 binding site
#wt3 arbitrarily chosen (coin flip) for comparisons and figures

###***3.1. Build the ggplot data structure***###
DATA.FIGURE.3 <- subset(DATA.ALL.ALTER, genotype %in% c("wt.3","m76","tata1","m66","tata2"))
DATA.FIGURE.3 <- droplevels(DATA.FIGURE.3)
plot.colors <- brewer.pal(5, 'Set3')

### Compare SHAM distributions
fig <- ggplot(subset(DATA.FIGURE.3,treatment=="sham"), aes(x=median.mean,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig <- fig + xlim(0, 1.20)

fig

A <- subset(DATA.FIGURE.3,treatment=="sham")

alter.genotypes <- c("wt.3", "m76","tata1","m66","tata2")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"median.mean"]
  print(mean(result))
}

sham <-subset(DATA.FIGURE.3,treatment=="sham")

boxplot(median.mean ~ genotype,
        data = sham
)

pw <- sham %>% 
  pairwise_t_test(
    median.mean ~ genotype, 
    comparisons = list(c("wt.3", "m76"), c("wt.3", "tata1"), c("wt.3", "m66"), c("wt.3", "tata2")),
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

boxplot(median.sd ~ genotype,
        data = sham
)

pw <- sham %>% 
  pairwise_t_test(
    median.sd ~ genotype, 
    comparisons = list(c("wt.3", "m76"), c("wt.3", "tata1"), c("wt.3", "m66"), c("wt.3", "tata2")),
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

### Compare EMS distributions
fig <- ggplot(subset(DATA.FIGURE.3,treatment=="ems"), aes(x=median.mean,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig <- fig + xlim(0, 1.20)

fig

A <- subset(DATA.FIGURE.3,treatment=="ems")

alter.genotypes <- c("wt.3", "m76","tata1","m66","tata2")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"median.mean"]
  print(mean(result))
}

pw <- ems %>% 
  pairwise_t_test(
    median.mean ~ genotype, 
    comparisons = list(c("wt.3", "m76"), c("wt.3", "tata1"), c("wt.3", "m66"), c("wt.3", "tata2")),
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

boxplot(median.sd ~ genotype,
        data = sham
)

pw <- ems %>% 
  pairwise_t_test(
    median.sd ~ genotype, 
    comparisons = list(c("wt.3", "m76"), c("wt.3", "tata1"), c("wt.3", "m66"), c("wt.3", "tata2")),
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

### Compare zscore distributions

fig <- ggplot(subset(DATA.FIGURE.3,treatment=="ems"), aes(x=zscore,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig

A <- subset(DATA.FIGURE.3,treatment=="ems")
alter.genotypes <- c("wt.3", "m76","tata1","m66","tata2")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(mean(result))
}

ems <-subset(DATA.FIGURE.3,treatment=="ems")

boxplot(zscore ~ genotype,
        data = ems
)

pw <- ems %>% 
  pairwise_t_test(
    zscore ~ genotype, 
    comparisons = list(c("wt.3", "m76"), c("wt.3", "tata1"), c("wt.3", "m66"), c("wt.3", "tata2")),
    pool.sd = FALSE,
    p.adjust.method = "bonferroni"
  )
pw

# t-tests of sham to ems for each genotype

str(DATA.FIGURE.3)

for (i in 1:length(alter.genotypes)){
  
  y <- subset(DATA.FIGURE.3,genotype==alter.genotypes[i])
  yt <- t_test(y, median.mean~treatment, p.adjust.method = "bonferroni") 
  print(alter.genotypes[i])
  print(yt)
}


for (i in 1:length(alter.genotypes)){
  
  y <- subset(DATA.FIGURE.3,genotype==alter.genotypes[i])
  yt <- t_test(y, median.sd~treatment, p.adjust.method = "bonferroni") 
  print(alter.genotypes[i])
  print(yt)
}


### Compare cis mutant mv to wt.1 (ref.0)

raw.dist.mv <- permutation.test(subset(DATA.ALL.ALTER,genotype=="wt.1"&treatment=="ems")[,"zscore"])

raw.dist.mv.sorted <-sort(raw.dist.mv) #must be sorted to calculate intervals

#write.csv(raw.dist.mv.sorted, file = "x.csv") # (save in whatever format you like)

# If you plan to use one reference mv distribution for all genotypes (as opposed to pulling 
# new distributions each time), run the above once and then do not overwrite these 
# objects for subsequent statistics and plotting.

raw.dist.mv.sorted <-read.csv("raw.dist.mv.sorted.csv", header = FALSE) # to pull from the included file
raw.dist.mv.sorted <-raw.dist.mv.sorted[,1]

int.lower<-raw.dist.mv.sorted[250] #for 2 tailed
int.upper<-raw.dist.mv.sorted[9750] #for 2 tailed

int.lower<-raw.dist.mv.sorted[500] #for 1 tailed
int.upper<-raw.dist.mv.sorted[9500] #for 1 tailed

# For figures, if all background plots should be identical (using the same 
# distribution), then the following plot is the base, with ablines to  
# be added for each section:

# base wt.1 distribution and intervals:

fig <- ggplot(data.frame(mv=raw.dist.mv.sorted),aes(mv))+
  geom_histogram(binwidth=0.1, alpha = 0.2)+
  theme_bw() 
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + geom_vline(xintercept=min(raw.dist.mv.sorted), color="black")
fig <- fig + geom_vline(xintercept=max(raw.dist.mv.sorted), color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))
fig <- fig + 
  xlim(0, 18)

fig


# add cis mutant strain comparisons:

alter.genotypes <- c("m76","tata1","m66","tata2")

A <- aggregate(zscore~genotype, data=subset(DATA.FIGURE.3,treatment=="ems"),
  FUN=var)
for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(length(raw.dist.mv.sorted[which(raw.dist.mv.sorted<result)])/length(raw.dist.mv.sorted))
}

fig <- ggplot(data.frame(mv=raw.dist.mv.sorted),aes(mv))
fig <- fig + geom_histogram(binwidth=0.1, alpha=0.2)
fig <- fig + geom_vline(xintercept=A[1,2], color="red")
fig <- fig + geom_vline(xintercept=A[2,2], color="blue")
fig <- fig + geom_vline(xintercept=A[3,2], color="purple")
fig <- fig + geom_vline(xintercept=A[4,2], color="green")
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + theme_classic()
fig <- fig + xlim(0, 18)
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig
A

### Compare cis mutant mc to wt.1 (ref.0)

raw.dist.mc <- permutation.test.mc(subset(DATA.ALL.ALTER,genotype=="wt.1"&treatment=="ems")[,"zscore"])

raw.dist.mc.sorted <-sort(raw.dist.mc) #for intervals

#write.csv(raw.dist.mc.sorted, file = "x.csv") # (save in whatever format you like)

# If you plan to use one reference mc distribution for all genotypes (as opposed to pulling 
# new distributions each time), run the above once and then do not overwrite these 
# objects for subsequent statistics and plotting.

raw.dist.mc.sorted <-read.csv("raw.dist.mc.sorted.csv", header = FALSE) # to pull from the included file
raw.dist.mc.sorted <-raw.dist.mc.sorted[,1]

int.lower<-raw.dist.mc.sorted[250] #for 2 tailed
int.upper<-raw.dist.mc.sorted[9750] #for 2 tailed

int.lower<-raw.dist.mc.sorted[500] #for 1 tailed
int.upper<-raw.dist.mc.sorted[9500] #for 1 tailed

# For figures, if all background plots should be identical (using the same 
# distribution), then the following plot is the base, with ablines to  
# be added for each section:

# base wt.1 distribution and intervals:

fig <- ggplot(data.frame(mc=raw.dist.mc.sorted),aes(mc))+
  geom_histogram(binwidth=0.002, alpha=0.2)+
  theme_bw() 
fig <- fig + geom_vline(xintercept=int.lower.mc, color="black")
fig <- fig + geom_vline(xintercept=int.upper.mc, color="black")
fig <- fig + geom_vline(xintercept=min(raw.dist.mc.sorted), color="black")
fig <- fig + geom_vline(xintercept=max(raw.dist.mc.sorted), color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))
fig <- fig + 
  xlim(-0.4, 0.3)

fig


# add cis mutant strain comparisons:

alter.genotypes <- c("m76","tata1","m66","tata2")

A <- aggregate(zscore~genotype, data=subset(DATA.FIGURE.3,treatment=="ems"),
               FUN=mc)
for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(length(raw.dist.mc.sorted[which(raw.dist.mc.sorted<result)])/length(raw.dist.mc.sorted))
}

fig <- ggplot(data.frame(mc=raw.dist.mc),aes(mc))
fig <- fig + geom_histogram(binwidth=0.002, alpha=0.2)
fig <- fig + geom_vline(xintercept=A[1,2], color="red")
fig <- fig + geom_vline(xintercept=A[2,2], color="blue")
fig <- fig + geom_vline(xintercept=A[3,2], color="purple")
fig <- fig + geom_vline(xintercept=A[4,2], color="green")
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))

fig
A


#########################################
#####*****Step 4: Trans mutants*****#####
#########################################

###***5.1. Build the ggplot data structure***###
DATA.FIGURE.5 <- subset(DATA.ALL.ALTER, genotype %in% c("wt.3","ade6","rap1","tye7","nam7"))
DATA.FIGURE.5 <- droplevels(DATA.FIGURE.5)
plot.colors <- brewer.pal(5, 'Set3')

### Compare SHAM distributions
fig <- ggplot(subset(DATA.FIGURE.5,treatment=="sham"), aes(x=median.mean,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig

A <- subset(DATA.FIGURE.5,treatment=="sham")
alter.genotypes <- c("wt.3", "ade6","rap1","tye7","nam7")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"median.mean"]
  print(mean(result))
}

str(DATA.FIGURE.5)
sham <-subset(DATA.FIGURE.5,treatment=="sham")

boxplot(median.mean ~ genotype,
        data = sham
)

pw <- sham %>% 
  pairwise_t_test(
    median.mean ~ genotype, pool.sd = FALSE,
    comparisons = list(c("wt.3", "ade6"), c("wt.3", "rap1"),  
                       c("wt.3", "tye7"),  c("wt.3", "nam7")),
    p.adjust.method = "bonferroni"
  )
pw


boxplot(median.sd ~ genotype,
        data = sham
)

pw <- sham %>% 
  pairwise_t_test(
    median.sd ~ genotype, pool.sd = FALSE,
    comparisons = list(c("wt.3", "ade6"), c("wt.3", "rap1"),  
                       c("wt.3", "tye7"),  c("wt.3", "nam7")),
    p.adjust.method = "bonferroni"
  )
pw

# t-tests of sham to ems for each genotype

str(DATA.FIGURE.5)

for (i in 1:length(alter.genotypes)){
  
  y <- subset(DATA.FIGURE.5,genotype==alter.genotypes[i])
  yt <- t_test(y, median.mean~treatment, p.adjust.method = "bonferroni") 
  print(alter.genotypes[i])
  print(yt)
}


for (i in 1:length(alter.genotypes)){
  
  y <- subset(DATA.FIGURE.5,genotype==alter.genotypes[i])
  yt <- t_test(y, median.sd~treatment, p.adjust.method = "bonferroni") 
  print(alter.genotypes[i])
  print(yt)
}


### Compare EMS distributions
fig <- ggplot(subset(DATA.FIGURE.5,treatment=="ems"), aes(x=median.mean,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig

A <- subset(DATA.FIGURE.5,treatment=="ems")
alter.genotypes <- c("wt.3", "ade6","rap1","tye7","nam7")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"median.mean"]
  print(mean(result))
}

pw <- ems %>% 
  pairwise_t_test(
    median.mean ~ genotype, pool.sd = FALSE,
    comparisons = list(c("wt.3", "ade6"), c("wt.3", "rap1"),  
                       c("wt.3", "tye7"),  c("wt.3", "nam7")),
    p.adjust.method = "bonferroni"
  )
pw

pw <- ems %>% 
  pairwise_t_test(
    median.sd ~ genotype, pool.sd = FALSE,
    comparisons = list(c("wt.3", "ade6"), c("wt.3", "rap1"),  
                       c("wt.3", "tye7"),  c("wt.3", "nam7")),
    p.adjust.method = "bonferroni"
  )
pw

### Compare zscore distributions
fig <- ggplot(subset(DATA.FIGURE.5,treatment=="ems"), aes(x=zscore,
  y=genotype, fill=genotype))
fig <- fig + geom_joy2(size=1)
fig <- fig + scale_fill_manual(values=plot.colors)
fig <- fig + theme_classic()
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig

A <- subset(DATA.FIGURE.5,treatment=="ems")
alter.genotypes <- c("wt.3", "ade6","rap1","tye7","nam7")

for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(mean(result))
}

ems <-subset(DATA.FIGURE.5,treatment=="ems")
str(ems)

boxplot(zscore ~ genotype,
        data = ems
)

pw <- ems %>% 
  pairwise_t_test(
    zscore ~ genotype, pool.sd = FALSE,
    comparisons = list(c("wt.3", "ade6"), c("wt.3", "rap1"),  
                       c("wt.3", "tye7"),  c("wt.3", "nam7")),
    p.adjust.method = "bonferroni"
  )
pw

### Compare trans mutant mv to wt.1 (ref.0)

raw.dist.mv <- permutation.test(subset(DATA.ALL.ALTER,genotype=="wt.1"&treatment=="ems")[,"zscore"])

raw.dist.mv.sorted <-sort(raw.dist.mv) #must be sorted to calculate intervals

#write.csv(raw.dist.mv.sorted, file = "x.csv") # (save in whatever format you like)

# If you plan to use one reference mv distribution for all genotypes (as opposed to pulling 
# new distributions each time), run the above once and then do not overwrite these 
# objects for subsequent statistics and plotting.

raw.dist.mv.sorted <-read.csv("raw.dist.mv.sorted.csv", header = FALSE) # to pull from the included file
raw.dist.mv.sorted <-raw.dist.mv.sorted[,1]

int.lower<-raw.dist.mv.sorted[250] #for 2 tailed
int.upper<-raw.dist.mv.sorted[9750] #for 2 tailed

int.lower<-raw.dist.mv.sorted[500] #for 1 tailed
int.upper<-raw.dist.mv.sorted[9500] #for 1 tailed

# For figures, if all background plots should be identical (using the same 
# distribution), then the following plot is the base, with ablines to  
# be added for each section:

# base wt.1 distribution and intervals:

fig <- ggplot(data.frame(mv=raw.dist.mv.sorted),aes(mv))+
  geom_histogram(binwidth=0.1, alpha = 0.2)+
  theme_bw() 
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + geom_vline(xintercept=min(raw.dist.mv.sorted), color="black")
fig <- fig + geom_vline(xintercept=max(raw.dist.mv.sorted), color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))
fig <- fig + 
  xlim(0, 18)

fig

# add trans mutant strain comparisons:

alter.genotypes <- c("ade6","rap1","tye7","nam7")

A <- aggregate(zscore~genotype, data=subset(DATA.FIGURE.5,treatment=="ems"),
  FUN=var)
for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(length(raw.dist.mv.sorted[which(raw.dist.mv.sorted<result)])/length(raw.dist.mv.sorted))
}

fig <- ggplot(data.frame(mv=raw.dist.mv.sorted),aes(mv))
fig <- fig + geom_histogram(binwidth=0.1, alpha=0.2)
fig <- fig + geom_vline(xintercept=9.3201, color="red")
fig <- fig + geom_vline(xintercept=4.3900, color="blue")
fig <- fig + geom_vline(xintercept=1.4783, color="purple")
fig <- fig + geom_vline(xintercept=3.6000, color="green")
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + theme_classic()
fig <- fig + xlim(0, 18)
fig <- fig + theme(
	axis.text.x=element_text(size=10),
	strip.background=element_rect(colour="white"))

fig

A

### Compare trans mutant mc to wt.1 (ref.0)

raw.dist.mc <- permutation.test.mc(subset(DATA.ALL.ALTER,genotype=="wt.1"&treatment=="ems")[,"zscore"])

raw.dist.mc.sorted <-sort(raw.dist.mc) #for intervals

#write.csv(raw.dist.mc.sorted, file = "x.csv") # (save in whatever format you like)

# If you plan to use one reference mc distribution for all genotypes (as opposed to pulling 
# new distributions each time), run the above once and then do not overwrite these 
# objects for subsequent statistics and plotting.

raw.dist.mc.sorted <-read.csv("raw.dist.mc.sorted.csv", header = FALSE) # to pull from the included file
raw.dist.mc.sorted <-raw.dist.mc.sorted[,1]

int.lower<-raw.dist.mc.sorted[250] #for 2 tailed
int.upper<-raw.dist.mc.sorted[9750] #for 2 tailed

int.lower<-raw.dist.mc.sorted[500] #for 1 tailed
int.upper<-raw.dist.mc.sorted[9500] #for 1 tailed

# For figures, if all background plots should be identical (using the same 
# distribution), then the following plot is the base, with ablines to  
# be added for each section:

# base wt.1 distribution and intervals:

fig <- ggplot(data.frame(mc=raw.dist.mc.sorted),aes(mc))+
  geom_histogram(binwidth=0.002, alpha=0.2)+
  theme_bw() 
fig <- fig + geom_vline(xintercept=int.lower.mc, color="black")
fig <- fig + geom_vline(xintercept=int.upper.mc, color="black")
fig <- fig + geom_vline(xintercept=min(raw.dist.mc.sorted), color="black")
fig <- fig + geom_vline(xintercept=max(raw.dist.mc.sorted), color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))
fig <- fig + 
  xlim(-0.4, 0.3)

fig


# add trans mutant strain comparisons:

alter.genotypes <- c("ade6","rap1","tye7","nam7")

A <- aggregate(zscore~genotype, data=subset(DATA.FIGURE.5,treatment=="ems"),
               FUN=mc)
for(i in alter.genotypes) {
  result <- subset(A, genotype==i)[,"zscore"]
  print(length(raw.dist.mc.sorted[which(raw.dist.mc.sorted<result)])/length(raw.dist.mc.sorted))
}

fig <- ggplot(data.frame(mc=raw.dist.mc.sorted),aes(mc))
fig <- fig + geom_histogram(binwidth=0.002, alpha=0.2)
fig <- fig + geom_vline(xintercept=A[1,2], color="red")
fig <- fig + geom_vline(xintercept=A[2,2], color="blue")
fig <- fig + geom_vline(xintercept=A[3,2], color="purple")
fig <- fig + geom_vline(xintercept=A[4,2], color="green")
fig <- fig + geom_vline(xintercept=int.lower, color="black")
fig <- fig + geom_vline(xintercept=int.upper, color="black")
fig <- fig + geom_vline(xintercept=min(raw.dist.mc.sorted), color="black")
fig <- fig + geom_vline(xintercept=max(raw.dist.mc.sorted), color="black")
fig <- fig + theme_classic()
fig <- fig + theme(
  axis.text.x=element_text(size=10),
  strip.background=element_rect(colour="white"))

fig

A

