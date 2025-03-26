# Scripts for analyzing scoring trans mutants
setwd("~/Projects/COMPENSATION/EMS7/Analysis")

########################
# 1. Loading Libraries #
########################
rm(list=ls())
options(warn=-1)

library(flowCore)
library(flowClust)
library(flowViz)
library(pcaPP)
library(mixtools)
library(plyr)
library(robustlmm)
library(plotrix)
library(lme4)
library(MASS)

box <- graphics::box

########################
# 1.1 Useful functions #
########################
CORRECT <- function(x, model, id) {
	# model: can use predict() on
	# id: column one wants to correct
	new.data <- data.frame(
		ROW = as.factor(x["ROW"]),
		RUN = as.factor(x["RUN"]),
		COLUMN = as.factor(x["COLUMN"]));
	ref.data <- data.frame(
		ROW = as.factor("A"),
		RUN = as.factor(1),
		COLUMN = as.factor(1));
	adjust <- predict(model, new.data) - predict(model, ref.data)
	return(as.numeric(x[id]) - adjust)
}

REMOVE.OUTLIER <- function(data, factor, id) {
	levels.names <- levels(data[, factor])
	for(i in levels.names) {
		sub.data <- subset(data, data[, factor] == i & data[, "STRAIN"] == "ctrl")
		L <- median(sub.data[, id]) - 4 * mad(sub.data[, id])
		H <- median(sub.data[, id]) + 4 * mad(sub.data[, id])
		deleted <- which(data[, factor] == i & data[, "STRAIN"] == "ctrl" & 
			(data[, id] < L | data[, id] > H))
		if(length(deleted) != 0) {
			data <- data[-deleted, ]
		}
	}
	return(droplevels(data))
}

######################################
# 2. Quality Control and Corrections #
######################################
parent.dir <- "~/Projects/COMPENSATION/EMS7/Analysis"
setwd(parent.dir)

classes <- c(rep("factor", 12), rep("integer", 4), rep("numeric", 23))
DATA <- read.table("Clean.Data2.txt", header=TRUE, sep="\t", colClasses=classes)

REF <- 0.905811693
NEG <- 0.519116913

DATA[,"log.RNA.MEDIAN"] <- log(((DATA[,"YFP.MEDIAN.FINAL"] - NEG) / (REF - NEG)) + 0.05)
DATA[,"log.RNA.SD"] <- DATA[,"YFP.SD.FINAL"] / ((DATA[,"YFP.MEDIAN.FINAL"] - NEG) + (REF - NEG)*0.05)

#########################################
# 2.1. Remove samples with small counts #
#########################################
DATA <- subset(DATA, COUNTS.FINAL > 1000)
DATA <- droplevels(DATA)

# Get all CONTROL samples
CONTROL <- subset(DATA, STRAIN == "ctrl")
RAW.DATA <- DATA
RAW.CONTROL <- CONTROL

##########################################
# 2.2. Plate correction based on CONTROL #
##########################################

########---2.2.1. Correct for FSC MEDIAN--########

CONTROL <- subset(DATA, STRAIN == "ctrl")

### Boxplot to check CONTROL distributions, outliers
#par(mfrow=c(1, 3))
#boxplot(FSC.MEDIAN.FINAL ~ ROW, data=CONTROL)
#boxplot(FSC.MEDIAN.FINAL ~ RUN, data=CONTROL)
#boxplot(FSC.MEDIAN.FINAL ~ COLUMN, data=CONTROL)

FSC.MEDIAN.CORRECT.MODEL <- rlmer(FSC.MEDIAN.FINAL ~ 1 + (1|RUN) + (1|ROW) + (1|COLUMN), data = CONTROL)
corrected <- apply(DATA, 1, CORRECT, FSC.MEDIAN.CORRECT.MODEL, "FSC.MEDIAN.FINAL")
DATA[, "FSC.MEDIAN.CORRECTED"] <- corrected

### Find out outliers by histogram
CONTROL <- subset(DATA, STRAIN == "ctrl")
#hist(CONTROL$FSC.MEDIAN.CORRECTED, breaks=50)
LOW <- with(CONTROL, median(FSC.MEDIAN.CORRECTED) - 4.5*mad(FSC.MEDIAN.CORRECTED))
HIGH <- with(CONTROL, median(FSC.MEDIAN.CORRECTED) + 4.5*mad(FSC.MEDIAN.CORRECTED))
DATA <- subset(DATA, STRAIN != "ctrl" | (STRAIN == "ctrl" & FSC.MEDIAN.CORRECTED > LOW & FSC.MEDIAN.CORRECTED < HIGH))

DATA <- droplevels(DATA)

### Then rerun the linear model
CONTROL <- subset(DATA, STRAIN == "ctrl")
FSC.MEDIAN.CORRECT.MODEL <- rlmer(FSC.MEDIAN.FINAL ~ 1 + (1|RUN) + (1|ROW) + (1|COLUMN), data = CONTROL)
corrected <- apply(DATA, 1, CORRECT, FSC.MEDIAN.CORRECT.MODEL, "FSC.MEDIAN.FINAL")
DATA[, "FSC.MEDIAN.CORRECTED"] <- corrected

### Check histogram again
CONTROL <- subset(DATA, STRAIN == "ctrl")
#hist(CONTROL$FSC.MEDIAN.CORRECTED, breaks=50)

CONTROL <- subset(DATA, STRAIN == "ctrl")

#par(mfrow=c(1, 3))
#boxplot(FSC.MEDIAN.CORRECTED ~ ROW, data=CONTROL)
#boxplot(FSC.MEDIAN.CORRECTED ~ RUN, data=CONTROL)
#boxplot(FSC.MEDIAN.CORRECTED ~ COLUMN, data=CONTROL)

########---2.2.2. Correct for YFP MEDIAN---########

### Boxplot to check CONTROL distributions, outliers
#par(mfrow=c(1, 3))
#boxplot(log.RNA.MEDIAN ~ ROW, data=CONTROL)
#boxplot(log.RNA.MEDIAN ~ RUN, data=CONTROL)
#boxplot(log.RNA.MEDIAN ~ COLUMN, data=CONTROL)

LRNA.MEDIAN.CORRECTED.MODEL <- rlmer(log.RNA.MEDIAN ~ 1 + (1|ROW) + (1|RUN) + (1+COLUMN), data=CONTROL)
corrected <- apply(DATA, 1, CORRECT, LRNA.MEDIAN.CORRECTED.MODEL, "log.RNA.MEDIAN")
DATA[, "log.RNA.MEDIAN.CORRECTED"] <- corrected

#DATA <- REMOVE.OUTLIER(DATA, "RUN", "log.RNA.MEDIAN.CORRECTED")
#DATA <- REMOVE.OUTLIER(DATA, "ROW", "log.RNA.MEDIAN.CORRECTED")

### Find out outliers by histogram
CONTROL <- subset(DATA, STRAIN == "ctrl")
#hist(CONTROL$log.RNA.MEDIAN.CORRECTED, breaks=50)
LOW <- with(CONTROL, median(log.RNA.MEDIAN.CORRECTED) - 4.5*mad(log.RNA.MEDIAN.CORRECTED))
HIGH <- with(CONTROL, median(log.RNA.MEDIAN.CORRECTED) + 4.5*mad(log.RNA.MEDIAN.CORRECTED))
DATA <- subset(DATA, STRAIN!="ctrl" | (STRAIN=="ctrl" & log.RNA.MEDIAN.CORRECTED > LOW & log.RNA.MEDIAN.CORRECTED < HIGH))

DATA <- droplevels(DATA)

### Rerun the linear model
CONTROL <- subset(DATA, STRAIN == "ctrl")
LRNA.MEDIAN.CORRECTED.MODEL <- rlmer(log.RNA.MEDIAN ~ 1 + (1|ROW) + (1|RUN) + (1+COLUMN), data=CONTROL)
corrected <- apply(DATA, 1, CORRECT, LRNA.MEDIAN.CORRECTED.MODEL, "log.RNA.MEDIAN")
DATA[, "log.RNA.MEDIAN.CORRECTED"] <- corrected

### Check histogram again
CONTROL <- subset(DATA, STRAIN == "ctrl")
#hist(CONTROL$log.RNA.MEDIAN.CORRECTED, breaks=50)

#par(mfrow=c(1, 3))
#boxplot(log.RNA.MEDIAN.CORRECTED ~ ROW, data=CONTROL)
#boxplot(log.RNA.MEDIAN.CORRECTED ~ RUN, data=CONTROL)
#boxplot(log.RNA.MEDIAN.CORRECTED ~ COLUMN, data=CONTROL)

DATA[, "YFP.MEDIAN.CORRECTED"] <- (exp(DATA[,"log.RNA.MEDIAN.CORRECTED"]) - 0.05) * (REF - NEG) + NEG

CONTROL <- subset(DATA, STRAIN == "ctrl")

########---2.2.3. Correct for SD on log scale---########

### Boxplot to check CONTROL distributions, outliers
#par(mfrow=c(1, 3))
#boxplot(log.RNA.SD ~ ROW, data=CONTROL)
#boxplot(log.RNA.SD ~ RUN, data=CONTROL)
#boxplot(log.RNA.SD ~ COLUMN, data=CONTROL)

LRNA.SD.CORRECTED.MODEL <- rlmer(log.RNA.SD ~ 1 + (1|ROW) + (1|RUN) + (1+COLUMN), data=CONTROL)
corrected <- apply(DATA, 1, CORRECT, LRNA.SD.CORRECTED.MODEL, "log.RNA.SD")
DATA[, "log.RNA.SD.CORRECTED"] <- corrected

### Find out outliers by histogram
CONTROL <- subset(DATA, STRAIN == "ctrl")
#hist(CONTROL$log.RNA.SD.CORRECTED, breaks=50)
LOW <- with(CONTROL, median(log.RNA.SD.CORRECTED) - 4.5*mad(log.RNA.SD.CORRECTED))
HIGH <- with(CONTROL, median(log.RNA.SD.CORRECTED) + 4.5*mad(log.RNA.SD.CORRECTED))
DATA <- subset(DATA, STRAIN!="ctrl" | (STRAIN=="ctrl" & log.RNA.SD.CORRECTED > LOW & log.RNA.SD.CORRECTED < HIGH))

DATA <- droplevels(DATA)

### Rerun the linear model
CONTROL <- subset(DATA, STRAIN == "ctrl")
LRNA.SD.CORRECTED.MODEL <- rlmer(log.RNA.SD ~ 1 + (1|ROW) + (1|RUN) + (1+COLUMN), data=CONTROL)
corrected <- apply(DATA, 1, CORRECT, LRNA.SD.CORRECTED.MODEL, "log.RNA.SD")
DATA[, "log.RNA.SD.CORRECTED"] <- corrected

### Check histogram again
CONTROL <- subset(DATA, STRAIN == "ctrl")
#hist(CONTROL$log.RNA.SD.CORRECTED, breaks=50)

#par(mfrow=c(1, 3))
#boxplot(log.RNA.SD.CORRECTED ~ ROW, data=CONTROL)
#boxplot(log.RNA.SD.CORRECTED ~ RUN, data=CONTROL)
#boxplot(log.RNA.SD.CORRECTED ~ COLUMN, data=CONTROL)

DATA[, "YFP.SD.CORRECTED"] <- DATA[,"log.RNA.SD.CORRECTED"] * ((DATA[,"YFP.MEDIAN.CORRECTED"] - NEG) + (REF - NEG)*0.05)

CONTROL <- subset(DATA, STRAIN == "ctrl")

########---2.2.4. Write the summary data---########
write.table(DATA, "DATA.PLATE.CORRECTED.txt", row.names=FALSE, sep="\t", quote=FALSE)


#####################################
# 3. Calculate all necessary values #
#####################################
parent.dir <- "~/Projects/COMPENSATION/EMS7/Analysis"
setwd(parent.dir)

classes <- c(rep("factor", 12), rep("integer", 4), rep("numeric", 23))
classes <- c(classes, rep("numeric", 7))
DATA <- read.table("DATA.PLATE.CORRECTED.txt", sep="\t", header=TRUE, colClasses=classes)

########---3.1. Remove NULL signals.---########
DATA.NULL <- subset(DATA, STRAIN == "null")

MEDIAN.MEDIAN.NULL <- median(DATA.NULL[, "YFP.MEDIAN.CORRECTED"])
SD.MEDIAN.NULL <- median(DATA.NULL[, "YFP.SD.CORRECTED"])

DATA[, "YFP.MEDIAN.NULL.REMOVED"] <- DATA[, "YFP.MEDIAN.CORRECTED"] - MEDIAN.MEDIAN.NULL
DATA[, "YFP.SD.NULL.REMOVED"] <- DATA[, "YFP.SD.CORRECTED"] - SD.MEDIAN.NULL

########---3.2. Convert to [relative-2-WT] scale---########

DATA.WT <- subset(DATA, STRAIN == 'sham.wt')

MEDIAN.MEDIAN.WT <- median(DATA.WT[, "YFP.MEDIAN.NULL.REMOVED"])

DATA[, "YFP.MEDIAN.RELATIVE.2WT"] <- DATA[, "YFP.MEDIAN.NULL.REMOVED"] / MEDIAN.MEDIAN.WT
DATA[, "YFP.SD.SCALED.2WT"] <- DATA[, "YFP.SD.NULL.REMOVED"] / MEDIAN.MEDIAN.WT
DATA[, "YFP.CV"] <- DATA[, "YFP.SD.SCALED.2WT"] / DATA[, "YFP.MEDIAN.RELATIVE.2WT"]

DATA.WT <- subset(DATA, STRAIN == 'sham.wt')

SD.SCALED.MEDIAN.WT <- median(DATA.WT[, "YFP.SD.SCALED.2WT"])
CV.MEDIAN.WT <- median(DATA.WT[, "YFP.CV"])

DATA[, "YFP.SD.RELATIVE.2WT"] <- DATA[, "YFP.SD.SCALED.2WT"] / SD.SCALED.MEDIAN.WT
DATA[, "YFP.CV.RELATIVE.2WT"] <- DATA[, "YFP.CV"] / CV.MEDIAN.WT

########---3.3. Convert to [relative-2-own] scale---########
DATA[, "YFP.MEDIAN.RELATIVE.2OWN"] <- DATA[, "YFP.MEDIAN.RELATIVE.2WT"]
DATA[, "YFP.SD.SCALED.2OWN"] <- DATA[, "YFP.SD.SCALED.2WT"]
DATA[, "YFP.SD.RELATIVE.2OWN"] <- DATA[, "YFP.SD.RELATIVE.2WT"]
DATA[, "YFP.CV.RELATIVE.2OWN"] <- DATA[, "YFP.CV.RELATIVE.2WT"]

# Work with m76 strains
M76.INDICES <- with(DATA, which(STRAIN == 'ems.m76' | STRAIN == 'sham.m76'))
DATA.M76 <- DATA[M76.INDICES, ]

DATA.SHAM.M76 <- subset(DATA.M76, STRAIN == 'sham.m76')

MEDIAN.MEDIAN.SHAM.M76 <- median(DATA.SHAM.M76[, "YFP.MEDIAN.NULL.REMOVED"])

DATA.M76[, "YFP.MEDIAN.RELATIVE.2OWN"] <- DATA.M76[, "YFP.MEDIAN.NULL.REMOVED"] / MEDIAN.MEDIAN.SHAM.M76
DATA.M76[, "YFP.SD.SCALED.2OWN"] <- DATA.M76[, "YFP.SD.NULL.REMOVED"] / MEDIAN.MEDIAN.SHAM.M76

DATA.SHAM.M76 <- subset(DATA.M76, STRAIN == 'sham.m76')

SD.SCALED.MEDIAN.SHAM.M76 <- median(DATA.SHAM.M76[, "YFP.SD.SCALED.2OWN"])
CV.MEDIAN.SHAM.M76 <- median(DATA.SHAM.M76[, "YFP.CV"])

DATA.M76[, "YFP.SD.RELATIVE.2OWN"] <- DATA.M76[, "YFP.SD.SCALED.2OWN"] / SD.SCALED.MEDIAN.SHAM.M76
DATA.M76[, "YFP.CV.RELATIVE.2OWN"] <- DATA.M76[, "YFP.CV"] / CV.MEDIAN.SHAM.M76

DATA[M76.INDICES, ] <- DATA.M76

# Work with tata strains
TATA.INDICES <- with(DATA, which(STRAIN == 'ems.tata' | STRAIN == 'sham.tata'))
DATA.TATA <- DATA[TATA.INDICES, ]

DATA.SHAM.TATA <- subset(DATA.TATA, STRAIN == 'sham.tata')

MEDIAN.MEDIAN.SHAM.TATA <- median(DATA.SHAM.TATA[, "YFP.MEDIAN.NULL.REMOVED"])

DATA.TATA[, "YFP.MEDIAN.RELATIVE.2OWN"] <- DATA.TATA[, "YFP.MEDIAN.NULL.REMOVED"] / MEDIAN.MEDIAN.SHAM.TATA
DATA.TATA[, "YFP.SD.SCALED.2OWN"] <- DATA.TATA[, "YFP.SD.NULL.REMOVED"] / MEDIAN.MEDIAN.SHAM.TATA

DATA.SHAM.TATA <- subset(DATA.TATA, STRAIN == 'sham.tata')

SD.SCALED.MEDIAN.SHAM.TATA <- median(DATA.SHAM.TATA[, "YFP.SD.SCALED.2OWN"])
CV.MEDIAN.SHAM.TATA <- median(DATA.SHAM.TATA[, "YFP.CV"])

DATA.TATA[, "YFP.SD.RELATIVE.2OWN"] <- DATA.TATA[, "YFP.SD.SCALED.2OWN"] / SD.SCALED.MEDIAN.SHAM.TATA
DATA.TATA[, "YFP.CV.RELATIVE.2OWN"] <- DATA.TATA[, "YFP.CV"] / CV.MEDIAN.SHAM.TATA

DATA[TATA.INDICES, ] <- DATA.TATA

write.table(DATA,"DATA.RELATIVE.txt",sep="\t",quote=FALSE,row.names=FALSE)



###########################
# 4. Get strain estimates #
###########################

parent.dir <- "~/Projects/COMPENSATION/EMS7/Analysis"
setwd(parent.dir)

classes <- c(rep("factor", 12), rep("integer", 4), rep("numeric", 23))
classes <- c(classes, rep("numeric", 7))
classes <- c(classes, rep("numeric", 11))
DATA <- read.table("DATA.RELATIVE.TXT", sep="\t", header=TRUE, colClasses=classes)

########---4.1. Remove outliers in replicates---########
#BASED ON FSC
for (i in 1:nrow(DATA))
{
	CUR <- subset(DATA, STRAIN == DATA[i,"STRAIN"])
	
	LOW <- median(CUR$FSC.MEDIAN.CORRECTED) - 5*mad(CUR$FSC.MEDIAN.CORRECTED)
	HIGH <- median(CUR$FSC.MEDIAN.CORRECTED) + 5*mad(CUR$FSC.MEDIAN.CORRECTED)
	
	if (DATA[i,"FSC.MEDIAN.CORRECTED"] > LOW & DATA[i,"FSC.MEDIAN.CORRECTED"] < HIGH)
	{
		DATA[i,"FSC.OUTLIER"] <- "NO"
	} else {
		DATA[i,"FSC.OUTLIER"] <- "YES"
	}
}

#BASED ON YFP MEAN relative to WT
for (i in 1:nrow(DATA))
{
	CUR <- subset(DATA, STRAIN == DATA[i,"STRAIN"])
	
	LOW <- median(CUR$YFP.MEDIAN.RELATIVE.2WT) - 5*mad(CUR$YFP.MEDIAN.RELATIVE.2WT)
	HIGH <- median(CUR$YFP.MEDIAN.RELATIVE.2WT) + 5*mad(CUR$YFP.MEDIAN.RELATIVE.2WT)
	
	if (DATA[i,"YFP.MEDIAN.RELATIVE.2WT"] > LOW & DATA[i,"YFP.MEDIAN.RELATIVE.2WT"] < HIGH)
	{
		DATA[i,"YFP.WT.OUTLIER"] <- "NO"
	} else {
		DATA[i,"YFP.WT.OUTLIER"] <- "YES"
	}
}

#BASED ON YFP MEAN relative to OWN
for (i in 1:nrow(DATA))
{
	CUR <- subset(DATA, STRAIN == DATA[i,"STRAIN"])
	
	LOW <- median(CUR$YFP.MEDIAN.RELATIVE.2OWN) - 5*mad(CUR$YFP.MEDIAN.RELATIVE.2OWN)
	HIGH <- median(CUR$YFP.MEDIAN.RELATIVE.2OWN) + 5*mad(CUR$YFP.MEDIAN.RELATIVE.2OWN)
	
	if (DATA[i,"YFP.MEDIAN.RELATIVE.2OWN"] > LOW & DATA[i,"YFP.MEDIAN.RELATIVE.2OWN"] < HIGH)
	{
		DATA[i,"YFP.OWN.OUTLIER"] <- "NO"
	} else {
		DATA[i,"YFP.OWN.OUTLIER"] <- "YES"
	}
}

########---4.2. Estimate strain expression values---########
DATA.NEW <- subset(DATA, YFP.WT.OUTLIER == "NO" & YFP.OWN.OUTLIER == "NO")

DATA.NEW[, "REPLICATE"] <- 1

RESULT <- aggregate(REPLICATE ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=sum)

RESULT$YFP.MEDIAN.R2WT.MEAN <- aggregate(
	YFP.MEDIAN.RELATIVE.2WT ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=mean)$YFP.MEDIAN.RELATIVE.2WT
RESULT$YFP.MEDIAN.R2WT.SD <- aggregate(
	YFP.MEDIAN.RELATIVE.2WT ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=sd)$YFP.MEDIAN.RELATIVE.2WT
RESULT$YFP.SD.R2WT.MEAN <- aggregate(
	YFP.SD.RELATIVE.2WT ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=mean)$YFP.SD.RELATIVE.2WT
RESULT$YFP.SD.R2WT.SD <- aggregate(
	YFP.SD.RELATIVE.2WT ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=sd)$YFP.SD.RELATIVE.2WT
RESULT$YFP.CV.R2WT.MEAN <- aggregate(
	YFP.CV.RELATIVE.2WT ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=mean)$YFP.CV.RELATIVE.2WT
RESULT$YFP.CV.R2WT.SD <- aggregate(
	YFP.CV.RELATIVE.2WT ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=sd)$YFP.CV.RELATIVE.2WT

RESULT$YFP.MEDIAN.R2OWN.MEAN <- aggregate(
	YFP.MEDIAN.RELATIVE.2OWN ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=mean)$YFP.MEDIAN.RELATIVE.2OWN
RESULT$YFP.MEDIAN.R2OWN.SD <- aggregate(
	YFP.MEDIAN.RELATIVE.2OWN ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=sd)$YFP.MEDIAN.RELATIVE.2OWN
RESULT$YFP.SD.R2OWN.MEAN <- aggregate(
	YFP.SD.RELATIVE.2OWN ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=mean)$YFP.SD.RELATIVE.2OWN
RESULT$YFP.SD.R2OWN.SD <- aggregate(
	YFP.SD.RELATIVE.2OWN ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=sd)$YFP.SD.RELATIVE.2OWN
RESULT$YFP.CV.R2OWN.MEAN <- aggregate(
	YFP.CV.RELATIVE.2OWN ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=mean)$YFP.CV.RELATIVE.2OWN
RESULT$YFP.CV.R2OWN.SD <- aggregate(
	YFP.CV.RELATIVE.2OWN ~ STRAIN + PLATE + POSITION, DATA.NEW, FUN=sd)$YFP.CV.RELATIVE.2OWN

RESULT <- subset(RESULT, STRAIN != "sham.wt" | (YFP.MEDIAN.R2WT.MEAN > 0.95 & YFP.SD.R2WT.MEAN > 0.8))
RESULT <- subset(RESULT, STRAIN != "sham.tata" | YFP.MEDIAN.R2WT.MEAN > 0.9)
RESULT <- subset(RESULT, REPLICATE >= 3)

########---4.3. Hierarchical Permutation Test---########
### Instructions:
### 1. First, get all SHAM data ~ PLATE + POSITION.
### 2. Second, iterate through MUTANT data ~ PLATE + POSITION.
### 3. Final output should be a data frame
RESULT$pvalues.median.r2wt <- 1
RESULT$pvalues.median.r2own <- 1
RESULT$pvalues.sd.r2wt <- 1
RESULT$pvalues.sd.r2own <- 1
RESULT$pvalues.cv.r2wt <- 1
RESULT$pvalues.cv.r2own <- 1

SHUFFLE <- function(null, WT.LIST) {
	index <- sample(1:length(WT.LIST), 1)
	sp <- WT.LIST[[index]]
	n1 <- length(null)
	n2 <- length(sp)
	dist.obs <- abs(mean(null) - mean(sp))

	sim.sample <- sample(c(null, sp), replace=FALSE)
	null.new <- sim.sample[1:n1]
	sp.new <- sim.sample[(n1+1):(n1+n2)]
	dist.sim <- abs(mean(null.new) - mean(sp.new))

	return(dist.obs - dist.sim)
}

N.PERM <- 20000
EMS.RESULT <- subset(RESULT, STRAIN %in% c("ems.m76", "ems.tata"))
N.EMS <- nrow(EMS.RESULT)
WT.M76.RESULT <- subset(RESULT, STRAIN == "sham.m76")
N.WT.M76 <- nrow(WT.M76.RESULT)
WT.TATA.RESULT <- subset(RESULT, STRAIN == "sham.tata")
N.WT.TATA <- nrow(WT.TATA.RESULT)

### For YFP.MEDIAN.RELATIVE.2WT
EMS.LIST <- vector("list", N.EMS)
WT.M76.LIST <- vector("list", N.WT.M76)
WT.TATA.LIST <- vector("list", N.WT.TATA)

for(i in 1:nrow(EMS.RESULT)) {
	EMS.LIST[[i]] <- subset(
		DATA, PLATE==EMS.RESULT[i, "PLATE"] & POSITION==EMS.RESULT[i, "POSITION"])$YFP.MEDIAN.RELATIVE.2WT
}
for(i in 1:nrow(WT.M76.RESULT)) {
	WT.M76.LIST[[i]] <- subset(
		DATA, PLATE==WT.M76.RESULT[i,"PLATE"] & POSITION==WT.M76.RESULT[i,"POSITION"])$YFP.MEDIAN.RELATIVE.2WT
}
for(i in 1:nrow(WT.TATA.RESULT)) {
	WT.TATA.LIST[[i]] <- subset(
		DATA, PLATE==WT.TATA.RESULT[i,"PLATE"] & POSITION==WT.TATA.RESULT[i,"POSITION"])$YFP.MEDIAN.RELATIVE.2WT
}

counter <- 1
for(i in 1:nrow(RESULT)) {
	if(RESULT[i, "STRAIN"] %in% c("sham.wt", "ctrl", "null", "sham.m76", "sham.tata")) {
		next
	}
	else if(RESULT[i, "STRAIN"]=="ems.m76") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.M76.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.median.r2wt"] <- p.val
		counter <- counter + 1
	}
	else if(RESULT[i, "STRAIN"]=="ems.tata") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.TATA.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.median.r2wt"] <- p.val
		counter <- counter + 1
	}
}

### For YFP.MEDIAN.RELATIVE.2OWN
EMS.LIST <- vector("list", N.EMS)
WT.M76.LIST <- vector("list", N.WT.M76)
WT.TATA.LIST <- vector("list", N.WT.TATA)

for(i in 1:nrow(EMS.RESULT)) {
	EMS.LIST[[i]] <- subset(
		DATA, PLATE==EMS.RESULT[i, "PLATE"] & POSITION==EMS.RESULT[i, "POSITION"])$YFP.MEDIAN.RELATIVE.2OWN
}
for(i in 1:nrow(WT.M76.RESULT)) {
	WT.M76.LIST[[i]] <- subset(
		DATA, PLATE==WT.M76.RESULT[i,"PLATE"] & POSITION==WT.M76.RESULT[i,"POSITION"])$YFP.MEDIAN.RELATIVE.2OWN
}
for(i in 1:nrow(WT.TATA.RESULT)) {
	WT.TATA.LIST[[i]] <- subset(
		DATA, PLATE==WT.TATA.RESULT[i,"PLATE"] & POSITION==WT.TATA.RESULT[i,"POSITION"])$YFP.MEDIAN.RELATIVE.2OWN
}

counter <- 1
for(i in 1:nrow(RESULT)) {
	if(RESULT[i, "STRAIN"] %in% c("sham.wt", "ctrl", "null", "sham.m76", "sham.tata")) {
		next
	}
	else if(RESULT[i, "STRAIN"]=="ems.m76") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.M76.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.median.r2own"] <- p.val
		counter <- counter + 1
	}
	else if(RESULT[i, "STRAIN"]=="ems.tata") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.TATA.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.median.r2own"] <- p.val
		counter <- counter + 1
	}
}

### For YFP.SD.RELATIVE.2WT
EMS.LIST <- vector("list", N.EMS)
WT.M76.LIST <- vector("list", N.WT.M76)
WT.TATA.LIST <- vector("list", N.WT.TATA)

for(i in 1:nrow(EMS.RESULT)) {
	EMS.LIST[[i]] <- subset(
		DATA, PLATE==EMS.RESULT[i, "PLATE"] & POSITION==EMS.RESULT[i, "POSITION"])$YFP.SD.RELATIVE.2WT
}
for(i in 1:nrow(WT.M76.RESULT)) {
	WT.M76.LIST[[i]] <- subset(
		DATA, PLATE==WT.M76.RESULT[i,"PLATE"] & POSITION==WT.M76.RESULT[i,"POSITION"])$YFP.SD.RELATIVE.2WT
}
for(i in 1:nrow(WT.TATA.RESULT)) {
	WT.TATA.LIST[[i]] <- subset(
		DATA, PLATE==WT.TATA.RESULT[i,"PLATE"] & POSITION==WT.TATA.RESULT[i,"POSITION"])$YFP.SD.RELATIVE.2WT
}

counter <- 1
for(i in 1:nrow(RESULT)) {
	if(RESULT[i, "STRAIN"] %in% c("sham.wt", "ctrl", "null", "sham.m76", "sham.tata")) {
		next
	}
	else if(RESULT[i, "STRAIN"]=="ems.m76") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.M76.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.sd.r2wt"] <- p.val
		counter <- counter + 1
	}
	else if(RESULT[i, "STRAIN"]=="ems.tata") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.TATA.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.sd.r2wt"] <- p.val
		counter <- counter + 1
	}
}

### For YFP.SD.RELATIVE.2OWN
EMS.LIST <- vector("list", N.EMS)
WT.M76.LIST <- vector("list", N.WT.M76)
WT.TATA.LIST <- vector("list", N.WT.TATA)

for(i in 1:nrow(EMS.RESULT)) {
	EMS.LIST[[i]] <- subset(
		DATA, PLATE==EMS.RESULT[i, "PLATE"] & POSITION==EMS.RESULT[i, "POSITION"])$YFP.SD.RELATIVE.2OWN
}
for(i in 1:nrow(WT.M76.RESULT)) {
	WT.M76.LIST[[i]] <- subset(
		DATA, PLATE==WT.M76.RESULT[i,"PLATE"] & POSITION==WT.M76.RESULT[i,"POSITION"])$YFP.SD.RELATIVE.2OWN
}
for(i in 1:nrow(WT.TATA.RESULT)) {
	WT.TATA.LIST[[i]] <- subset(
		DATA, PLATE==WT.TATA.RESULT[i,"PLATE"] & POSITION==WT.TATA.RESULT[i,"POSITION"])$YFP.SD.RELATIVE.2OWN
}

counter <- 1
for(i in 1:nrow(RESULT)) {
	if(RESULT[i, "STRAIN"] %in% c("sham.wt", "ctrl", "null", "sham.m76", "sham.tata")) {
		next
	}
	else if(RESULT[i, "STRAIN"]=="ems.m76") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.M76.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.sd.r2own"] <- p.val
		counter <- counter + 1
	}
	else if(RESULT[i, "STRAIN"]=="ems.tata") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.TATA.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.sd.r2own"] <- p.val
		counter <- counter + 1
	}
}

### For YFP.CV.RELATIVE.2WT
EMS.LIST <- vector("list", N.EMS)
WT.M76.LIST <- vector("list", N.WT.M76)
WT.TATA.LIST <- vector("list", N.WT.TATA)

for(i in 1:nrow(EMS.RESULT)) {
	EMS.LIST[[i]] <- subset(
		DATA, PLATE==EMS.RESULT[i, "PLATE"] & POSITION==EMS.RESULT[i, "POSITION"])$YFP.CV.RELATIVE.2WT
}
for(i in 1:nrow(WT.M76.RESULT)) {
	WT.M76.LIST[[i]] <- subset(
		DATA, PLATE==WT.M76.RESULT[i,"PLATE"] & POSITION==WT.M76.RESULT[i,"POSITION"])$YFP.CV.RELATIVE.2WT
}
for(i in 1:nrow(WT.TATA.RESULT)) {
	WT.TATA.LIST[[i]] <- subset(
		DATA, PLATE==WT.TATA.RESULT[i,"PLATE"] & POSITION==WT.TATA.RESULT[i,"POSITION"])$YFP.CV.RELATIVE.2WT
}

counter <- 1
for(i in 1:nrow(RESULT)) {
	if(RESULT[i, "STRAIN"] %in% c("sham.wt", "ctrl", "null", "sham.m76", "sham.tata")) {
		next
	}
	else if(RESULT[i, "STRAIN"]=="ems.m76") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.M76.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.cv.r2wt"] <- p.val
		counter <- counter + 1
	}
	else if(RESULT[i, "STRAIN"]=="ems.tata") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.TATA.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.cv.r2wt"] <- p.val
		counter <- counter + 1
	}
}

### For YFP.CV.RELATIVE.2OWN
EMS.LIST <- vector("list", N.EMS)
WT.M76.LIST <- vector("list", N.WT.M76)
WT.TATA.LIST <- vector("list", N.WT.TATA)

for(i in 1:nrow(EMS.RESULT)) {
	EMS.LIST[[i]] <- subset(
		DATA, PLATE==EMS.RESULT[i, "PLATE"] & POSITION==EMS.RESULT[i, "POSITION"])$YFP.CV.RELATIVE.2OWN
}
for(i in 1:nrow(WT.M76.RESULT)) {
	WT.M76.LIST[[i]] <- subset(
		DATA, PLATE==WT.M76.RESULT[i,"PLATE"] & POSITION==WT.M76.RESULT[i,"POSITION"])$YFP.CV.RELATIVE.2OWN
}
for(i in 1:nrow(WT.TATA.RESULT)) {
	WT.TATA.LIST[[i]] <- subset(
		DATA, PLATE==WT.TATA.RESULT[i,"PLATE"] & POSITION==WT.TATA.RESULT[i,"POSITION"])$YFP.CV.RELATIVE.2OWN
}

counter <- 1
for(i in 1:nrow(RESULT)) {
	if(RESULT[i, "STRAIN"] %in% c("sham.wt", "ctrl", "null", "sham.m76", "sham.tata")) {
		next
	}
	else if(RESULT[i, "STRAIN"]=="ems.m76") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.M76.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.cv.r2own"] <- p.val
		counter <- counter + 1
	}
	else if(RESULT[i, "STRAIN"]=="ems.tata") {
		result <- sapply(1:N.PERM, FUN=function(x) SHUFFLE(EMS.LIST[[counter]], WT.TATA.LIST))
		p.val <- (length(which(result<0)) + 1) / (length(which(result!=0)) + 1)
		RESULT[i, "pvalues.cv.r2own"] <- p.val
		counter <- counter + 1
	}
}

write.table(RESULT,"STRAIN.ESTIMATES.txt",sep="\t",quote=FALSE,row.names=FALSE)



