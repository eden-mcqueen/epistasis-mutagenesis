#Clear memory
rm(list=ls())
options(warn=-1)


#####################
#1-LOADING LIBRARIES#
#####################

library(flowCore)
library(flowClust)
library(flowViz)
library(plotrix)
library(nlme)
library(MethComp)
library(outliers)
library(pcaPP)

library(reshape2)
library(MASS)
library(ggplot2)
library(Hmisc)
library(fBasics)
library(lawstat)
library(fitdistrplus)
library(mixtools)
library(vioplot)
library(gplots)
library(RColorBrewer)
library(calibrate)

box <- graphics::box


#####################
#2-LOADING FUNCTIONS#
#####################

#Set source directory
source.dir <- "~/Projects/COMPENSATION/EMS7/Scripts/"
setwd(source.dir)

source("Cleaning.Functions.2.R")



##############
#3-CLEAN DATA#
##############

#Set working directory
parent.dir <- "~/Projects/COMPENSATION/EMS7/"
setwd(parent.dir)

#Load experiment setup
classes <- c("character", "character", "character", "character", "character", "character",
	"integer", "character", "integer", "integer", "integer", "integer")
SETUP <- read.table("Analysis/TEMPLATE.csv",sep=",",header=TRUE,as.is=TRUE,colClasses=classes)

#Load paths of FCS files
FILENAMES <- list.files("FCS-WELLS",pattern=".fcs",recursive=TRUE,include.dirs=TRUE)
SETUP[,"FILENAMES"] <- paste("FCS-WELLS", FILENAMES, sep="/")
names(SETUP) <- toupper(names(SETUP))

#Determine Hard Gates
GATES <- GATE.CALIB(SETUP[1, "FILENAMES"])

#Clean Data
#SETUP <- head(SETUP, 200)
Output <- apply(SETUP,1,CLEANING,GATES=GATES)

OUTPUT <- as.data.frame(Output[[1]])

OUTPUT[2:nrow(SETUP),1:ncol(OUTPUT)] <- NA

for (i in 2:nrow(SETUP))
{
	if (is.null(Output[[i]]))
	{} else {
		OUTPUT[i,] <- Output[[i]]
	}
}

write.table(OUTPUT,"Analysis/Experiment.Output2.txt",row.names=FALSE,sep="\t")

CLEAN <- cbind.data.frame(SETUP[,1:(ncol(SETUP)-1)],OUTPUT)

write.table(CLEAN,"Analysis/Clean.Data2.txt",row.names=FALSE,sep="\t")


#########################
#4-QUALITY CONTROL PLOTS#
#########################


#Set working directory
parent.dir <- "~/Projects/COMPENSATION/EMS7/"
setwd(parent.dir)

SETUP <- read.table("Analysis/TEMPLATE.csv",header=TRUE,as.is=TRUE,sep=",",colClasses=classes)
FILENAMES <- list.files("FCS-WELLS",pattern=".fcs",recursive=TRUE,include.dirs=TRUE)
SETUP[,"INITIAL.DATA"] <- paste("FCS-WELLS", FILENAMES, sep="/")
names(SETUP) <- toupper(names(SETUP))
CLEAN.DATA <- paste("CLEAN.DATA/","Day",SETUP[,"DAY"],"_Rep",SETUP[,"REP"],"_Plate",SETUP[,"PLATE"],"_Well",SETUP[,"POSITION"],".txt",sep="")

SETUP <- cbind.data.frame(SETUP,CLEAN.DATA)
SETUP$INITIAL.DATA <- as.character(SETUP$INITIAL.DATA)
SETUP$CLEAN.DATA <- as.character(SETUP$CLEAN.DATA)

#CHECK CONTROL SAMPLES
PLOT <- SETUP[which(SETUP$STRAIN == "ctrl"),]
apply(PLOT,1,FCS.PLOT)
