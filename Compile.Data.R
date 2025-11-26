##########################################################
### This script is used to compile the data into one   ###    
### file for Analysis 1.13                             ###
###                                                    ###
### Date: 2025-02-27                                   ###
### Version: 1.0                                       ###
### Author: BY and EWM                                 ###
##########################################################

###***0.1. Set working directory and clean working space***###
rm(list=ls())
options(warn=-1)


parent.dir <- "C:/Users/edenm/OneDrive/Desktop/Projects/Current_Projects/Epistasis_paper/Epistasis Mutagenesis Paper/Epistasis Mutagenesis Paper/Scripts and Data"
output.dir <- "C:/Users/edenm/OneDrive/Desktop/Projects/Current_Projects/Epistasis_paper/Epistasis Mutagenesis Paper/Epistasis Mutagenesis Paper/Scripts and Data"

library(dplyr)

###***1.0. Load all datasets***###

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
###* Columns in original files have a less intuitive naming system. This step
###* renames them for interpretabilty.

DATA.0 <- subset(DATA.0.LOAD, CLASS!="NULL", select=c("CLASS",
                                                      "YFP.MEDIAN.RELATIVE.MEAN","YFP.MEDIAN.RELATIVE.SD",
                                                      "YFP.SD.RELATIVE.MEAN","YFP.SD.RELATIVE.SD"))

selected.cols <- c("STRAIN",
                   "YFP.MEDIAN.R2OWN.MEAN","YFP.MEDIAN.R2OWN.SD",
                   "YFP.SD.R2OWN.MEAN","YFP.SD.R2OWN.SD")


new.names <- c("strain",
               "median.mean","median.sd",
               "sd.mean","sd.sd")


colnames(DATA.0) <- new.names
DATA.1 <- subset(DATA.1.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.1) <- new.names
DATA.2 <- subset(DATA.2.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.2) <- new.names
DATA.3 <- subset(DATA.3.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.3) <- new.names
DATA.4 <- subset(DATA.4.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.4) <- new.names
DATA.5 <- subset(DATA.5.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.5) <- new.names
DATA.6 <- subset(DATA.6.LOAD, !(STRAIN %in% c("ctrl","null")), select=selected.cols); colnames(DATA.6) <- new.names

###***1.2. Label different genotypes***###
###* The same plate configuration file was reused for each flow run, so every 
###* file above appears to consist of the same set of genotypes, even though
###* each run was different.
###* The code below converts the labels to the correct genotype for each run.

DATA.0$genotype <- "wt.1"

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

###***1.3. Assemble the datasets***###
###* Compiles the data into one master file with the corrected genotype labels

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

treatment <- c()
for(t in DATA.ALL$strain) {
  if(t %in% c("EMS","ems.m76","ems.tata")) { treatment <- c(treatment,"ems") }
  else if(t %in% c("SHAM","sham.wt","sham.m76","sham.tata")) { treatment <- c(treatment,"sham") }
}

DATA.ALL$treatment <- treatment

DATA.ALL$genotype <- as.factor(DATA.ALL$genotype)
DATA.ALL$expids <- as.factor(DATA.ALL$expids)
DATA.ALL$treatment <- as.factor(DATA.ALL$treatment)

DATA.ALL <- droplevels(DATA.ALL)

write.csv(DATA.ALL, file = "SUMMARY.PER.ISOLATE.csv", row.names = FALSE) 
