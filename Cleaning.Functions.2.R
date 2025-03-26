
  ##########################
### 1-ROTATION OF FCS DATA ###
  ##########################
   
ROT <- function(x,Rotation){
	Result <- Rotation%*%x
	return(Result)
}

#--------------------------------------------------------------------------------------------------------------------------------#

  ###########################
### 2-HARD GATE CALIBRATION ###
  ###########################

GATE.CALIB <- function(x,logFSC.A.MIN,logFSC.A.MAX,logFSC.H.MIN,logFSC.H.MAX,FSC.A_FSC.H.MIN,FSC.A_FSC.H.MAX,Width.MIN,Width.MAX) {

	if (missing(logFSC.A.MIN)) {logFSC.A.MIN <- 5.1}
	if (missing(logFSC.A.MAX)) {logFSC.A.MAX <- 6.5}
	if (missing(logFSC.H.MIN)) {logFSC.H.MIN <- 5.2}
	if (missing(logFSC.H.MAX)) {logFSC.H.MAX <- 6.8}
	if (missing(FSC.A_FSC.H.MIN)) {FSC.A_FSC.H.MIN <- 0.88}
	if (missing(FSC.A_FSC.H.MAX)) {FSC.A_FSC.H.MAX <- 0.94}
	if (missing(Width.MIN)) {Width.MIN <- 30}
	if (missing(Width.MAX)) {Width.MAX <- 80}
	
	Merge.Frame <- read.FCS(x,transformation=FALSE,alter.names=TRUE)
		
	##############################	
	##Log transformation of data##
	##############################
	
	logTrans <- logTransform(transformationId="log10-transformation",logbase=10,r=1,d=1)
	Merge.Frame <- transform(Merge.Frame,`logFSC.A`=logTrans(`FSC.A`))
	Merge.Frame <- transform(Merge.Frame,`logFSC.H`=logTrans(`FSC.H`))
	Merge.Frame <- transform(Merge.Frame,`logFL1.A`=logTrans(`FL1.A`))
	Merge.Frame <- transform(Merge.Frame,`logFL1.H`=logTrans(`FL1.H`))
	
	####################################
	##Calculate phenotypes of interest##
	####################################
	
	Data.Fluo <- as.data.frame(exprs(Merge.Frame))
	Data.Fluo[Data.Fluo == 0] <- 1
	
	Phenotype3 <- Data.Fluo[,"logFL1.A"]/Data.Fluo[,"logFSC.A"]
	Phenotype3 <- as.matrix(Phenotype3)
	colnames(Phenotype3) <- "FL1/FSC"
	Merge.Frame <- cbind2(Merge.Frame, Phenotype3)
	
	Phenotype4 <- (Data.Fluo[,"logFSC.A"])/(Data.Fluo[,"logFSC.H"])
	Phenotype4 <- as.matrix(Phenotype4)
	colnames(Phenotype4) <- "FSC.A/FSC.H"
	Merge.Frame <- cbind2(Merge.Frame, Phenotype4)		
	
	PlotAll <- as.data.frame(exprs(Merge.Frame))
	
	
	#############################
	###Quick plot to set gates###
	#############################
	
	# #Quick plot to set gates
	# quartz(height=14,width=14)
	par(mfrow=c(2,2))
	
	plot(PlotAll[,"Width"],PlotAll[,"logFSC.A"],pch=20,cex=0.5,col="#00000022",xlim=c(20,100),ylim=c(4,7))
	abline(v=Width.MIN)
	abline(v=Width.MAX)
	abline(h=logFSC.A.MIN)
	abline(h=logFSC.A.MAX)
	
	plot(PlotAll[,"logFSC.A"],PlotAll[,"logFL1.A"],pch=20,cex=0.5,col="#00000022",xlim=c(4,7),ylim=c(1.5,6.5))
	abline(v=logFSC.A.MIN)
	abline(v=logFSC.A.MAX)
	
	plot(PlotAll[,"logFSC.H"],PlotAll[,"logFSC.A"],pch=20,cex=0.4,col="#00000022",xlim=c(3,8),ylim=c(4,7))
	abline(v=logFSC.H.MIN)
	abline(v=logFSC.H.MAX)
	abline(h=logFSC.A.MIN)
	abline(h=logFSC.A.MAX)
	
	plot(PlotAll[,"logFSC.A"],PlotAll[,"FSC.A/FSC.H"],pch=20,cex=0.4,col="#00000022",xlim=c(3,8))
	abline(v=logFSC.A.MIN)
	abline(v=logFSC.A.MAX)
	abline(h=FSC.A_FSC.H.MIN)
	abline(h=FSC.A_FSC.H.MAX)
	
	
	 
	OUTPUT <- c(logFSC.A.MIN,logFSC.A.MAX,logFSC.H.MIN,logFSC.H.MAX,FSC.A_FSC.H.MIN,FSC.A_FSC.H.MAX,Width.MIN,Width.MAX)
	names(OUTPUT) <- c("logFSC.A.MIN","logFSC.A.MAX","logFSC.H.MIN","logFSC.H.MAX","FSC.A_FSC.H.MIN","FSC.A_FSC.H.MAX","Width.MIN","Width.MAX")

	return(OUTPUT)
}


#--------------------------------------------------------------------------------------------------------------------------------#

  #######################
### 3-CLEANING FCS DATA ###
  #######################


CLEANING <- function(x,GATES) {
	
	Merge.Frame <- read.FCS(x["FILENAMES"],transformation=FALSE,alter.names=TRUE)
	
	OUTPUT <- data.frame(matrix(nrow=1,ncol=length(c("COUNTS.INITIAL",	"COUNTS.GATES",	"COUNTS.SINGLES",	"COUNTS.FINAL", "FSC.KURTOSIS",	"WIDTH",	"FSC.MEDIAN.INITIAL",	"FSC.MAD.INITIAL",	"FL1.MEDIAN.INITIAL",	"FL1.MAD.INITIAL",	"YFP.MEDIAN.INITIAL",	"YFP.MAD.INITIAL", "YFP.SD.INITIAL",	"INTERCEPT.INITIAL",	"SLOPE.INITIAL", "THETA", "YFP.MEDIAN.ROT","YFP.MAD.ROT","YFP.SD.ROT",	"FSC.MEDIAN.FINAL",	"FSC.MAD.FINAL",	"YFP.MEDIAN.FINAL",	"YFP.MAD.FINAL", "YFP.SD.FINAL","log.YFP.MEDIAN","log.YFP.MAD","log.YFP.SD")
)))
	colnames(OUTPUT) <- c("COUNTS.INITIAL",	"COUNTS.GATES",	"COUNTS.SINGLES",	"COUNTS.FINAL", "FSC.KURTOSIS",	"WIDTH",	"FSC.MEDIAN.INITIAL",	"FSC.MAD.INITIAL",	"FL1.MEDIAN.INITIAL",	"FL1.MAD.INITIAL",	"YFP.MEDIAN.INITIAL",	"YFP.MAD.INITIAL", "YFP.SD.INITIAL",	"INTERCEPT.INITIAL",	"SLOPE.INITIAL", "THETA", "YFP.MEDIAN.ROT","YFP.MAD.ROT","YFP.SD.ROT",	"FSC.MEDIAN.FINAL",	"FSC.MAD.FINAL",	"YFP.MEDIAN.FINAL",	"YFP.MAD.FINAL", "YFP.SD.FINAL","log.YFP.MEDIAN","log.YFP.MAD","log.YFP.SD")
		
	OUTPUT["COUNTS.INITIAL"] <- nrow(exprs(Merge.Frame))
	
	if (OUTPUT["COUNTS.INITIAL"] > 1000)	{
		
		##############################	
		##Log transformation of data##
		##############################
		
		Start.exp <- exprs(Merge.Frame)
		Start.exp[,"FL1.A"] <- Start.exp[,"FL1.A"] + 10
		
		Merge.Frame <- new("flowFrame",Start.exp)
		
		Merge.Frame <- transform(Merge.Frame,`logFSC.A`=logTrans(`FSC.A`))
		Merge.Frame <- transform(Merge.Frame,`logFSC.H`=logTrans(`FSC.H`))
		Merge.Frame <- transform(Merge.Frame,`logFL1.A`=logTrans(`FL1.A`))
		Merge.Frame <- transform(Merge.Frame,`logFL1.H`=logTrans(`FL1.H`))
		
		####################################
		##Calculate phenotypes of interest##
		####################################
		
		Data.Fluo <- as.data.frame(exprs(Merge.Frame))
		Data.Fluo[Data.Fluo == 0] <- 1
		
		Phenotype3 <- Data.Fluo[,"logFL1.A"]/Data.Fluo[,"logFSC.A"]
		Phenotype3 <- as.matrix(Phenotype3)
		colnames(Phenotype3) <- "YFP.INITIAL"
		Merge.Frame <- cbind2(Merge.Frame, Phenotype3)
		
		Phenotype4 <- (Data.Fluo[,"logFSC.A"])/(Data.Fluo[,"logFSC.H"])
		Phenotype4 <- as.matrix(Phenotype4)
		colnames(Phenotype4) <- "FSC.A/FSC.H"
		Merge.Frame <- cbind2(Merge.Frame, Phenotype4)		
		
		PlotAll <- as.data.frame(exprs(Merge.Frame))	
	
		############
		#Hard Gates# 
		############
		
		rectGate <- rectangleGate(filterId="Noise Removal","logFSC.A"=c(GATES["logFSC.A.MIN"],GATES["logFSC.A.MAX"]), "logFSC.H"=c(GATES["logFSC.H.MIN"],GATES["logFSC.H.MAX"]), "FSC.A/FSC.H"=c(GATES["FSC.A_FSC.H.MIN"],GATES["FSC.A_FSC.H.MAX"]),"Width"=c(GATES["Width.MIN"],GATES["Width.MAX"]),"YFP.INITIAL"=c(min(PlotAll[,"YFP.INITIAL"]),max(PlotAll[,"YFP.INITIAL"])),"FL1.A"=c(11,max(PlotAll[,"FL1.A"])))
		
		Hard.Gates <- Subset(Merge.Frame, rectGate)
		Hard.Gates.exp <- exprs(Hard.Gates)

		OUTPUT["COUNTS.GATES"] <- nrow(Hard.Gates.exp)
	
		####################
		#Doublet Hard Gates# 
		####################
		
		Doublet.Model <- PCAgrid(cbind(Hard.Gates.exp[,"logFSC.A"],Hard.Gates.exp[,"FSC.A/FSC.H"]),k=2,method="sd",scores=TRUE,center="median")
		
		Scores <- Doublet.Model$scores
		
		Distri <- normalmixEM2comp(Scores[,2],sigsqrd=c(0.0022,0.0068)^2,mu=c(-0.0013,0.0086),lambda=c(0.56,0.44))
		
		Lambda <- Distri$lambda
		Mu <- Distri$mu
		Sigma <- Distri$sigma 
		
		Order <- c(which(Mu == min(Mu)),which(Mu == max(Mu)))
		
		Lambda <- Lambda[Order]
		Mu <- Mu[Order]
		Sigma <- Sigma[Order] 
		
		#Good cluster
		f <- function(x) dnorm(x,m=Mu[1],sd=Sigma[1])*Lambda[1]-dnorm(x,m=Mu[2],sd=Sigma[2])*Lambda[2]
		Threshold <- try(uniroot(f,interval=c(Mu[1],Mu[2]+Sigma[2]))$root)
		
		# if (Threshold > 0.1)
		# {
			# Threshold <- 0.075
		# }
			
		#Remove big cells based on FSC.A/FSC.H
		Position <- which(Scores[,2] < Threshold)
		Doublet.Gates.exp <- Hard.Gates.exp[Position,]
		Doublet.Gates <- new("flowFrame",Doublet.Gates.exp)
		
		#Remove cells with extreme FSC.A 
		DOUBLETS <- Doublet.Gates.exp[,"logFSC.A"]
		
		MEDIAN <- median(DOUBLETS)
		MAD <- mad(DOUBLETS)
		LOW <- MEDIAN - 2*MAD
		HIGH <- MEDIAN + 2*MAD
	
		NEW.MEDIAN <- MEDIAN
		OLD.MEDIAN <- median(DOUBLETS[which(DOUBLETS > LOW & DOUBLETS < HIGH)])
		
		while (abs(NEW.MEDIAN-OLD.MEDIAN) > 0.001)
		{
			NEW.DOUBLETS <- DOUBLETS[which(DOUBLETS > LOW & DOUBLETS < HIGH)]
			OLD.MEDIAN <- NEW.MEDIAN
			NEW.MEDIAN <- median(NEW.DOUBLETS)
			NEW.MAD <- mad(NEW.DOUBLETS)
			
			LOW <- NEW.MEDIAN - 2*NEW.MAD
			HIGH <- NEW.MEDIAN + 2*NEW.MAD
		}

		rectGate <- rectangleGate(filterId="Outliers logFSC.A","logFSC.A"=c(LOW,HIGH))
		
		Final.Doublets <- Subset(Doublet.Gates, rectGate)
		
		OUTPUT["FSC.KURTOSIS"] <- kurtosis(DOUBLETS)

		#################
		#Singles Cluster#
		#################
		
		Doublet.filter <- flowClust(Final.Doublets,varNames=c("logFSC.H","logFSC.A"),K=1,B=50,min.count=1000,nu.est=2,trans=0,seed=10,z.cutoff=0,level=0.9,tol=1e-4)
		Well.pop <- split(Final.Doublets,Doublet.filter,population=list(sc1=1))
		Well.Doublet <- Well.pop$sc1
		
		Doublets.exp <- as.data.frame(exprs(Well.Doublet))
						
		OUTPUT["COUNTS.SINGLES"] <- nrow(Doublets.exp)
		
		# plot(Hard.Gates.exp[,"logFSC.A"],Hard.Gates.exp[,"logFSC.H"],pch=20,col="#00000066")
		# points(Doublets.exp[,"logFSC.A"],Doublets.exp[,"logFSC.H"],pch=20,col="#FF000099")
		
		#############
		#Fluo filter#
		#############
	
		#REMOVE OUTLIERS FL1/FSC
		FLUO <- Doublets.exp[,"YFP.INITIAL"]
		
		MEDIAN <- median(FLUO)
		MAD <- mad(FLUO)
		LOW <- MEDIAN - 4*MAD
		HIGH <- MEDIAN + 4*MAD
	
		NEW.MEDIAN <- MEDIAN
		OLD.MEDIAN <- median(FLUO[which(FLUO > LOW & FLUO < HIGH)])
		
		while (abs(NEW.MEDIAN-OLD.MEDIAN) > 0.001)
		{
			NEW.FLUO <- FLUO[which(FLUO > LOW & FLUO < HIGH)]
			OLD.MEDIAN <- NEW.MEDIAN
			NEW.MEDIAN <- median(NEW.FLUO)
			NEW.MAD <- mad(NEW.FLUO)
			
			LOW <- NEW.MEDIAN - 4*NEW.MAD
			HIGH <- NEW.MEDIAN + 4*NEW.MAD
		}
	
		rectGate <- rectangleGate(filterId="Outliers FL1/FSC Removal","YFP.INITIAL"=c(LOW,HIGH))
		
		Hard.Fluo <- Subset(Well.Doublet, rectGate)
	
		Gate.Fluo <- flowClust(Hard.Fluo, varNames=c("logFSC.A","logFL1.A"),K=1,B=50,min.count=1000,nu.est=1,trans=0,z.cutoff=0.5,seed=10,tol=1e-5,nu=1.5,level=0.98)
		Well.pop <- split(Hard.Fluo,Gate.Fluo,population=list(sc1=1))
		Well.Fluo <- Well.pop$sc1
		Fluo.exp <- as.data.frame(exprs(Well.Fluo))

		#SAVE DATA
		OUTPUT["COUNTS.FINAL"] <- nrow(Fluo.exp)
		OUTPUT["WIDTH"] <- median(Fluo.exp[,"Width"])
		OUTPUT["FSC.MEDIAN.INITIAL"] <- median(Fluo.exp[,"logFSC.A"])
		OUTPUT["FSC.MAD.INITIAL"] <- mad(Fluo.exp[,"logFSC.A"])
		OUTPUT["FL1.MEDIAN.INITIAL"] <- median(Fluo.exp[,"logFL1.A"])
		OUTPUT["FL1.MAD.INITIAL"] <- mad(Fluo.exp[,"logFL1.A"])	
		OUTPUT["YFP.MEDIAN.INITIAL"] <- median(Fluo.exp[,"YFP.INITIAL"])
		OUTPUT["YFP.MAD.INITIAL"] <- mad(Fluo.exp[,"YFP.INITIAL"])
		OUTPUT["YFP.SD.INITIAL"] <- sd(Fluo.exp[,"YFP.INITIAL"])
		
		# plot(Doublets.exp$logFSC.A,Doublets.exp$logFL1.A,pch=20,col="#00000066")
		# points(Fluo.exp$logFSC.A,Fluo.exp$logFL1.A,pch=20,col="#FF000099")		
		
		###############################################################
		#####Remove correlation between logFSC.A and FL1.A/FSC.A#######
		###############################################################
		
		#1-Define orthogonal regression
		Intercept <- c()
		Slope <- c()
		Theta <- c()
		
		Fluo.Model <- PCAgrid(cbind(Fluo.exp[,"logFSC.A"],Fluo.exp[,"YFP.INITIAL"]),k=2,method="sd",scores=FALSE,center="median")
		
		#2-Center of rotation
		x.center <- Fluo.Model$center[1]
		y.center <- Fluo.Model$center[2]
				
		#3-Initial Intercept and Slope
		Slope[1] <- Fluo.Model$loadings[2,1] / Fluo.Model$loadings[1,1]
		Intercept[1] <- Fluo.Model$center[2] - Slope[1]*Fluo.Model$center[1]
				
		#4-Calculate angle of rotation
		a <- c(x.center-0,y.center-Intercept[1]) #Vector from Intercept to Centroid
		b <- c(x.center-0,y.center-y.center) #Vector with slope 0 through Centroid
				
		Theta[1] <- acos(sum(a*b)/(sqrt(sum(a*a))*sqrt(sum(b*b)))) #Angle between 2 vectors
				
		if (Slope[1] < 0)
		{
			Theta[1] <- -Theta[1]
		}		
				
		#5-Define rotation matrix
		Rotation <- matrix(c(cos(Theta[1]),-sin(Theta[1]),sin(Theta[1]),cos(Theta[1])),ncol=2,nrow=2)
	
		#6-Transform Data
		Coord <- t(as.matrix(Fluo.exp[,c("logFSC.A","YFP.INITIAL")]))
		
		Coord[1,] <- Coord[1,] - x.center
		Coord[2,] <- Coord[2,] - y.center
		
		Result <- ROT(x=Coord,Rotation=Rotation)
		
		Result[1,] <- Result[1,] + x.center
		Result[2,] <- Result[2,] + y.center
		
		#7-Keep record of rotated values
		Fluo.exp[,"FSC.FINAL"] <- Result[1,]
		Fluo.exp[,"YFP.ROT"] <- Result[2,]
		
		
		###########################################################################
		#Apply correction for linear relation between fluorescence and mRNA levels#
		###########################################################################
				
		REF <- 0.905811693
		NEG <- 0.519116913
				
		for (j in 1:nrow(Fluo.exp))
		{
			Fluo.exp[j,"YFP.FINAL"] <- (exp((Fluo.exp[j,"YFP.ROT"]-0.905274742)*log(10)/0.294448097) - 0.05) * (REF - NEG) + NEG
		}
	
		OUTPUT["INTERCEPT.INITIAL"] <- Intercept[1]
		OUTPUT["SLOPE.INITIAL"] <- Slope[1]
		OUTPUT["THETA"] <- Theta[1]
		OUTPUT["YFP.MEDIAN.ROT"] <- median(Fluo.exp[,"YFP.ROT"])
		OUTPUT["YFP.MAD.ROT"] <- mad(Fluo.exp[,"YFP.ROT"])
		OUTPUT["YFP.SD.ROT"] <- sd(Fluo.exp[,"YFP.ROT"])
		OUTPUT["FSC.MEDIAN.FINAL"] <- median(Fluo.exp[,"FSC.FINAL"])
		OUTPUT["FSC.MAD.FINAL"] <- mad(Fluo.exp[,"FSC.FINAL"])	
		OUTPUT["YFP.MEDIAN.FINAL"] <- median(Fluo.exp[,"YFP.FINAL"])
		OUTPUT["YFP.MAD.FINAL"] <- mad(Fluo.exp[,"YFP.FINAL"])	
		OUTPUT["YFP.SD.FINAL"] <- sd(Fluo.exp[,"YFP.FINAL"])	
		OUTPUT["log.YFP.MEDIAN"] <- median(log(Fluo.exp[,"YFP.FINAL"]))
		OUTPUT["log.YFP.MAD"] <- mad(log(Fluo.exp[,"YFP.FINAL"]))	
		OUTPUT["log.YFP.SD"] <- sd(log(Fluo.exp[,"YFP.FINAL"]))	
		
		#Save clean data with rotation correction
		Data <- Fluo.exp[,c("Width","Time","logFSC.A","logFSC.H","logFL1.A","logFL1.H","YFP.INITIAL","YFP.ROT","FSC.FINAL","YFP.FINAL")]
		write.table(Data,file=paste("CLEAN.DATA/","Day",x["DAY"],"_Rep",x["REP"],"_Plate",x["PLATE"],"_Well",x["POSITION"],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)

		return(OUTPUT)
		
	} else {
		OUTPUT[] <- NA
		OUTPUT["COUNTS.INITIAL"] <- nrow(exprs(Merge.Frame))
		
		Data <- as.data.frame(matrix(data=NA,ncol=10))
		colnames(Data) <- c("Width","Time","logFSC.A","logFSC.H","logFL1.A","logFL1.H","YFP.INITIAL","YFP.ROT","FSC.FINAL","YFP.FINAL")
		write.table(Data,file=paste("CLEAN.DATA/","Day",x["DAY"],"_Rep",x["REP"],"_Plate",x["PLATE"],"_Well",x["POSITION"],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
		
		return(OUTPUT)
	}
	
}


#--------------------------------------------------------------------------------------------------------------------------------#

  ###########################
### 4-QUALITY CONTROL PLOTS ###
  ###########################


FCS.PLOT <- function(x) {
	
	INITIAL.DATA <- read.FCS(x["INITIAL.DATA"],transformation=FALSE,alter.names=TRUE)
	CLEAN.DATA <-  read.table(x["CLEAN.DATA"],header=TRUE,as.is=TRUE)

	if(nrow(CLEAN.DATA) < 100) {
		return(-1)
	}

	##############################	
	##Log transformation of data##
	##############################
		
	Merge.Frame <- transform(INITIAL.DATA,`logFSC.A`=logTrans(`FSC.A`))
	Merge.Frame <- transform(Merge.Frame,`logFSC.H`=logTrans(`FSC.H`))
	Merge.Frame <- transform(Merge.Frame,`logFL1.A`=logTrans(`FL1.A`))
	Merge.Frame <- transform(Merge.Frame,`logFL1.H`=logTrans(`FL1.H`))
		
	####################################
	##Calculate phenotypes of interest##
	####################################
		
	Data.Fluo <- as.data.frame(exprs(Merge.Frame))
	Data.Fluo[Data.Fluo == 0] <- 1
		
	Phenotype3 <- Data.Fluo[,"logFL1.A"]/Data.Fluo[,"logFSC.A"]
	Phenotype3 <- as.matrix(Phenotype3)
	colnames(Phenotype3) <- "YFP.INITIAL"
	Merge.Frame <- cbind2(Merge.Frame, Phenotype3)
		
	Phenotype4 <- (Data.Fluo[,"logFSC.A"])/(Data.Fluo[,"logFSC.H"])
	Phenotype4 <- as.matrix(Phenotype4)
	colnames(Phenotype4) <- "FSC.A/FSC.H"
	Merge.Frame <- cbind2(Merge.Frame, Phenotype4)		
		
	PLOT.ALL <- as.data.frame(exprs(Merge.Frame))

	NAME <- paste("D",x["DAY"],"P",x["PLATE"],"R",x["REP"],"_",x["POSITION"],"_",x["STRAIN"],sep="")
	

	pdf(paste("CLEANING.PLOTS/",NAME,".pdf",sep=""))

	par(mfrow=c(2,2))

	plot(PLOT.ALL$logFSC.A,PLOT.ALL$logFSC.H,pch=20,cex=0.5,col="#00000033",xlim=c(4.3,6.5),ylim=c(5,7),xlab="logFSC.A",ylab="logFSC.H",main=NAME)
	points(CLEAN.DATA$logFSC.A,CLEAN.DATA$logFSC.H,pch=20,cex=0.5,col="#FF000066",main=NAME)

	plot(PLOT.ALL$logFSC.A,PLOT.ALL$logFL1.A,pch=20,cex=0.5,col="#00000033",xlim=c(4.3,6.5),ylim=c(1.5,6.5),xlab="logFSC.A",ylab="logFL1.A",main=NAME)
	points(CLEAN.DATA$logFSC.A,CLEAN.DATA$logFL1.A,pch=20,cex=0.5,col="#FF000066",main=NAME)

	plot(CLEAN.DATA$logFSC.A,CLEAN.DATA$YFP.INITIAL,pch=20,cex=0.5,col="#FF0000AA",xlab="logFSC.A",ylab="FL1/FSC",main=NAME)
	abline(lm(CLEAN.DATA$YFP.INITIAL~CLEAN.DATA$logFSC.A),lty=2,col="red",lwd=2.5)
	points(CLEAN.DATA$FSC.FINAL,CLEAN.DATA$YFP.ROT,pch=20,cex=0.5,col="#00FF00AA",main=NAME)
	abline(lm(CLEAN.DATA$YFP.ROT~CLEAN.DATA$FSC.FINAL),lty=2,col="green",lwd=2.5)
	
	hist(CLEAN.DATA$YFP.FINAL,breaks=50,main=NAME,xlab="FL1/FSC")

	dev.off()

}
	
	
#--------------------------------------------------------------------------------------------------------------------------------#

logTrans <- logTransform(transformationId="log10-transformation",logbase=10,r=1,d=1)

