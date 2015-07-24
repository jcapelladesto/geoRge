PuInc <- function (XCMSet=NULL,ULtag=NULL,Ltag=NULL,separator=NULL,sep.pos=NULL, fc.threshold=1.2, p.value.threshold=0.05,
PuInc.int.lim=5000,...){


XCMSet <- xset3

X1 <- groupval(XCMSet, value = "maxo") # geoRge has only been used on "maxo". I suppose it works for others too.
D1 <- data.frame(t(X1))
colnames(D1) <- as.character(1:nrow(X1))

classv <- as.factor(XCMSet@phenoData$class) # sample classes
classv <- as.factor(gsub("CELL_","",classv))
# If you have Quality Control samples. We do not need them.
D1 <- D1[grep("QC",rownames(D1),invert=T),] 
classv <- as.factor(as.character(classv)[grep("QC",as.character(classv),invert=T)])

##' This is a way to extract the condition classes automatically from the phenoData vector
##' 
##' At the moment it requires that the sample classes have a structure like this:
##' Unlabelled/Labelled tag (For example: "Glc12" / "Glc13") +
##' a separator (For example: "_") +
##' sample condition (For example: "WT" / "KO")
##' 
##' This part could be improved forcing the user to use factor vectors. 
##' Any suggestions are welcome. 

ULtag <- "Glc12"
Ltag <- "Glc13"
separator <- "_"
sep.pos <- "f" # front "b" back

conditions <- sapply(levels(classv),USE.NAMES=F,function(x){
	if(sep.pos=="f"){
	a <- strsplit(x,split=paste(ULtag,separator,sep=""))[[1]][2]
	}else{
	a <- strsplit(x,split=paste(separator,ULtag,sep=""))[[1]][1]
	}
	return (a)
})
conditions <- unique(conditions[-which(is.na(conditions))])
#conditions <- conditions[c(2,4)]
## Filter features that are not found in at least NUM samples of one of the condition groups ##

# filtsamps <- unlist(sapply(1:ncol(D1),function(x){
# 	if(any(XCMSet@groups[x,grep(Ltag,colnames(XCMSet@groups))]>=3)){
# 		return(x)
# 	}
# }))

filtsamps <- 1:ncol(D1)

#D1 <- D1[,filtsamps]
D1[is.na(D1)] <- 200 #500

meanintensities <- apply(D1, 2, function(x) tapply(x, classv,mean))
meanintensities <- meanintensities[grep(ULtag,rownames(meanintensities)),]
# meanintensities <- meanintensities[c(2,4),]
filtsampsint <- filtsamps[apply(meanintensities, 2, function(x) all(x<5000))]

#plot(x=seq(500,20000,500),y=sapply(seq(500,20000,500),function(z) {length(filtsamps[apply(meanintensities, 2, function(x) all(x<z))])}))

Xn <- groupval(XCMSet, value = "sn") 
Dn <- data.frame(t(Xn))
Dn <- Dn[rownames(D1),]
colnames(Dn) <- as.character(1:ncol(Dn))
meannoise <- apply(Dn, 2, function(x) tapply(x, classv,mean))
meannoise[is.na(meannoise)] <- 3
mn <- quantile(meannoise,0.2, na.rm=T) #quintela 0.33
filtsampsnoise  <- which(apply(meannoise, 2, function(x) any(x>mn))==T) #8270

# meannoise <- sapply(xset3@groupidx,function(x) mean(xset3@peaks[x,"sn"],na.rm=T))
# filtsampsnoise <- which(meannoise>quantile(meannoise,0.2,na.rm=T)) # 8014

### 7662 common
#setdiff(filtsampsnoise.cutre,filtsampsnoise.strict)
#setdiff(filtsampsnoise.strict,filtsampsnoise.cutre)

# feature 8154 (cutre ok) = strict is more accurate :) + 608 specific that would be lost otherwise

#### First step: Compare Labelled vs Unlabelled ####

# Welch's test
pvalues <- sapply(conditions,function(y){
  apply(D1[,filtsampsint],2,function (x){
    a <- try(t.test(x[which(classv==paste(ULtag,y,sep=separator))],
                    x[which(classv==paste(Ltag,y,sep=separator))],var.equal=F)$p.value,silent=T)
    if(is(a,"try-error")|is.na(a)){a <- 1}
    return(a)
  })
})
pvalues <- data.frame(pvalues)
colnames(pvalues) <- conditions

# Fold change (This fold-change calculation is not commonly used*)
fc.test <- sapply(conditions,function(y){
  apply(D1[,filtsampsint],2,function (x){ 
    ulm <- mean(x[which(classv==paste(ULtag,y,sep=separator))])
    labm <- mean(x[which(classv==paste(Ltag,y,sep=separator))]) 
    FC <- labm/ulm
    FC2 <- -ulm/labm
    FC[FC<1] <- FC2[FC<1]
    return(FC)
  })
})
fc.test <- data.frame(fc.test)
colnames(fc.test) <- conditions

##' Ratio between Labelled and Unlabelled samples for a feature to be an incorporation candidate.
##' This value depends on how strict you want to be
##' The higher the threshold, the less incorporations you may find

fc.threshold <- 1.2

##' P-value threshold for the Welch's test. 0.05 as default.
p.value.threshold <- 0.05

# Returns the features that are putative incorporations and the conditions where
# the incorporation was detected

compinc <- sapply(1:nrow(pvalues),function(x){
  t <- which(pvalues[x,]<p.value.threshold)
  t2 <- which(fc.test[x,]>fc.threshold)
  if(length(t)==0|length(t2)==0){
  	return()
  }else{
  	t <- t[is.element(t, t2)]
  if(ncol(pvalues)==1){ #If only condition is given
  	if(length(t)>0){
  	return(conditions)
  	}else{
  	return()
  	}
  }else{
  	return(conditions[t])
  }
}
})

compinc <- lapply(1:length(compinc),function(x)paste(compinc[[x]],collapse=";"))
compinc <- unlist(compinc)
names(compinc) <- colnames(D1[,filtsampsint])
compinc  <- compinc[which(compinc!="")]

res_inc <- XCMSet@groups[as.numeric(names(compinc)),c("mzmed","rtmed","rtmin","rtmax")]
rownames(res_inc) <- 1:nrow(res_inc)
return(list())
}
