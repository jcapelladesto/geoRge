#' Detect putative Isotope pairs
#'
#' @param XCMSet the xcmsSet with labelled and unlabelled samples.
#' @param XCMSmode="maxo" Which type of peak intensity to use
#' @param ULtag part of the samplename for unlabelled sample
#' @param Ltag part of the samplename for unlabelled sample
#' @param separator="_" separator in samplename
#' @param sep.pos unsure.
#' @param fc.threshold=1.2 Fold-change threshhold between XXX
#' @param p.value.threshold=0.05 p-Value threshhold between XXX
#' @param PuInc.int.lim=5000 Some intensity threshhold
#' @return Potential Isotope pairs
#' @export
PuInc_seeker <-
function (XCMSet=NULL, XCMSmode="maxo", ULtag=NULL, Ltag=NULL, separator="_", sep.pos=NULL, fc.threshold=1.2, 
p.value.threshold=0.05, PuInc.int.lim=5000, ...){

X1 <- groupval(XCMSet, value = XCMSmode) # geoRge has only been tested on "maxo".
D1 <- data.frame(t(X1))
colnames(D1) <- as.character(1:nrow(X1))

classv <- as.factor(XCMSet@phenoData$class) # sample classes (use separate folders per each group when running XCMS)

##' This is a way to extract the condition classes automatically from the phenoData vector
##' 
##' At the moment it requires that the sample classes have a structure like this:
##' Unlabelled/Labelled tag (For example: "Glc12" / "Glc13") +
##' a separator (For example: "_") +
##' sample condition (For example: "WT" / "KO")
##' 
##' This part could be improved forcing the user to use factor vectors. 
##' Any suggestions are welcome. 

conditions <- levels(classv)
if(sep.pos=="f"){
	conditions <- gsub(paste(ULtag,separator,sep=""),"",conditions)[-grep(paste(Ltag,separator,sep=""),conditions)]
	}else{
	conditions <- gsub(paste(separator,ULtag,sep=""),"",conditions)[-grep(paste(separator,Ltag,sep=""),conditions)]
}

filtsamps <- 1:ncol(D1)

meanintensities <- apply(D1, 2, function(x) tapply(x, classv, mean))
meanintensities <- meanintensities[grep(ULtag, rownames(meanintensities)), ]
filtsampsint <- filtsamps[apply(meanintensities, 2, function(x) all(x < PuInc.int.lim))]

#### First step: Compare Labelled vs Unlabelled ####

# Welch's test
pvalues <- sapply(conditions, function(y) {
  apply(D1[ ,filtsampsint], 2, function (x) {
    a <- try(t.test(x[intersect(grep(ULtag,classv),grep(y,classv))],
                    x[intersect(grep(Ltag,classv),grep(y,classv))], var.equal=F)$p.value, silent=T)
    if (is(a, "try-error")|is.na(a)) {a <- 1}
    return(a)
  })
})
pvalues <- data.frame(pvalues)
colnames(pvalues) <- conditions

# Fold change (This fold-change calculation is not commonly used*)
fc.test <- sapply(conditions, function(y) {
  apply(D1[ ,filtsampsint], 2, function (x) { 
    ulm <- mean(x[intersect(grep(ULtag,classv),grep(y,classv))])
    labm <- mean(x[intersect(grep(Ltag,classv),grep(y,classv))]) 
    FC <- labm/ulm
    FC2 <- (-(ulm/labm))
    FC[FC<1] <- FC2[FC<1]
    return(FC)
  })
})
fc.test <- data.frame(fc.test)
colnames(fc.test) <- conditions

##' Ratio between Labelled and Unlabelled samples for a feature to be an incorporation candidate.
##' This value depends on how strict you want to be
##' The higher the threshold, the less incorporations you may find

# Returns the features that are putative incorporations and the conditions where
# the incorporation was detected

compinc <- sapply(1:nrow(pvalues),function(x) {
  t <- which(pvalues[x, ] < p.value.threshold)
  t2 <- which(fc.test[x, ] > fc.threshold)
  if (length(t) == 0 | length(t2) == 0) {
  	return()
  } else {
  	t <- t[is.element(t, t2)]
  if (ncol(pvalues) == 1) { #If only condition is given
  	if (length(t) > 0) {
  	return(conditions)
  	} else {
  	return()
  	}
  } else {
  	return(conditions[t])
  }
}
})

compinc <- lapply(1:length(compinc), function(x) paste(compinc[[x]], collapse=";"))
compinc <- unlist(compinc)
names(compinc) <- colnames(D1[ ,filtsampsint])
compinc  <- compinc[which(compinc != "")]

res_inc <- XCMSet@groups[as.numeric(names(compinc)), c("mzmed", "rtmed", "rtmin", "rtmax")]
rownames(res_inc) <- 1:nrow(res_inc)

return(
list("PuInc"=res_inc, "PuInc_conditions"=compinc, "pvalue"=pvalues, "foldchange"=fc.test)
)
}
