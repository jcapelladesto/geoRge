#' @import xcms
#' @importClassesFrom Rcpp "C++Object"
NULL # This is required so that roxygen knows where the first manpage starts

#' Detect Putative Incorporation features
#'
#' Use Welch's T-test to detect features upregulated in labelled samples
#' @param XCMSet xcmsSet object of samples processed by xcms workflow
#' @param XCMSmode Type of quantification value for the features, as in \code{xcms::groupval(XCMSet)}
#' @param ULtag Unlabelled sample name tag, for example "C12"
#' @param Ltag Labelled sample name tag, for example "C13"
#' @param separator Character used as separator in sample names, as in \code{xcms::sampclass(XCMSet)}
#' @param sep.pos.front \code{TRUE} if position of \code{ULtag} and \code{Ltag} is found \emph{previous} to the condition names. Check \code{xcms::sampclass(XCMSet)}
#' @param fc.threshold Fold-change threshold for detection of upregulated features associated to Putative incorporation 
#' @param p.value.threshold p-Value threshhold between for detection of upregulated features associated to Putative incorporation
#' @param PuInc.int.lim Intensity limit in unlabelled sample for the selection of features to be tested as possible Putative incorporation 
#' @return A list containing: \describe{
#'   \item{PuInc}{Putative Incorporation Features matrix}
#'   \item{PuInc_conditions}{Experimetal conditions (as in \code{xcms::sampclass(XCMSet)}) in which PuInc were detected}
#'   \item{pvalue}{Matrix of p-values from the tested features}
#'   \item{foldchange}{Matrix of fold-changes from the tested features}  
#'   \item{params}{Inputed arguments to be inherited by \code{\link{basepeak_finder}}}
#' }
#' @examples
#' data(mtbls213)
#' xcms::sampclass(mtbls213)
#' PuInc_seeker(XCMSet=mtbls213,ULtag="CELL_Glc12",Ltag="CELL_Glc13",
#' sep.pos.front=TRUE ,fc.threshold=1.5,p.value.threshold=.05,	PuInc.int.lim = 4000)
#' @export
PuInc_seeker <-
function (XCMSet=NULL, XCMSmode="maxo", ULtag=NULL, Ltag=NULL, separator="_", sep.pos.front=TRUE, fc.threshold=1.2, 
p.value.threshold=0.05, PuInc.int.lim=NULL){

X1 <- xcms::groupval(XCMSet, value = XCMSmode) # geoRge has only been tested on "maxo".
D1 <- data.frame(t(X1))
colnames(D1) <- as.character(1:nrow(X1))

classv <- as.factor(xcms::sampclass(XCMSet)) # sample classes (use separate folders per each group when running XCMS or set them before running geoRge)
cond_class_split <- sapply(as.character(classv),function(x){
	pos <- 1
	if(sep.pos.front){pos <- pos+1}
	strsplit(x,split="-")[[1]][pos]})
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
if(sep.pos.front){
	conditions <- gsub(paste(ULtag,separator,sep=""),"",conditions)[-grep(paste(Ltag,separator,sep=""),conditions)]
	}else{
	conditions <- gsub(paste(separator,ULtag,sep=""),"",conditions)[-grep(paste(separator,Ltag,sep=""),conditions)]
}

filtsamps <- 1:ncol(D1)

meanintensities <- apply(D1, 2, function(x) tapply(x, classv,mean))
meanintensities <- as.data.frame(meanintensities[grep(ULtag,rownames(meanintensities)),]) ## changed 1 condition

if (length(conditions)==1){ ## changed 1 condition
	filtsampsint <- filtsamps[apply(meanintensities, 1, function(x) all(x< PuInc.int.lim))]
}else{
	filtsampsint <- filtsamps[apply(meanintensities, 2, function(x) all(x< PuInc.int.lim))]
}

#### First step: Compare Labelled vs Unlabelled ####

# Welch's test
pvalues <- sapply(conditions, function(y) {
  apply(D1[ ,filtsampsint], 2, function (x) {
    a <- try(t.test(x[intersect(grep(Ltag,classv),which(cond_class_split==x))],
                    x[intersect(grep(ULtag,classv),which(cond_class_split==x))], var.equal=F)$p.value, silent=T)
    if (is(a, "try-error")|is.na(a)) {a <- 1}
    return(a)
  })
})
pvalues <- data.frame(pvalues)
colnames(pvalues) <- conditions

# Fold change (This fold-change calculation is not commonly used*)
fc.test <- sapply(conditions, function(y) {
  apply(D1[ ,filtsampsint], 2, function (x) { 
    ulm <- mean(x[intersect(grep(ULtag,classv),which(cond_class_split==x))])
    labm <- mean(x[intersect(grep(Ltag,classv),which(cond_class_split==x))]) 
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
list("PuInc"=res_inc, "PuInc_conditions"=compinc, "pvalue"=pvalues, "foldchange"=fc.test,
		 "params" = list("XCMSmode"=XCMSmode, "ULtag"=ULtag, "Ltag"=Ltag, "separator"=separator,
		 "sep.pos.front"=sep.pos.front,"conditions"=conditions))
)
}
