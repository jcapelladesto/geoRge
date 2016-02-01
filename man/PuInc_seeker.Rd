\name{PuInc_seeker}
\alias{PuInc_seeker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PuInc_seeker(XCMSet = NULL, XCMSmode = "maxo", ULtag = NULL, Ltag = NULL, separator = "_", sep.pos = NULL, fc.threshold = 1.2, p.value.threshold = 0.05, PuInc.int.lim = 5000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XCMSet}{
%%     ~~Describe \code{XCMSet} here~~
}
  \item{XCMSmode}{
%%     ~~Describe \code{XCMSmode} here~~
}
  \item{ULtag}{
%%     ~~Describe \code{ULtag} here~~
}
  \item{Ltag}{
%%     ~~Describe \code{Ltag} here~~
}
  \item{separator}{
%%     ~~Describe \code{separator} here~~
}
  \item{sep.pos}{
%%     ~~Describe \code{sep.pos} here~~
}
  \item{fc.threshold}{
%%     ~~Describe \code{fc.threshold} here~~
}
  \item{p.value.threshold}{
%%     ~~Describe \code{p.value.threshold} here~~
}
  \item{PuInc.int.lim}{
%%     ~~Describe \code{PuInc.int.lim} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (XCMSet = NULL, XCMSmode = "maxo", ULtag = NULL, Ltag = NULL, 
    separator = "_", sep.pos = NULL, fc.threshold = 1.2, p.value.threshold = 0.05, 
    PuInc.int.lim = 5000, ...) 
{
    X1 <- groupval(XCMSet, value = XCMSmode)
    D1 <- data.frame(t(X1))
    colnames(D1) <- as.character(1:nrow(X1))
    classv <- as.factor(XCMSet@phenoData$class)
    conditions <- levels(classv)
    if (sep.pos == "f") {
        conditions <- gsub(paste(ULtag, separator, sep = ""), 
            "", conditions)[-grep(paste(Ltag, separator, sep = ""), 
            conditions)]
    }
    else {
        conditions <- gsub(paste(separator, ULtag, sep = ""), 
            "", conditions)[-grep(paste(separator, Ltag, sep = ""), 
            conditions)]
    }
    filtsamps <- 1:ncol(D1)
    meanintensities <- apply(D1, 2, function(x) tapply(x, classv, 
        mean))
    meanintensities <- meanintensities[grep(ULtag, rownames(meanintensities)), 
        ]
    filtsampsint <- filtsamps[apply(meanintensities, 2, function(x) all(x < 
        PuInc.int.lim))]
    pvalues <- sapply(conditions, function(y) {
        apply(D1[, filtsampsint], 2, function(x) {
            a <- try(t.test(x[intersect(grep(ULtag, classv), 
                grep(y, classv))], x[intersect(grep(Ltag, classv), 
                grep(y, classv))], var.equal = F)$p.value, silent = T)
            if (is(a, "try-error") | is.na(a)) {
                a <- 1
            }
            return(a)
        })
    })
    pvalues <- data.frame(pvalues)
    colnames(pvalues) <- conditions
    fc.test <- sapply(conditions, function(y) {
        apply(D1[, filtsampsint], 2, function(x) {
            ulm <- mean(x[intersect(grep(ULtag, classv), grep(y, 
                classv))])
            labm <- mean(x[intersect(grep(Ltag, classv), grep(y, 
                classv))])
            FC <- labm/ulm
            FC2 <- (-(ulm/labm))
            FC[FC < 1] <- FC2[FC < 1]
            return(FC)
        })
    })
    fc.test <- data.frame(fc.test)
    colnames(fc.test) <- conditions
    compinc <- sapply(1:nrow(pvalues), function(x) {
        t <- which(pvalues[x, ] < p.value.threshold)
        t2 <- which(fc.test[x, ] > fc.threshold)
        if (length(t) == 0 | length(t2) == 0) {
            return()
        }
        else {
            t <- t[is.element(t, t2)]
            if (ncol(pvalues) == 1) {
                if (length(t) > 0) {
                  return(conditions)
                }
                else {
                  return()
                }
            }
            else {
                return(conditions[t])
            }
        }
    })
    compinc <- lapply(1:length(compinc), function(x) paste(compinc[[x]], 
        collapse = ";"))
    compinc <- unlist(compinc)
    names(compinc) <- colnames(D1[, filtsampsint])
    compinc <- compinc[which(compinc != "")]
    res_inc <- XCMSet@groups[as.numeric(names(compinc)), c("mzmed", 
        "rtmed", "rtmin", "rtmax")]
    rownames(res_inc) <- 1:nrow(res_inc)
    return(list(PuInc = res_inc, PuInc_conditions = compinc, 
        pvalue = pvalues, foldchange = fc.test))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
