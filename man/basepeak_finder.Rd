\name{basepeak_finder}
\alias{basepeak_finder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
basepeak_finder(PuIncR = NULL, XCMSmode = "maxo", XCMSet = NULL, ULtag = NULL, Ltag = NULL, separator = "_", sep.pos = NULL, UL.atomM = NULL, L.atomM = NULL, ppm.s = NULL, rt.win.min = 1, Basepeak.minInt = NULL, Basepeak.Percentage = 0.7, noise.quant = 0.2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PuIncR}{
%%     ~~Describe \code{PuIncR} here~~
}
  \item{XCMSmode}{
%%     ~~Describe \code{XCMSmode} here~~
}
  \item{XCMSet}{
%%     ~~Describe \code{XCMSet} here~~
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
  \item{UL.atomM}{
%%     ~~Describe \code{UL.atomM} here~~
}
  \item{L.atomM}{
%%     ~~Describe \code{L.atomM} here~~
}
  \item{ppm.s}{
%%     ~~Describe \code{ppm.s} here~~
}
  \item{rt.win.min}{
%%     ~~Describe \code{rt.win.min} here~~
}
  \item{Basepeak.minInt}{
%%     ~~Describe \code{Basepeak.minInt} here~~
}
  \item{Basepeak.Percentage}{
%%     ~~Describe \code{Basepeak.Percentage} here~~
}
  \item{noise.quant}{
%%     ~~Describe \code{noise.quant} here~~
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
function (PuIncR = NULL, XCMSmode = "maxo", XCMSet = NULL, ULtag = NULL, 
    Ltag = NULL, separator = "_", sep.pos = NULL, UL.atomM = NULL, 
    L.atomM = NULL, ppm.s = NULL, rt.win.min = 1, Basepeak.minInt = NULL, 
    Basepeak.Percentage = 0.7, noise.quant = 0.2, ...) 
{
    X1 <- groupval(XCMSet, value = XCMSmode)
    D1 <- data.frame(t(X1))
    colnames(D1) <- as.character(1:nrow(X1))
    filtsamps <- 1:ncol(D1)
    classv <- as.factor(XCMSet@phenoData$class)
    res_inc <- PuIncR[["PuInc"]]
    compinc <- PuIncR[["PuInc_conditions"]]
    max_atoms <- function(mass, uatommass) {
        floor((mass/uatommass))
    }
    xgroup <- cbind(XCMSet@groups[filtsamps, c("mzmed", "rtmed")], 
        t(D1))
    mass_diff <- L.atomM - UL.atomM
    Xn <- groupval(XCMSet, value = "sn")
    Dn <- data.frame(t(Xn))
    Dn <- Dn[rownames(D1), ]
    colnames(Dn) <- as.character(1:ncol(Dn))
    meannoise <- apply(Dn, 2, function(x) tapply(x, classv, mean))
    meannoise[is.na(meannoise)] <- 3
    mn <- quantile(meannoise, noise.quant, na.rm = T)
    filtsampsnoise <- which(apply(meannoise, 2, function(x) any(x > 
        mn)) == T)
    george <- lapply(rownames(res_inc), function(y) {
        isot <- sapply(1:max_atoms(res_inc[y, "mzmed"], L.atomM), 
            function(x) {
                res_inc[y, "mzmed"] - (x * mass_diff)
            })
        isot <- sort(isot)
        rt_d <- (res_inc[y, "rtmax"] - res_inc[y, "rtmin"])
        if (rt_d < rt.win.min) {
            rt_d <- rt.win.min
        }
        rt_range <- c(res_inc[y, "rtmed"] - rt_d, res_inc[y, 
            "rtmed"] + rt_d)
        isot.match <- lapply(isot, function(x) {
            mass_range <- c(x - ppm.s * (x/1e+06), x + ppm.s * 
                (x/1e+06))
            a <- names(which(xgroup[filtsampsnoise, "mzmed"] >= 
                mass_range[1] & xgroup[filtsampsnoise, "mzmed"] <= 
                mass_range[2]))
            b <- names(which(xgroup[filtsampsnoise, "rtmed"] >= 
                rt_range[1] & xgroup[filtsampsnoise, "rtmed"] <= 
                rt_range[2]))
            r <- intersect(a, b)
            return(r)
        })
        isot.match <- unlist(isot.match)
        inc <- names(compinc)[as.numeric(y)]
        isotv <- c(isot.match, inc)
        if (length(isotv) < 2) {
            return()
        }
        else {
            iso_int <- xgroup[isotv, ]
            inc_condition <- compinc[as.numeric(y)]
            mi <- apply(iso_int[, 3:ncol(iso_int)], 1, function(x) tapply(x, 
                classv, mean))
            mi <- mi[, -which(colnames(mi) == inc)]
            cond <- sort(unique(unlist(strsplit(as.character(unique(inc_condition)), 
                split = ";"))))
            if (is.vector(mi)) {
                pos <- sapply(cond, USE.NAMES = F, simplify = T, 
                  function(x) {
                    mi12 <- mi[grep(ULtag, names(mi))]
                    mi12 <- mi12[grep(x, names(mi12))]
                    if ((mi12 > Basepeak.minInt)) {
                      pos <- isot.match
                      return(pos)
                    }
                    else {
                      return()
                    }
                  })
                pos <- unique(unlist(pos))
            }
            else {
                pos <- sapply(cond, USE.NAMES = F, simplify = T, 
                  function(x) {
                    mi12 <- mi[grep(ULtag, rownames(mi)), ]
                    mi12 <- mi12[grep(x, rownames(mi12)), ]
                    if (any(mi12 > Basepeak.minInt)) {
                      mi12 <- mi12[which(mi12 > Basepeak.minInt)]
                      pos <- names(which(mi12 > (Basepeak.Percentage * 
                        max(mi12))))
                      return(pos)
                    }
                    else {
                      return()
                    }
                  })
                pos <- unlist(pos)
                if (is.matrix(pos)) {
                  pos <- as.vector(pos)
                }
                pos <- unique(unlist(pos))
            }
            if (length(pos) < 1) {
                pos <- NULL
            }
            if (is.null(pos)) {
                return()
            }
            else {
                pos <- as.character(pos)
                if (length(pos) > 1) {
                  r <- lapply(length(pos):1, function(x) {
                    pos1 <- setdiff(pos, pos[x])
                    iso_int <- iso_int[which(iso_int[, "mzmed"] >= 
                      min(iso_int[pos[x], "mzmed"])), ]
                    if (any(rownames(iso_int) \%in\% pos1)) {
                      iso_int <- iso_int[-which(rownames(iso_int) \%in\% 
                        pos1), ]
                    }
                    else {
                      iso_int <- iso_int
                    }
                    inc_id <- rep(paste(y, x, sep = "."), times = nrow(iso_int))
                    feature_id <- rownames(iso_int)
                    atoms <- round(((iso_int[, "mzmed"]) - min(iso_int[, 
                      "mzmed"]))/mass_diff, digits = 1)
                    atoms <- sort(unname(atoms))
                    inc_condition <- compinc[feature_id]
                    inc_condition[is.na(inc_condition)] <- ""
                    res <- cbind(as.data.frame(inc_id), as.data.frame(feature_id), 
                      iso_int[, 1:2], as.data.frame(atoms), as.data.frame(inc_condition), 
                      iso_int[, -(1:2)])
                    rownames(res) <- NULL
                    return(res)
                  })
                  res <- do.call("rbind", r)
                  rownames(res) <- NULL
                  return(res)
                }
                else {
                  iso_int <- iso_int[which(iso_int[, "mzmed"] >= 
                    min(iso_int[pos, "mzmed"])), ]
                  inc_id <- rep(y, times = nrow(iso_int))
                  feature_id <- rownames(iso_int)
                  atoms <- round(((iso_int[, "mzmed"]) - min(iso_int[, 
                    "mzmed"]))/mass_diff, digits = 1)
                  atoms <- sort(unname(atoms))
                  inc_condition <- compinc[feature_id]
                  inc_condition[is.na(inc_condition)] <- ""
                  res <- cbind(as.data.frame(inc_id), as.data.frame(feature_id), 
                    iso_int[, 1:2], as.data.frame(atoms), as.data.frame(inc_condition), 
                    iso_int[, -(1:2)])
                  rownames(res) <- NULL
                  return(res)
                }
            }
        }
    })
    georgedf <- do.call("rbind", george)
    return(georgedf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
