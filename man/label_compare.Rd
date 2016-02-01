\name{label_compare}
\alias{label_compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
label_compare(geoRgeR = NULL, XCMSmode = "maxo", PuIncR = NULL, XCMSet = NULL, ULtag = NULL, Ltag = NULL, separator = "_", sep.pos = NULL, UL.atomM = NULL, L.atomM = NULL, ppm.s = NULL, rt.win.min = 1, control.cond = NULL, fc.vs.Control = 1, p.value.vs.Control = 0.05, Show.bp = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geoRgeR}{
%%     ~~Describe \code{geoRgeR} here~~
}
  \item{XCMSmode}{
%%     ~~Describe \code{XCMSmode} here~~
}
  \item{PuIncR}{
%%     ~~Describe \code{PuIncR} here~~
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
  \item{control.cond}{
%%     ~~Describe \code{control.cond} here~~
}
  \item{fc.vs.Control}{
%%     ~~Describe \code{fc.vs.Control} here~~
}
  \item{p.value.vs.Control}{
%%     ~~Describe \code{p.value.vs.Control} here~~
}
  \item{Show.bp}{
%%     ~~Describe \code{Show.bp} here~~
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
function (geoRgeR = NULL, XCMSmode = "maxo", PuIncR = NULL, XCMSet = NULL, 
    ULtag = NULL, Ltag = NULL, separator = "_", sep.pos = NULL, 
    UL.atomM = NULL, L.atomM = NULL, ppm.s = NULL, rt.win.min = 1, 
    control.cond = NULL, fc.vs.Control = 1, p.value.vs.Control = 0.05, 
    Show.bp = T, ...) 
{
    georgedf <- geoRgeR
    X1 <- groupval(XCMSet, value = XCMSmode)
    D1 <- data.frame(t(X1))
    colnames(D1) <- as.character(1:nrow(X1))
    filtsamps <- 1:ncol(D1)
    classv <- as.factor(XCMSet@phenoData$class)
    xgroup <- cbind(XCMSet@groups[filtsamps, c("mzmed", "rtmed")], 
        t(D1))
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
    mass_diff <- L.atomM - UL.atomM
    percent.incorp <- lapply(unique(georgedf$inc_id), function(y) {
        inc_id_features <- georgedf[which(georgedf$inc_id == 
            y), ]
        inc_id_int <- inc_id_features[, 7:ncol(inc_id_features)]
        rts <- inc_id_features$rtmed
        rt_range <- c(min(rts), max(rts))
        inc_isot <- max(inc_id_features$atoms) + 1
        isot_m <- inc_id_features[1, "mzmed"] + (inc_isot * mass_diff)
        isot_id <- lapply(isot_m, function(x) {
            mass_range <- c(x - ppm.s * (x/1e+06), x + ppm.s * 
                (x/1e+06))
            a <- which(xgroup[, "mzmed"] >= mass_range[1] & xgroup[, 
                "mzmed"] <= mass_range[2])
            b <- which(xgroup[, "rtmed"] >= rt_range[1] & xgroup[, 
                "rtmed"] <= rt_range[2])
            r <- intersect(a, b)
            return(r)
        })
        isot_id <- unlist(isot_id)
        if (length(isot_id) < 1) {
            all_id <- xgroup[as.character(inc_id_features$feature_id), 
                ]
        }
        else {
            isot_id <- isot_id[1]
            all_id <- c(as.character(inc_id_features$feature_id), 
                isot_id)
            all_id <- xgroup[all_id, ]
        }
        all_id_int <- all_id[, 3:ncol(all_id)]
        inc_percent <- sapply(conditions, function(x) {
            inc_id_intL <- inc_id_int[, intersect(grep(Ltag, 
                colnames(inc_id_int)), grep(x, colnames(inc_id_int)))]
            all_id_intL <- all_id_int[, intersect(grep(Ltag, 
                colnames(all_id_int)), grep(x, colnames(all_id_int)))]
            inc_cal <- sapply(1:ncol(inc_id_intL), function(x) {
                (inc_id_intL[, x]/sum(all_id_intL[, x])) * 100
            })
            return(inc_cal)
        })
        colnames(inc_percent) <- conditions
        atoms <- inc_id_features$atoms
        rownames(inc_percent) <- rep(atoms, length.out = nrow(inc_percent))
        mean_inc <- sapply(conditions, USE.NAMES = T, simplify = T, 
            function(x) {
                sapply(atoms, function(z) {
                  inc_p_v <- inc_percent[which(rownames(inc_percent) == 
                    z), x]
                  inc_p_m <- mean(inc_p_v)
                  return(inc_p_m)
                })
            })
        colnames(mean_inc) <- paste0(colnames(mean_inc), "_MEAN")
        sd_inc <- sapply(conditions, USE.NAMES = T, simplify = T, 
            function(x) {
                sapply(atoms, function(z) {
                  inc_p_v <- inc_percent[which(rownames(inc_percent) == 
                    z), x]
                  inc_p_s <- sd(inc_p_v)
                  return(inc_p_s)
                })
            })
        colnames(sd_inc) <- paste0(colnames(sd_inc), "_SD")
        noncontrol <- setdiff(conditions, control.cond)
        pvals <- sapply(noncontrol, function(x) {
            sapply(atoms, function(y) {
                a <- try(t.test(inc_percent[which(rownames(inc_percent) == 
                  y), which(conditions == x)], inc_percent[which(rownames(inc_percent) == 
                  y), which(conditions == control.cond)], var.equal = T)$p.value, 
                  silent = T)
                if (is(a, "try-error")) {
                  a <- 1
                }
                return(a)
            })
        })
        fct <- sapply(noncontrol, simplify = T, function(x) {
            sapply(1:nrow(mean_inc), function(y) {
                case <- mean_inc[y, which(conditions == x)]
                control <- mean_inc[y, which(conditions == control.cond)]
                FC <- case/control
                FC2 <- (-(control/case))
                FC[FC < 1] <- FC2[FC < 1]
                names(FC) <- NULL
                return(FC)
            })
        })
        comp <- sapply(1:nrow(pvals), function(x) {
            t <- which(pvals[x, ] < p.value.vs.Control)
            if (length(t) == 0) {
                t <- ""
                return(t)
            }
            else {
                up <- which(fct[x, names(t)] > fc.vs.Control)
                down <- which(fct[x, names(t)] < (-fc.vs.Control))
                if (length(up) != 0) {
                  names(t)[up] <- paste("UP", names(t)[up], sep = "_")
                }
                if (length(down) != 0) {
                  names(t)[down] <- paste("DOWN", names(t)[down], 
                    sep = "_")
                }
                return(names(t))
            }
        })
        if (is.matrix(comp)) {
            comp <- sapply(1:ncol(comp), function(x) {
                a <- paste(comp[, x], collapse = ";")
            })
        }
        else {
            comp <- lapply(1:length(comp), function(x) paste(comp[[x]], 
                collapse = ";"))
            comp <- unlist(comp)
        }
        comp[1] <- "Base peak"
        r <- data.frame(Comparison = comp, mean_inc, sd_inc)
        if (!Show.bp) {
            r[1, ] <- rep("Base Peak", times = ncol(r))
        }
        return(r)
    })
    percent.incorpdf <- do.call("rbind", percent.incorp)
    return(percent.incorpdf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
