\name{database_query}
\alias{database_query}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
database_query(geoRgeR = NULL, adducts = NULL, db = NULL, ppm.db = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geoRgeR}{
%%     ~~Describe \code{geoRgeR} here~~
}
  \item{adducts}{
%%     ~~Describe \code{adducts} here~~
}
  \item{db}{
%%     ~~Describe \code{db} here~~
}
  \item{ppm.db}{
%%     ~~Describe \code{ppm.db} here~~
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
function (geoRgeR = NULL, adducts = NULL, db = NULL, ppm.db = 10) 
{
    georgedf <- geoRgeR
    DatabaseSearch <- function(input, db, ppm.db, adducts) {
        input <- input
        search_vect <- (input + adducts$AdductMass)
        h <- lapply(search_vect, function(y) {
            mass_range <- c(y - ppm.db * (y/1e+06), y + ppm.db * 
                (y/1e+06))
            suppressWarnings(a <- which(as.numeric(db$Monoisotopic.Mass) >= 
                mass_range[1] & as.numeric(db$Monoisotopic.Mass) <= 
                mass_range[2]))
        })
        h2 <- sapply(1:length(h), function(z) {
            y <- h[[z]]
            if (length(y) > 0) {
                r <- paste(db$Name[y], sep = "", collapse = ";")
            }
            else {
                r <- "No hit"
            }
            return(r)
        })
        return(h2)
    }
    database.hits <- t(apply(georgedf[, c("mzmed", "atoms")], 
        1, function(y) {
            if (as.numeric(y["atoms"]) != 0) {
                return(rep("", times = nrow(adducts)))
            }
            else {
                input <- as.numeric(as.character(y["mzmed"]))
                return(DatabaseSearch(input, db, ppm.db, adducts))
            }
        }))
    colnames(database.hits) <- adducts$Adduct
    return(database.hits)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
