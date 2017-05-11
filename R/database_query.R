#' Find Monoisotopic Peaks in a Custom Metabolite Database
#' 
#' Look for hits for mz of monoisotopic peaks in a database in .CSV format
#' @param geoRgeR Result of basepeak_finder
#' @param adducts List of adducts in the current ion mode
#' @param db Database in .CSV format (see \code{system.file("extdata/ExampleDatabase.csv", package="geoRge")})
#' @param ppm.db ppm value to search in Database
#' @export

database_query <-
function(geoRgeR=NULL, adducts=NULL, db=NULL, ppm.db=10) {

georgedf <- geoRgeR$geoRge
#### search base peaks in a database ####

##' This function requires having a database that can be read in R as a data frame. Also
##' geoRge provides txt files (that are read as data frame) with the adducts masses that 
##' is used for 'adducts' argument. The 'db' data frame must contain at least 2 columns:
##' Name
##' Monoisotopic.Mass

## Function to perform the search
DatabaseSearch <- function(input, db, ppm.db, adducts) {
	input <- input
	search_vect <- (input + adducts$AdductMass)
	h <- lapply(search_vect, function(y) {
		mass_range <- c(y - ppm.db * (y/1e6),y + ppm.db * (y/1e6))
		suppressWarnings(a <- which(as.numeric(db$Monoisotopic.Mass) >= mass_range[1] & 
					as.numeric(db$Monoisotopic.Mass) <= mass_range[2]))
	})

	h2 <- sapply(1:length(h), function(z) {
		y <- h[[z]]	
		if (length(y) > 0) {
			r <- paste(db$Name[y], sep="", collapse=";")
		} else {
			r <- "No hit"
		}
		return(r)
	})
	
	return(h2)
}

## Search results
database.hits <- t(apply(georgedf[,c("mzmed","atoms")], 1, function(y) {
	if(as.numeric(y["atoms"]) != 0) {
		return(rep("", times=nrow(adducts)))
	} else {
		input <- as.numeric(as.character(y["mzmed"]))
		return(DatabaseSearch(input, db, ppm.db, adducts))
	}
})
)
colnames(database.hits) <- adducts$Adduct
return(database.hits)
}
