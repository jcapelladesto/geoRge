#' Find MS features in metabolite database
#' 
#' Longer Description of database_query
#' @param geoRgeR result of basepeak_finder
#' @param adducts list of adducts in the current ion mode
#' @param db Database in form of a CSV
#' @param ppm.db=10 default search PPM in Database
#' @export

database_query <-
function(geoRgeR=NULL, adducts=NULL, db=NULL, ppm.db=10) {

georgedf <- geoRgeR
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
