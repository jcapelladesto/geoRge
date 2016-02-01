basepeak_finder <-
function(PuIncR=NULL, XCMSmode="maxo", XCMSet=NULL, ULtag=NULL, Ltag=NULL, separator="_", sep.pos=NULL,
UL.atomM=NULL, L.atomM=NULL, ppm.s=NULL, rt.win.min=1, Basepeak.minInt=NULL, Basepeak.Percentage=0.7, noise.quant=0.2,  ...) {

X1 <- groupval(XCMSet, value = XCMSmode) # geoRge has only been used on "maxo". I suppose it works for others too.
D1 <- data.frame(t(X1))
colnames(D1) <- as.character(1:nrow(X1))
filtsamps <- 1:ncol(D1)

classv <- as.factor(XCMSet@phenoData$class) # sample classes (use separate folders per group when running XCMS)

res_inc <- PuIncR[["PuInc"]]
compinc <- PuIncR[["PuInc_conditions"]]

#### Second step: Seek the base peaks for each putative incorporation ####

# The ppm window that is used when seeking base peaks. Must be slightly higher than the
# mass accuracy of your instrument.

##' The retention time window at which the isotopologues may elute. We have
##' observed that sometimes the features are fit in a retention time window
##' narrower than expected.


##' When more than one feature is a candidate base peak.
##' A second threshold is used: a Percentage of the max mean intensity. 
##' This filters unexpected intense isotopic pattern peaks.


# Estimate max number of atoms for a given mass
max_atoms <- function(mass, uatommass) {
	# mass = incorporation mass
	# uatommass = mass of the unlabeled atom
	floor((mass/uatommass))
}

xgroup <- cbind(XCMSet@groups[filtsamps, c("mzmed", "rtmed")], t(D1)) 

mass_diff <- L.atomM - UL.atomM

## Signal-to-noise filtering
Xn <- groupval(XCMSet, value = "sn") 
Dn <- data.frame(t(Xn))
Dn <- Dn[rownames(D1), ]
colnames(Dn) <- as.character(1:ncol(Dn))
meannoise <- apply(Dn, 2, function(x) tapply(x, classv, mean))
meannoise[is.na(meannoise)] <- 3 
mn <- quantile(meannoise, noise.quant, na.rm=T) 
filtsampsnoise  <- which(apply(meannoise, 2, function(x) any(x > mn))==T)

george <- lapply(rownames(res_inc), function(y) {
	
	# Calculate aprox mass of the isotopologues
	isot <- sapply(1:max_atoms(res_inc[y, "mzmed"], L.atomM), function(x) {res_inc[y, "mzmed"]-(x * mass_diff)})
	isot <- sort(isot)
	
	# Set retention time window
	rt_d <- (res_inc[y, "rtmax"]-res_inc[y, "rtmin"])
	if (rt_d<rt.win.min) {rt_d <- rt.win.min}
	rt_range <- c(res_inc[y, "rtmed"] - rt_d, res_inc[y, "rtmed"] + rt_d)
	
	# Search isotopologues in feature matrix
	isot.match <- lapply(isot, function(x) { 
		mass_range <- c(x - ppm.s*(x/1e6), x + ppm.s*(x/1e6))
		a <- names(which(xgroup[filtsampsnoise, "mzmed"]>=mass_range[1] & xgroup[filtsampsnoise, "mzmed"]<=mass_range[2]))
		b <- names(which(xgroup[filtsampsnoise, "rtmed"]>=rt_range[1] & xgroup[filtsampsnoise, "rtmed"]<=rt_range[2]))
		r <- intersect(a,b)
		return(r)
	})
	isot.match <- unlist(isot.match)
	
	# Retrieve the feature that is the incorporation detected 
	inc <- names(compinc)[as.numeric(y)]
	
	# Vector of isotopologues
	isotv  <- c(isot.match, inc)
	
	if (length(isotv) < 2) { # Less than 2 isotopologues = no incorporation
		return()
		
	} else {
		
		iso_int <- xgroup[isotv,]
		
		inc_condition <- compinc[as.numeric(y)] 
		
		mi <- apply(iso_int[ ,3:ncol(iso_int)], 1, function(x) tapply(x, classv, mean))
		mi <- mi[ ,-which(colnames(mi) == inc)]
		
		cond <- sort(unique(unlist(strsplit(as.character(unique(inc_condition)), split=";"))))
		
		# Solve which of the isotopologues found is the base peak
		
		if (is.vector(mi)){ # If there is only 1 base peak candidate
			
			pos <- sapply(cond, USE.NAMES=F, simplify=T, function(x) {
				mi12 <- mi[grep(ULtag, names(mi))]
				mi12 <- mi12[grep(x, names(mi12))]
				if ((mi12 > Basepeak.minInt)){
					pos <- isot.match
					return(pos)
				} else {
					return()
				}				
			})
			
			pos <- unique(unlist(pos))
			
		} else { # If there is +1 base peak candidates
			
			pos <- sapply(cond, USE.NAMES=F, simplify=T, function(x) {
				mi12 <- mi[grep(ULtag, rownames(mi)), ]
				mi12 <- mi12[grep(x, rownames(mi12)), ]
				if (any(mi12 > Basepeak.minInt)) {
					mi12 <- mi12[which(mi12 > Basepeak.minInt)]
					pos <- names(which(mi12 > (Basepeak.Percentage*max(mi12))))
					return(pos)
				} else {
					return()
				}				
			})
			pos <- unlist(pos)
			
			if (is.matrix(pos)) {pos <- as.vector(pos)}
			pos <- unique(unlist(pos))
		}
		
		if (length(pos) < 1) {pos <- NULL}
		
		# If at least one feature can be considered a base peak
		if (is.null(pos)) {
			return()
			
		} else { 
			pos <- as.character(pos)
			
			if (length(pos) > 1) {# If more than one feature can be considered a base peak
				
				r <- lapply(length(pos):1,function(x) {
					pos1 <- setdiff(pos, pos[x])
					iso_int <- iso_int[which(iso_int[,"mzmed"] >= min(iso_int[pos[x],"mzmed"])), ]
					if (any(rownames(iso_int)%in%pos1)) {
						iso_int <- iso_int[-which(rownames(iso_int)%in%pos1), ]
					} else {
						iso_int <- iso_int
					}
					
					inc_id <- rep(paste(y,x,sep="."), times=nrow(iso_int))							 
					feature_id <- rownames(iso_int)
					atoms <- round(((iso_int[,"mzmed"]) - min(iso_int[,"mzmed"])) / mass_diff, digits=1)
					atoms <- sort(unname(atoms))
					
					inc_condition <- compinc[feature_id]
					inc_condition[is.na(inc_condition)] <- ""
					
					res <- cbind(as.data.frame(inc_id),as.data.frame(feature_id),iso_int[,1:2],as.data.frame(atoms),
							as.data.frame(inc_condition),iso_int[,-(1:2)])
					rownames(res) <- NULL
					return(res)
				})
				res <- do.call("rbind", r)
				rownames(res) <- NULL
				return(res)
				
			} else { # One feature can be considered a base peak
				
				iso_int <- iso_int[which(iso_int[,"mzmed"]>=min(iso_int[pos,"mzmed"])), ]
				inc_id <- rep(y, times=nrow(iso_int))							 
				feature_id <- rownames(iso_int)
				atoms <- round(((iso_int[,"mzmed"]) - min(iso_int[,"mzmed"])) / mass_diff,digits=1)
				atoms <- sort(unname(atoms))
				
				inc_condition <- compinc[feature_id]
				inc_condition[is.na(inc_condition)] <- ""
				
				res <- cbind(as.data.frame(inc_id),as.data.frame(feature_id),iso_int[,1:2],as.data.frame(atoms),
								 as.data.frame(inc_condition),iso_int[,-(1:2)])
				rownames(res) <- NULL
				return(res)
			}
		}
	}
}
)

georgedf <- do.call("rbind",george) # This should be a data.frame with 6 + Number of samples columns

##' Column 1 (inc_id): ID of the incorporation detected in the first step. 
##' Column 2 (feature_id): Row index of the feature in xcmsSet@groups 
##' Column 3 & 4 (mzmed / rtmed): Fitting values of the feature 
##' Column 5 (atoms):Number of atoms labelled in the feature 
##' Column 6 (inc_condition): Sample condition where the incorporation was detected.
##' 						If it was detected in more than one condition, they are separated by ";".
return(georgedf)
}
