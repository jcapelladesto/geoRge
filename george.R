PuInc_seeker <- function (XCMSet=NULL, XCMSmode="maxo", ULtag=NULL, Ltag=NULL, separator="_", sep.pos=NULL, fc.threshold=1.2, 
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

basepeak_finder <- function(PuIncR=NULL, XCMSmode="maxo", XCMSet=NULL, ULtag=NULL, Ltag=NULL, separator="_", sep.pos=NULL,
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


#### Calculate the percentage of incorporation and compare it versus a "control" condition for each "inc_id" ####
label_compare <- function(geoRgeR=NULL, XCMSmode="maxo", PuIncR=NULL, XCMSet=NULL, ULtag=NULL, Ltag=NULL, separator="_", sep.pos=NULL,
UL.atomM=NULL, L.atomM=NULL, ppm.s=NULL, rt.win.min=1, control.cond=NULL, fc.vs.Control=1, p.value.vs.Control=0.05, Show.bp= T, ...) {

georgedf <- geoRgeR

X1 <- groupval(XCMSet, value = XCMSmode) # geoRge has only been used on "maxo". I suppose it works for others too.
D1 <- data.frame(t(X1))
colnames(D1) <- as.character(1:nrow(X1))
filtsamps <- 1:ncol(D1)
classv <- as.factor(XCMSet@phenoData$class) # sample classes (use separate folders per group when running XCMS)
xgroup <- cbind(XCMSet@groups[filtsamps,c("mzmed", "rtmed")], t(D1)) 

conditions <- sapply(levels(classv), USE.NAMES=F, function(x) {
	if (sep.pos=="f") {
	a <- strsplit(x, split = paste(ULtag, separator, sep=""))[[1]][2]
	} else {
	a <- strsplit(x, split = paste(separator, ULtag, sep=""))[[1]][1]
	}
	return (a)
})
conditions <- unique(conditions[-which(is.na(conditions))])

mass_diff <- L.atomM - UL.atomM

percent.incorp <- lapply(unique(georgedf$inc_id), function(y) {
	
	inc_id_features <- georgedf[which(georgedf$inc_id==y), ] 
	inc_id_int <- inc_id_features[ ,7:ncol(inc_id_features)]
	
	rts <- inc_id_features$rtmed  
	rt_range <- c(min(rts), max(rts))
	
	inc_isot <- max(inc_id_features$atoms)+1
	isot_m <- inc_id_features[1, "mzmed"] + (inc_isot * mass_diff) # mass of isotope of the incorporation
	
	isot_id <- lapply(isot_m,function(x) { # seek isotope of the incorporation in raw data
		mass_range <- c(x - ppm.s*(x/1e6), x + ppm.s*(x / 1e6))
		a <- which(xgroup[,"mzmed"] >= mass_range[1] & xgroup[,"mzmed"] <= mass_range[2])  
		b <- which(xgroup[,"rtmed"] >= rt_range[1] & xgroup[,"rtmed"] <= rt_range[2])
		r <- intersect(a,b)
		return(r)
	})
	
	isot_id <- unlist(isot_id)
	
	if(length(isot_id)<1) {
		all_id <- xgroup[as.character(inc_id_features$feature_id), ] 
	} else {
		isot_id <- isot_id[1]
		all_id <- c(as.character(inc_id_features$feature_id), isot_id)
		all_id <- xgroup[all_id, ]  
	}
	
	all_id_int <- all_id[ ,3:ncol(all_id)] # intensities
	
	inc_percent <- sapply(conditions,function(x){
		inc_id_intL <- inc_id_int[,intersect(grep(Ltag,colnames(inc_id_int)),grep(x,colnames(inc_id_int)))]
		all_id_intL <- all_id_int[,intersect(grep(Ltag,colnames(all_id_int)),grep(x,colnames(all_id_int)))]
		
		inc_cal <- sapply(1:ncol(inc_id_intL), function(x) {(inc_id_intL[ ,x] / sum(all_id_intL[ ,x]))*100})
		return(inc_cal)
	})
	colnames(inc_percent) <- conditions
	
	atoms <- inc_id_features$atoms
	rownames(inc_percent) <- rep(atoms, length.out=nrow(inc_percent)) 
	
	mean_inc <- sapply(conditions, USE.NAMES=T, simplify=T, function(x) {
		sapply(atoms, function(z) {
			inc_p_v <- inc_percent[which(rownames(inc_percent) == z), x]
			inc_p_m <- mean(inc_p_v)  
			return(inc_p_m)
		})
	})
	colnames(mean_inc) <- paste0(colnames(mean_inc),"_MEAN")
	
	sd_inc <- sapply(conditions, USE.NAMES=T, simplify=T, function(x) {
		sapply(atoms,function (z) {
			inc_p_v <- inc_percent[which(rownames(inc_percent)==z),x]
			inc_p_s <- sd(inc_p_v)  
			return(inc_p_s)
		})
	})
	colnames(sd_inc) <- paste0(colnames(sd_inc), "_SD")
	
	noncontrol <- setdiff(conditions, control.cond)
	
	pvals <- sapply(noncontrol, function (x) {
		sapply(atoms, function(y) {
			a <- try(t.test(inc_percent[which(rownames(inc_percent) == y), which(conditions == x)],
					inc_percent[which(rownames(inc_percent) == y), which(conditions == control.cond)],
					var.equal=T)$p.value,silent=T)
			if(is(a,"try-error")) {a <- 1}
			return(a)
		})
	})
	
	fct <- sapply(noncontrol, simplify=T, function (x) {
		sapply(1:nrow(mean_inc), function (y) { 
			case <- mean_inc[y, which(conditions == x)]
			control <- mean_inc[y, which(conditions == control.cond)]
			FC <- case/control
			FC2 <- (-(control/case))
			FC[FC<1] <- FC2[FC<1]
			names(FC) <- NULL
			return(FC)
		})
	})
	
	comp <- sapply(1:nrow(pvals),function(x) {
		t <- which(pvals [x, ] < p.value.vs.Control)
		if(length(t) == 0) {
		t <- ""
		return(t)
	} else {
		up <- which(fct[x, names(t)] > fc.vs.Control)
		down <- which(fct[x, names(t)] < (-fc.vs.Control))
		if(length(up)!=0) {
		names(t)[up] <- paste("UP", names(t)[up], sep="_")
		}
		if(length(down)!=0) {
		names(t)[down] <- paste("DOWN", names(t)[down], sep="_")
		}
		return(names(t))
		}
	})
	
	if(is.matrix(comp)) {
		comp <- sapply(1:ncol(comp), function(x) {
			a <- paste(comp[ ,x], collapse=";")
		})
	} else {
		comp <- lapply(1:length(comp), function(x) paste(comp[[x]], collapse=";"))
		comp <- unlist(comp)
	}
	
	comp[1] <- "Base peak"
	
	r <- data.frame("Comparison" = comp, mean_inc, sd_inc)
	if (!Show.bp) {r[1, ] <- rep("Base Peak", times=ncol(r))}
	return(r)
})

percent.incorpdf <- do.call("rbind",percent.incorp)

return(percent.incorpdf)

}

database_query <- function(geoRgeR=NULL, adducts=NULL, db=NULL, ppm.db=10) {

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
