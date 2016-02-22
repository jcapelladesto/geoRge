#' Calculate and compare labelling values between different conditions
#' 
#' Longer Description of label_compare
#' @param geoRgeR result of basepeak_finder
#' @param XCMSet The xcmsSet with labelled and unlabelled samples
#' @param XCMSmode="maxo" Type of quantification value for the features, as in `xcms::groupval(XCMSet)`
#' @param ULtag Unlabelled sample name tag, for example "C12"
#' @param Ltag Labelled sample name tag, for example "C13"
#' @param separator="_" separator in samplename
#' @param sep.pos position of the label/unlabelled tags in phenoData. Check `xcms::sampclass(XCMSet)`
#' @param UL.atomM Mass of the Unlabelled atom used in labelling experiments
#' @param L.atomM Mass of the Labelled atom used in labelling experiments
#' @param ppm.s ppm window to use to search the monoisotopic peak
#' @param rt.win.min Minimum retention time window in which the isotopologues are expected to coelute
#' @param Basepeak.minInt Minimum Value for the candidate monoisotopic pea
#' @param Basepeak.Percentage If more than one high-intensity base peak is found, a percentage of the highest one is used to avoid the assignation of high natural abundant peaks belonging to a monoisotopic peak 
#' @param noise.quant Expected quantile of peaks to be considered as background noise
#' @param control.cond Condition tag to be used as reference for the Welch t.test statistics
#' @param fc.vs.Control=1.0 Default fold-change value to applied to report changes between different conditions
#' @param p.value.vs.Control=0.05 Default p-value to applied to report changes between different conditions
#' @param Show.bp Boolean that switches if the monoisotopic labelling percentages should be reported in the output table
#' @return Return label_comparison
#' @export


label_compare <-
function(geoRgeR=NULL, XCMSSet=NULL, XCMSmode="maxo", PuIncR=NULL, XCMSet=NULL, ULtag=NULL, Ltag=NULL, separator="_", sep.pos=NULL,
UL.atomM=NULL, L.atomM=NULL, ppm.s=NULL, rt.win.min=1, control.cond=NULL, fc.vs.Control=1, p.value.vs.Control=0.05, Show.bp= T, ...) {

georgedf <- geoRgeR

X1 <- groupval(XCMSet, value = XCMSmode) # geoRge has only been used on "maxo". I suppose it works for others too.
D1 <- data.frame(t(X1))
colnames(D1) <- as.character(1:nrow(X1))
filtsamps <- 1:ncol(D1)
classv <- as.factor(XCMSet@phenoData$class) # sample classes (use separate folders per group when running XCMS)
xgroup <- cbind(XCMSet@groups[filtsamps,c("mzmed", "rtmed")], t(D1))

conditions <- levels(classv)
if(sep.pos=="f"){
	conditions <- gsub(paste(ULtag,separator,sep=""),"",conditions)[-grep(paste(Ltag,separator,sep=""),conditions)]
	}else{
	conditions <- gsub(paste(separator,ULtag,sep=""),"",conditions)[-grep(paste(separator,Ltag,sep=""),conditions)]
}


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
