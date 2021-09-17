#' Find Monoisotopic Peaks for enriched features
#' 
#' Find the monoisotopic peak of the putative isotopologues detected by PuInc_seeker()
#' @param PuIncR Result of PuInc_seeker
#' @param XCMSet The xcmsSet with labelled and unlabelled samples
#' @param UL.atomM Mass of the Unlabelled atom used in labelling experiments
#' @param L.atomM Mass of the Labelled atom used in labelling experiments
#' @param ppm.s ppm window to use to search the monoisotopic peak
#' @param rt.win.min Minimum retention time window in which the isotopologues are expected to coelute
#' @param Basepeak.minInt Minimum value of intensity/Area for the candidate to be a candidate monoisotopic peak
#' @param Basepeak.Percentage If more than one high-intensity base peak is found, a percentage of the highest one is used to avoid the assignation of high natural abundant peaks belonging to a monoisotopic peak 
#' @return A list containing: \describe{
#' \item{geoRge}{Dataframe containing results for features for which a base peak was found}
#' \item{params}{Inputed arguments to be inherited by \code{\link{label_compare}}}}
#' @export
basepeak_finder <- function(PuIncR=NULL, XCMSet=NULL, UL.atomM=NULL, L.atomM=NULL, ppm.s=NULL,
             rt.win.min=1, Basepeak.minInt=NULL, Basepeak.Percentage=0.7, noise.quant=0.0) {
        
        if(!(class(PuIncR)=="list") | !(names(PuIncR)[1]=="PuInc")){
            stop("PuIncR does not match an expected PuIncR object")
        }
        res_inc <- PuIncR[["PuInc"]]
        compinc <- PuIncR[["PuInc_conditions"]]
        puinc_params <- PuIncR[["params"]]
        XCMSmode <- puinc_params$XCMSmode
        ULtag <- puinc_params$ULtag
        conditions <- puinc_params$conditions
        mass_diff <- L.atomM - UL.atomM

        xdata <- getData(XCMSet,XCMSmode)
        D1 <- xdata$D1
        classv <- xdata$classv
        featdef <- xdata$featdef
        xgroup <- cbind(featdef[,c("mzmed", "rtmed")], D1)
        
        ## Second step: Seek the base peaks for each putative incorporation #
        
        # The ppm window that is used when seeking base peaks. Must be slightly higher than the
        # mass accuracy of your instrument.
        
        # The retention time window at which the isotopologues may elute. We have
        # observed that sometimes the features are fit in a retention time window
        # narrower than expected.
        
        
        # When more than one feature is a candidate base peak.
        # A second threshold is used: a Percentage of the max mean intensity. 
        # This filters unexpected intense isotopic pattern peaks.
        
        george <- lapply(rownames(res_inc), function(y) {
            isot <- sapply(1:max_atoms(res_inc[y, "mzmed"], L.atomM), 
                           function(x) {
                               res_inc[y, "mzmed"] - (x * mass_diff)
                           })
            isot <- sort(isot)
            rt_range <- c(res_inc[y, "rtmed"] - rt.win.min,
                          res_inc[y, "rtmed"] + rt.win.min)
            b <- which(xgroup[, "rtmed"] >= rt_range[1] & 
                           xgroup[, "rtmed"] <= rt_range[2])
            isot.match <- lapply(isot, function(x) {
                mass_range <- c(x - (ppm.s * (x/1e+06)), 
                                x + (ppm.s * (x/1e+06)))
                a <- which(xgroup[, "mzmed"] >= 
                               mass_range[1] & xgroup[, "mzmed"] <= 
                               mass_range[2])
                
                r <- intersect(a, b)
                return(r)
            })
            isot.match <- rownames(xgroup)[unlist(isot.match)]
            inc_condition <- compinc[y]
            cond <- sort(unique(unlist(strsplit(as.character(unique(inc_condition)),
                                                split = ";"))))
            
            if ((length(isot.match)==0)) {return()}
            isotlog <- sapply(isot.match,function(j){
                iso_int <- xgroup[c(j,y), ]
                mi <- t(apply(iso_int[, 3:ncol(iso_int)], 1, function(x)
                    tapply(x,classv, mean)))
                mi <- mi[,grep(ULtag, colnames(mi)),drop=F]
                any(sapply(cond,function(o){
                    o <- mi[,grep(o,colnames(mi))]
                    o[1]>o[2]}))
            })
            
            isotv <- c(isot.match[which(isotlog)], y)
            if ((length(isotv) < 2)) {
                return()
            } else {
                iso_int <- xgroup[isotv, ]
                mi <- t(apply(iso_int[, 3:ncol(iso_int)], 1, function(x) 
                    tapply(x, classv, mean)))
                mi <- mi[-which(rownames(mi) == y),grep(ULtag, colnames(mi)) ,drop=F]
                pos <- lapply(cond,
                              function(x) {
                                  mi12 <- mi[,which(conditions==x),drop=T]
                                  if (any(mi12 > Basepeak.minInt)) {
                                      mibp <- which(mi12 > Basepeak.minInt)
                                      pos <- which(mi12 > 
                                                       (Basepeak.Percentage * max(mi12[mibp])))
                                      return(pos)
                                  }else {
                                      return()
                                  }
                              })
                pos <- unique(unlist(pos))
                if (length(pos) == 0) {
                    pos <- NULL
                }
                if (is.null(pos)) {
                    return()
                } else {
                    if (length(pos) > 1) {
                        r <- lapply(length(pos):1, function(x) {
                            pos1 <- setdiff(pos, pos[x])
                            iso_int <- iso_int[which(iso_int[, "mzmed"] >= 
                                                         iso_int[pos[x], "mzmed"]), ]
                            if (any(rownames(iso_int) %in% pos1)) {
                                iso_int <- iso_int[-which(rownames(iso_int) 
                                                          %in% pos1), ]
                            }
                            inc_id <- rep(paste(y, x, sep = "."), times = nrow(iso_int))
                            feature_id <- rownames(iso_int)
                            atoms <- iso_int[, "mzmed"] - iso_int[1,"mzmed"]
                            atoms <- round(atoms/mass_diff, digits = 1)
                            inc_condition <- compinc[feature_id]
                            inc_condition[is.na(inc_condition)] <- ""
                            res <- cbind(as.data.frame(inc_id), as.data.frame(feature_id), 
                                         iso_int[, 1:2], as.data.frame(atoms), 
                                         as.data.frame(inc_condition), iso_int[, 
                                                                               -(1:2)])
                            rownames(res) <- NULL
                            res <- removeDuplicates(res,mass_diff)
                            return(res)
                        })
                        res <- do.call("rbind", r)
                        rownames(res) <- NULL
                    }else {
                        iso_int <- iso_int[which(iso_int[, "mzmed"] >= 
                                                     min(iso_int[pos, "mzmed"])), ]
                        inc_id <- rep(y, times = nrow(iso_int))
                        feature_id <- rownames(iso_int)
                        atoms <- iso_int[, "mzmed"] - iso_int[1,"mzmed"]
                        atoms <- round(atoms/mass_diff, digits = 1)
                        inc_condition <- compinc[feature_id]
                        inc_condition[is.na(inc_condition)] <- ""
                        res <- cbind(as.data.frame(inc_id), as.data.frame(feature_id), 
                                     iso_int[, 1:2], as.data.frame(atoms), as.data.frame(inc_condition), 
                                     iso_int[, -(1:2)])
                        rownames(res) <- NULL
                        res <- removeDuplicates(res,mass_diff)
                    }
                    return(res)
                }
            }
        })
        georgedf <- do.call("rbind", george)
  
        ##Column 1 (inc_id): ID of the incorporation detected in the first step. 
        ##Column 2 (feature_id): Row index of the feature in xcmsSet-groups
        ##Column 3 & 4 (mzmed / rtmed): Fitting values of the feature 
        ##Column 5 (atoms):Number of atoms labelled in the feature 
        ##Column 6 (inc_condition): Sample condition where the incorporation was detected.
        ##						If it was detected in more than one condition, they are separated by ";".
        ##					
        
        puinc_params <- append(puinc_params, values = list(UL.atomM = UL.atomM, 
                                                           L.atomM = L.atomM))
        return(list("geoRge" = georgedf,
                    "params" = puinc_params))
    }


#' Refine basepeak_finder Results
#' 
#' Performs a cleaning step of geoRge results based on number of isotopologues and correlation
#' @param PuIncR Result of basepeak_finder
#' @param minAtoms Minimum number of atoms of an inc_id cluster
#' @param maxAtoms Maximum number of atoms of an inc_id cluster
#' @param minCorr Minimum \link{stats::cor} correlation of M0 and Mn isotopologues along labelled samples.
#' @note Depending on the number of samples, the correlation value may be oversensitive, be careful with this argument.
#' @return A refined geoRge data.frame object. 
#' @export
refinegeoRge <- function(finderR, minAtoms=3, maxAtoms=50, minCorr=0.25){
    geoRge <- finderR$geoRge
    L.atomM <- finderR$params$L.atomM
    UL.atomM <- finderR$params$UL.atomM
    Ltag <- finderR$params$Ltag
    mass_diff <- L.atomM - UL.atomM
    
    u_inc_id <- unique(geoRge$inc_id)
    res <- lapply(u_inc_id,function(x){
        r <- geoRge[which(geoRge$inc_id==x),]
        ratoms <- r$atoms[nrow(r)]
        if(ratoms>maxAtoms){return()}
        if(ratoms<minAtoms){return()}
        cv <- sapply(strsplit(r$inc_condition[nrow(r)],";")[[1]],function(i){
            a0 <- as.numeric(r[1,intersect(grep(i,colnames(r)),
                                           grep(Ltag,colnames(r)))])
            aI <- as.numeric(r[nrow(r),intersect(grep(i,colnames(r)),
                                                 grep(Ltag,colnames(r)))]) 
            cor(a0,aI)
        })
        if(any(cv>minCorr)){
            return(r)
        }else{
            return()
        }
    })
    res <- do.call("rbind",res)
    
    u_inc_id <- unique(res$feature_id[which(res$atoms==0)])
    res <- lapply(u_inc_id,function(x){
        # print(x)
        x <- unique(res$inc_id[which(res$feature_id==x)])
        if(length(x)>1){
            # icount <- sapply(x, function(i) length(which(res$inc_id==i)))
            # x <- x[which(icount==max(icount))]
            # if(length(x)>1){ x <- x[length(x)]}
            u_feat <- unique(unlist(res$feature_id[res$inc_id%in%x]))
            u_feat <- sort(u_feat)
            idx <- sapply(u_feat,function(j){which(res$feature_id==j)[1]})
            idx <- unique(unlist(idx))
            res <- res[idx,]
            res$inc_id <- res$feature_id[nrow(res)]
            res <- removeDuplicates(res, mass_diff)
            # res$inc_id <- res$feature_id[nrow(res)]
            return(res)
        }else{
            return(res[which(res$inc_id==x),])
        }
    })
    inc_ids <- sapply(res,function(x) x$inc_id[1])
    dup_ids <- which((duplicated(inc_ids)+duplicated(inc_ids,fromLast = T))>0)
    dupdone <- c()
    for(i in dup_ids){
        dup <- inc_ids[i]
        if(dup%in%dupdone){next}
        dupidx <- which(inc_ids==dup)
        for(j in seq_along(dupidx)){
            r <- res[[dupidx[j]]]
            r$inc_id <- paste(r$inc_id[1],j, sep = ".")
            res[[dupidx[j]]] <- r
            dupdone <- append(dupdone,dup)
        }
    }
    res <- do.call("rbind",res)
    rownames(res) <- NULL
    return(list("geoRge" = res,
                "params" = finderR$params))
}

# Estimate max number of atoms for a given mass
max_atoms <- function(mass, uatommass) {
    # mass = incorporation mass
    # uatommass = mass of the unlabeled atom
    floor((mass/uatommass))
}

# Remove duplicated isotopologues in geoRgedf list elements
removeDuplicates <- function(res,mass_diff){
    if (any(duplicated(res$atoms))) {
        dup <- unique(res$atoms[which(duplicated(res$atoms))])
        for (i in dup) {
            iix <- which(res$atoms == i)
            temp <- res[iix, ]
            if (sum(temp$inc_condition != "") == 1) {
                j <- which(temp$inc_condition == "")
                j <- iix[j]
                res <- res[-j, ]
            }  else {
                mz0 <- (res$mzmed[1] + (mass_diff * i))
                ppmerr <- abs(temp$mzmed - mz0) * 1e+06/mz0
                j <- (1:length(ppmerr))[-which.min(ppmerr)]
                j <- iix[j]
                res <- res[-j, ]
            }
        }
    }
    return(res)
}

