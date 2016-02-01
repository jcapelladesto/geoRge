## Example script for geoRge usage

## Check our metabolights repository for the mzXML files used in the example below (http://www.ebi.ac.uk/metabolights/MTBLS213)
## and place it in a directory "../MTBLS213/"
## Use XCMS package (https://bioconductor.org/packages/release/bioc/html/xcms.html) for peak picking, alignment and grouping

library(xcms)

recalc=FALSE

if (recalc) {
    setwd("../MTBLS213/")
    
    xset <- xcmsSet(method="centWave", ppm=30, peakwidth=c(5,20))
    sampclass(xset) <- substring(sampnames(xset), 1, 21)

    xset2 <- retcor(xset,method="obiwarp", profStep=0.1) 
    xset3 <- group(xset2, mzwid=0.0065, minfrac=0.5, bw= 4)
    xset3 <- fillPeaks(xset3)
    setwd("../geoRge/")

    save(xset3, file="mtbls213.Rda")

} else {
    load("mtbls213.Rda")
}

## Parameters such as ppm, peakwidth, mzwid, minfrac and bw are not strict. They depend on data acquisition and methodology.
## Please find the best parameters for your own data using your own experience or by trial and error

source("george.R")

s1 <- PuInc_seeker(XCMSet=xset3,ULtag="CELL_Glc12",Ltag="CELL_Glc13",sep.pos="f")

s2 <- basepeak_finder(PuIncR=s1,XCMSet=xset3,ULtag="CELL_Glc12",Ltag="CELL_Glc13",
	sep.pos="f",UL.atomM=12.0,L.atomM=13.003355,
	ppm.s=6.5,Basepeak.minInt=2000)

negative <- read.table("./adducts_negative.txt",header=T,stringsAsFactors=F)
db <- read.csv("./ExampleDatabase.csv",header=T,stringsAsFactors=F,fill=T)
hits <- database_query(geoRgeR = s2, adducts = negative, db = db)



