\name{mtbls213}
\alias{mtbls213}
\docType{data}
\title{xcmsSet for MTBLS213
}
\description{
LC/MS data from a human retinal pigment epithelium cell line (ARPE-19) grown on normal and high glucose concentrations that mimics diabetic retinopathy conditions in vitro.
}
\usage{data("mtbls213")}
\format{
  The format is:
Formal class 'xcmsSet' [package "xcms"] with 14 slots
  ..@ peaks           : num [1:205594, 1:11] 169 169 169 169 173 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:11] "mz" "mzmin" "mzmax" "rt" ...
  ..@ groups          : num [1:14607, 1:11] 100 100 100 100 101 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:11] "mzmed" "mzmin" "mzmax" "rtmed" ...
  ..@ groupidx        :List of 14607
  .. ..$ : int [1:12] 4940 14904 26709 36721 67803 77949 88441 152196 159043 186197 ...
[...]
  .. ..$ : int [1:12] 5704 16678 27596 37577 48785 58589 68787 79316 90029 101538 ...
  .. .. [list output truncated]
  ..@ filled          : int [1:80267] 125328 125329 125330 125331 125332 125333 125334 125335 125336 125337 ...
  ..@ phenoData       :'data.frame':	12 obs. of  1 variable:
  .. ..$ class: Factor w/ 4 levels "CELL_Glc12_05mM_Normo",..: 1 1 1 2 2 2 3 3 3 4 ...
  ..@ rt              :List of 2
  .. ..$ raw      :List of 12
  .. .. ..$ : num [1:5023] 2.23 2.48 2.73 2.98 3.23 ...
[...]
  .. .. ..$ : num [1:5024] 2.15 2.4 2.65 2.9 3.15 ...
  .. ..$ corrected:List of 12
  .. .. ..$ : num [1:5023] 2.21 2.46 2.72 2.97 3.22 ...
[...]
  .. .. ..$ : num [1:5024] 2.21 2.46 2.72 2.97 3.22 ...
  ..@ filepaths       : chr [1:12] "/home/sneumann/src/R/geoRgeProject/MTBLS213/CELL_Glc12_05mM_Normo_04.mzXML" "/home/sneumann/src/R/geoRgeProject/MTBLS213/CELL_Glc12_05mM_Normo_05.mzXML" "/home/sneumann/src/R/geoRgeProject/MTBLS213/CELL_Glc12_05mM_Normo_06.mzXML" "/home/sneumann/src/R/geoRgeProject/MTBLS213/CELL_Glc12_25mM_Normo_16.mzXML" ...
  ..@ profinfo        :List of 2
  .. ..$ method: chr "bin"
  .. ..$ step  : num 0.1
  ..@ dataCorrection  : int(0) 
  ..@ polarity        : chr(0) 
  ..@ progressInfo    :List of 12
  .. ..$ group.density          : num 0
  .. ..$ group.mzClust          : num 0
  .. ..$ group.nearest          : num 0
  .. ..$ findPeaks.centWave     : num 0
  .. ..$ findPeaks.massifquant  : num 0
  .. ..$ findPeaks.matchedFilter: num 0
  .. ..$ findPeaks.MS1          : num 0
  .. ..$ findPeaks.MSW          : num 0
  .. ..$ retcor.obiwarp         : num 1
  .. ..$ retcor.peakgroups      : num 0
  .. ..$ fillPeaks.chrom        : num 0
  .. ..$ fillPeaks.MSW          : num 0
  ..@ progressCallback:function (progress)  
  ..@ mslevel         : num(0) 
  ..@ scanrange       : num(0) 
}
\details{
Obtained by 

setwd("../MTBLS213/")
    
    xset <- xcmsSet(method="centWave", ppm=30, peakwidth=c(5,20))
    sampclass(xset) <- substring(sampnames(xset), 1, 21)

    xset2 <- retcor(xset,method="obiwarp", profStep=0.1) 
    xset3 <- group(xset2, mzwid=0.0065, minfrac=0.5, bw= 4)
    xset3 <- fillPeaks(xset3)

}
\source{
http://www.ebi.ac.uk/metabolights/MTBLS213
}
\references{
"geoRge: a computational tool for stable isotope labelling detection in LC/MS-based untargeted metabolomics"

Jordi Capellades, Miriam Navarro, Sara Samino, Marta Garcia-Ramirez, Cristina Hernandez, Rafael Simo, Maria Vinaixa and Oscar Yanes

Analytical Chemistry. 2015.
}
\examples{
data(mtbls213)
xset3
}
\keyword{datasets}
