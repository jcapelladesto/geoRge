\name{geoRge-package}
\alias{geoRge-package}
\alias{geoRge}
\docType{package}
\title{
geoRge: a computational tool for stable isotope labelling detection in LC/MS-based untargeted metabolomics
}
\description{
geoRge is a package designed for analyzing untargeted LC/MS data from stable isotope-labeling experiment by comparing isotopic distributions in the mass spectra. Isotopically enriched compounds change their isotopic distribution as compared to unlabeled compounds. This is directly reflected in a number of new m/z peaks and/or higher intensity peaks in the mass spectra of labeled samples relative to the unlabeled equivalents.
}
\details{
A total of four functions are currently provided. The package also includes data sources such as \code{\link{mtbls213}}, a database in CSV format and short adducts list for both positive and negative modes. 
\describe{
 \item{\code{\link{PuInc_seeker}}}{Detect Putative Incorporation features}
 \item{\code{\link{basepeak_finder}}}{Find Monoisotopic Peaks for enriched features}
 \item{\code{\link{label_compare}}}{Calculate and Compare Labelling Values between Conditions}
 \item{\code{\link{database_query}}}{Find Monoisotopic Peaks in a Custom Metabolite Database}
}
}
\author{
Jordi Capellades, Miriam Navarro, Sara Samino, Marta Garcia-Ramirez, Cristina Hernandez, Rafael Simo, Maria Vinaixa and Oscar Yanes. (And valuable help from Steffen Neumann)

Maintainer: Jordi Capellades <j.capellades.to@gmail.com>
}
\references{
"geoRge: a computational tool for stable isotope labelling detection in LC/MS-based untargeted metabolomics"

Jordi Capellades, Miriam Navarro, Sara Samino, Marta Garcia-Ramirez, Cristina Hernandez, Rafael Simo, Maria Vinaixa and Oscar Yanes

Analytical Chemistry. 2016.
}
\keyword{ package }
