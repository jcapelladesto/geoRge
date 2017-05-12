# geoRge

"geoRge: a computational tool for stable isotope labelling detection in LC/MS-based untargeted metabolomics"

Jordi Capellades, Miriam Navarro, Sara Samino, Marta Garcia-Ramirez, Cristina Hernandez, Rafael Simo, Maria Vinaixa and Oscar Yanes

Analytical Chemistry. 2015.

--------------------------------------------------------------------------------------------------------------------------

To Install from R console:

````
install.packages("devtools", dependencies=TRUE)
library(devtools) 

install_github("jcapelladesto/geoRge")
library(geoRge) 
````

--------------------------------------------------------------------------------------------------------------------------
To test the example data from the package you need to run the following:
````
s1 <- PuInc_seeker(XCMSet=mtbls213,ULtag="CELL_Glc12",Ltag="CELL_Glc13",
sep.pos.front=TRUE ,fc.threshold=1.5,p.value.threshold=.05,PuInc.int.lim = 4000)

s2 <- basepeak_finder(PuIncR = s1, XCMSet = mtbls213, UL.atomM=12.0,L.atomM=13.003355,
	ppm.s=6.5,Basepeak.minInt=2000)
````
--------------------------------------------------------------------------------------------------------------------------
## geoRge 1.0 

+ __Same functions, less arguments__ I have updated the functions so they are less tiring to use. Check `help(function)` to see the new arguments and function usage

+ __New documentation__ Now all functions are properly documented

+ __`label_compare()`__ has been patched

--------------------------------------------------------------------------------------------------------------------------
If you have any problem, please raise an issue ticket (https://github.com/jcapelladesto/geoRge/issues). 
