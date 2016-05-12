# Visualization of the data presented in "Measurements of intrahost viral diversity are extremely sensitive to systematic errors in variant calling"

This directory contains a Shiny application and web page that runs a summary of our validation analysis presented in *McCrone and Lauring, 2016*. 

The Shiny app is run in R, and allows the user to choose serval variables that affect variant calling :
* Whether or not to include pcr duplicate reads
* Which p value correction method to use (Bonferroni or Benjamini Hochberg)
* Should the DeepSNV error model be a binomial model or a betabinomial.
* If betabinomial should the test be onesided or twosided (two sided is recommend with large amounts of PCR amiplication) 
* What p value threshold to use
* What mapping qualilty threshold to apply (if any)
* What Phred score threshold to apply (if any)
* What read position threshold should be applied (if any)
* What frequency threshold to apply (if any)
* Should only the coding regions be included in analysis
* Should the ends of the genomic segements be excluded from analysis (if so how much should be trimmed)



One of the draw backs to the shiny app is that it can take a litte work to get it up and running. We have made a lighter version of the app that is hosted here http://lauringlab.github.io/benchmarking_shiny/. This version of only deals with the data presented in the paper. 
Here you can choose : 
* What mapping qualilty threshold to apply (if any)
* What Phred score threshold to apply (if any)
* What read position threshold should be applied (if any)

Feel free to try it out before downloading the Shiny version.  Please note that only variants with p<0.01 are included in the app so the ROC may appear slightly truncated compared to what is seen in the paper figures.


## Dependencies 
The large csv files that are used in this analysis are curated by git-lfs for large files. You may need to install git-lfs to get the files. Instructions for this can be found at https://git-lfs.github.com.

Once you have git-lfs installed and have cloned this repository you will need to get the large data files that are used by the app. You can do this by running the following commands from inside the repository directory.

```bash
git lfs fetch
git lfs checkout 
```
The files should be downloaded and updated.

### R
* Biostrings
* rCharts
* pROC
* ggplot2
* reshape2
* shiny

See below for help installing.

To run the app open R and set the current directory to this one ("benchmarking_shiny").

To investigate the data used in Figures 3,4 and 5 data set use the following commands in R.
```R
require(shiny)
runApp("./20_mut_data_set")
```


Play around with the setting to see if you can find better parameters. The defualt settings are those present in the manuscript. There are plots in the Remaining variants tab that characterizes the variants that are left after applying the cut offs.


### R dependencies
To run the app you will need several R libraries.  Some are tricky and don't come from CRAN. To install them in R use

```R
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("devtools")

require(devtools)
install_github("rCharts", "ramnathv")

```
The others are easy.  I they should be 

```R
install.packages("pROC")
install.packages("ggplot2")
install.packages("plyr")
install.packages("reshape2")
```
*Note : The save output button doesn't work yet. I'll remove it soon. -JT*

The app has been tested and works under the following sessionInfo()
```
sessionInfo()
R version 3.2.3 (2015-12-10)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] reshape2_1.4.1      Biostrings_2.36.4   XVector_0.8.0      
 [4] IRanges_2.2.9       S4Vectors_0.6.6     BiocGenerics_0.14.0
 [7] pROC_1.8            ggplot2_2.1.0       plyr_1.8.3         
[10] rCharts_0.4.5       shiny_0.13.0       

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.4      magrittr_1.5     whisker_0.3-2    zlibbioc_1.14.0 
 [5] munsell_0.4.3    colorspace_1.2-6 xtable_1.8-2     lattice_0.20-33 
 [9] R6_2.1.2         stringr_1.0.0    tools_3.2.3      grid_3.2.3      
[13] gtable_0.2.0     htmltools_0.3.5  yaml_2.1.13      digest_0.6.9    
[17] RJSONIO_1.3-0    mime_0.4         labeling_0.3     stringi_1.0-1   
[21] scales_0.4.0     jsonlite_0.9.19  httpuv_1.3.3    
```
