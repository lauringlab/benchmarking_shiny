# Visualization of data in *McCrone and Lauring, 2016*

This directory is a shiny app that runs a summary of our validation analysis presented in *McCrone and Lauring, 2016*. 

The large csv files that are used in this analysis are curated by git-lfs for large files. You may need to install git-lfs to get the files. Instructions for this can be found at *https://git-lfs.github.com*


## Dependencies 
## R
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
*Note : I am working to make the app even faster and perhaps host it online -JT*
