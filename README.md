# Visualization of *McCrone and Lauring, 2016*

This directory is a shiny app that runs a summary of our validation analysis presented in *McCrone and Lauring, 2016*. 

The large csv files that are used in this analysis are curated by git-lfs for large files. You may need to install git-lfs to get the files. Instructions for this can be found at *https://git-lfs.github.com*

```
git clone https://github.com/lauringlab/benchmarking_shiny.git
```


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
To run the app open R and set the current directory to this one ("benchmarking_shiny").

Two datasets were used in the paper. In one two influenza strains were mixed (PR8 and WSN33) to investigate the data of this data set in R run the following commands.

```R
library(shiny)
runApp("./PR8-WSN33")
```

An html page should open with the app. For several reasons we didn't tinker with this data set too much : 1) PR8 and WSN33 are more diverse than normally seen in patient samples. 2) The titer of these samples was much larger than that seen in patient samples 3) The samples were made in a manner that avoided and PCR biasing.

In the more interesting data set we mixed 20 point mutants in wt at differing frequencies and titers. To investigate this data set use the following commands in R.
```R
require(shiny)
runApp("./20_mut_data_set")
```


Play around with the setting to see if you can find better parameters. The defualt settings are those present in the manuscript. There are plots in the Remaining variants tab that characterizes the variants that are left after applying the cut offs.

