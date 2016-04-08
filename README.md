# Visualization of *McCrone and Lauring, 2016*

This directory is a shiny app that runs a summary of our validation analysis presentd in *McCrone and Lauring, 2016*. 

To clone this to your computer use 

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
To run the app open R and set the current directory to this one ("benchmarking_shiny")

Then just type 

```R
library(shiny)
runApp("./")
```

An html page should open with the app.

Play around with the setting to see if you can find better parameters. The defualt settings are those present in the manuscript. There are plots in the Remaining variants tab that characterizes the variants that are left after applying the cut offs.

