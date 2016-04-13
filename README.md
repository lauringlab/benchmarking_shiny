# Visualization of the data presented in "Measurements of intrahost viral diversity are extremely sensitive to systematic errors in variant calling"

This directory is a shiny app that runs a summary of our validation analysis presented in *McCrone and Lauring, 2016*. 

One of the draw backs to the shiny app is that it can take a litte work to get it up and running. We have made a lighter version of the app that is hosted here *http://lauringlab.github.io/benchmarking_shiny/*. This version of the app only deals with the data presented in the paper, it does not provide means for changing as many variables as the shiny app. Give it a look. Currently each plot can only get 500 views a day. I'm not sure if that will be a problem but I'm working on a work around. Also please note that only variants with p<0.01 are included in the app so the ROC may appear slightly truncated compared to what is seen in the paper figures.

## Dependencies 
The large csv files that are used in this analysis are curated by git-lfs for large files. You may need to install git-lfs to get the files. Instructions for this can be found at *https://git-lfs.github.com*.

Once you have git-lfs installed and have cloned or downloaded this repository you will need to get the large data files that are used by the app. You can do this by running the following commands from inside the repository directory.

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
*Note : I am working to make the app even faster and perhaps host it online. Also the save output button doesn't work yet. -JT*
