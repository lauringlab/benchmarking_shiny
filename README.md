# benchmarking_shiny
A shiny app that runs a summary of the benchmarking analysis.  The data used includes duplicate reads (which is legit because of our super deep coverage) and estimates dispersion using a two sided alternative hypothesis due to the high level of expected PCR interference.

To clone this to your computer use 

```
git clone https://github.com/lauringlab/benchmarking_shiny.git
```
you should be prompted for your username and password.

To run the app you will need several R libraries.  Some are tricky and don't come from CRAN. To install them in R use

```R
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("devtools")

require(devtools)
install_github('rCharts', 'ramnathv’)

```
The others are easy.  I they should be 

```R
install.packages("pROC")
install.packages("ggplot2")
install.packages("plyr")
install.packages("reshape2")
```
If I forgot one in my haste and you get an error saying that a package was not found please install it and try again.

To run the app open R and set the current directory to this one ("benchmarking_shiny")

Then just type 

```R
library(shiny)
runApp("./")
```

a web page should open with the app.

Play around with the setting to see if you can find better parameters.  When you find ones you like click the save output button.  It will log the results in a table on the progress tab.  For now I'm just averaging TP and FP accross all samples in the detailed table on the front. If you think a sample is schewing the result you can remove it from the check boxes on the overiew tab.  Also there are plots in the Remaining variants tab that characterizes the variants that are left after applying the cut offs to help set more stringent ones if you would like.

