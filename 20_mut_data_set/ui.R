## ui.R
require(shiny)
require(rCharts)
require(plyr)
require(ggplot2)
shinyUI(fluidPage(
  titlePanel("Benchmarking variant caller"),
  
  tabsetPanel(
    tabPanel("Overview",
      fluidRow(
        column(6,
               h4("ROC"),
               showOutput("myChart", "nvd3")),
        column(6,
               h4("Output table"),
               tableOutput('table'))
      ),
      
      hr(),
      fluidRow(
        column(3,
               sliderInput(inputId = "MapQ",
                            label = "Mean MapQ cutoff",
                            min=0,max=42,value=0),
               sliderInput(inputId = "Phred",
                           label = "Mean Phred cutoff",
                           min=30,max=42,value=30),
               
               numericInput(inputId = "freq.var",
                            label = "Frequency cutoff",
                            value=0),
               numericInput(inputId = "p.val",
                            label = "p value cutoff",
                            value=0.01)
        ),
        column(3,
               radioButtons("method",
                            label= "P value correction",
                            choices=list("Bonferroni" = "bon", 
                                         "BH" = "BH"),
                            selected='bon'),
               radioButtons("dups",
                            label= " Read Duplicates",
                            choices=list("with duplicates" = "with", 
                                         "remove duplicates" = "no"),
                            selected="no"),
               radioButtons("disp",
                            label= "DeepSNV dispersion",
                            choices=list("Binomial" = "bin", 
                                         "Betabin one sided" = "one.sided",
                                         "Betabin two sided" = "two.sided"),
                            selected="two.sided")
               
        ),
        
        column(3,
               radioButtons("gc_roc",
                            label= "Genome copy input for ROC",
                            choices=list("10^5" = 5, 
                                         "10^4" = 4, "10^3" = 3),
                            selected = 5),
               checkboxGroupInput("gc", 
                                  label = "Genome copy table and plots", 
                                  choices = list("10^5" = 5, 
                                                 "10^4" = 4, "10^3" = 3),
                                  selected = 5),
               checkboxGroupInput("exp.freq", 
                                  label = "Expected frequency for table and plots", 
                                  choices = list("5%"=0.05,"2%"=0.02,"1%"=0.01,"0.5%"=0.005,"0.2%"=0.002),
                                  selected = c(0.05,0.02,0.01,0.005,0.002))
        ),
        column(3,
               checkboxInput("coding",
                             label = "Only Coding region",
                             value=F),
               numericInput(inputId = "trim",
                            label = "Trim segment ends",
                            value=0),
               sliderInput("pos",
                           label="Read position cut off",
                           min = 0, max = 125, value = c(0, 125)),
               actionButton("save","Save Output"),
               submitButton(text = "Apply Changes", icon = NULL, width = NULL)
        )
      )
    ),
    tabPanel("Remaining Variants",
             fluidRow(
               column(12,
                      #plotOutput("samp.dis"),
                      plotOutput("freq"),
                      plotOutput("mean.pos"),
                      plotOutput("mean.qual")
                      #plotOutput("mean.mapq"),
                      #plotOutput("position")
               )
             )  
    )
  )  
))