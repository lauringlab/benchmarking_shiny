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
               numericInput(inputId = "MapQ",
                            label = "Mean MapQ cutoff",
                            value=20),
               numericInput(inputId = "freq.var",
                            label = "Frequency cutoff",
                            value=0.002)
        ),
        column(3,
               radioButtons("method",
                            label= "P value correction",
                            choices=list("Bonferroni" = "bon", 
                                         "BH" = "BH"),
                            selected='bon'),
               numericInput(inputId = "p.val",
                            label = "p value cutoff",
                            value=0.01)
        ),
        
        column(3,
               radioButtons("gc_roc",
                            label= "Genome copy input for ROC",
                            choices=list("10^5" = "5", 
                                         "10^4" = "4", "10^3" = "3"),
                            selected = "5"),
               checkboxGroupInput("gc", 
                                  label = "Genome copy table and plots", 
                                  choices = list("10^5" = "5", 
                                                 "10^4" = "4", "10^3" = "3","10^2"="2"),
                                  selected = c("5","4")),
               checkboxGroupInput("exp.freq", 
                                  label = "Expected frequency for table and plots", 
                                  choices = list("5%"=0.05,"2%"=0.02,"1%"=0.01,"0.5%"=0.005,"0.2%"=0.002,"0.1%"=0.001),
                                  selected = c(0.05,0.02,0.01,0.005))
               ),
        column(3,
               radioButtons("q",
                            label="Phred cut off",
                            choices=list("25"="25","30"="30"),
                            selected="25"),
               sliderInput("pos",
                            label="Read position cut off",
                           min = 0, max = 125, value = c(32, 94)),
               actionButton("save","Save Output")
        )
      )
    ),
    tabPanel("Remaining Variants",
      fluidRow(
        column(12,
               plotOutput("samp.dis"),
               plotOutput("freq"),
               plotOutput("mean.pos"),
               plotOutput("mean.phred"),
               plotOutput("mean.mapq"),
               plotOutput("position")
               )
      )  
    ),
    tabPanel("Progess",
        fluidRow(
          column(12,
                 tableOutput("saved.table")
          )
        )  
    )
  )  
))