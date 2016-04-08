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
               h4("DeepSNV ROC"),
               showOutput("myChart", "nvd3")),
        #column(6,
        #       h4("Output table"),
        #       tableOutput('table'))
        column(6,
               h4(" Lofreq ROC"),
               showOutput("myChart2", "nvd3"))      
        ),
      
      hr(),
      fluidRow(
        column(3,
               numericInput(inputId = "MapQ",
                            label = "Mean MapQ cutoff",
                            value=-1),
               numericInput(inputId = "freq.var",
                            label = "Frequency cutoff",
                            value=0.002),
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
                                         "remove duplicates" = "no")),
               radioButtons("disp",
                            label= "DeepSNV dispersion",
                            choices=list("Binomial" = "bin", 
                                         "Betabin one sided" = "one.sided",
                                         "Betabin two sided" = "two.sided"))

        ),
        
        column(3,
               radioButtons("lofreq_dups",
                            label= " Lofreq Read Duplicates",
                            choices=list("with duplicates" = "bam", 
                                         "remove duplicates" = "removed")),
               checkboxGroupInput("exp.freq", 
                                  label = "Expected frequency for table and plots", 
                                  choices = list("5%"=0.05,"2.5%"=0.025,"1.25%"=0.0125,"0.6%"=0.0063,"0.16%"=0.0016),
                                  selected = c(0.05,0.025,0.0125,0.0063))
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
                           min = 0, max = 250, value = c(62, 188)),
               actionButton("save","Save Output"),
               submitButton(text = "Apply Changes", icon = NULL, width = NULL)
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