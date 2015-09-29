## server.r
require(shiny)
require(rCharts)
require(plyr)
require(ggplot2)

source("./roc_functions.R")
pat<-"([0-9]+)_([0-9]+)"

sum.df<-read.csv("./processed_data/all.removed.mapq_vcf.vcf_csv.csv",stringsAsFactors = F)

table.save<-data.frame(samp=numeric(0),mean.fp=numeric(0),mean.tp=numeric(0),p.val=numeric(0),MapQ=numeric(0),freq=numeric(0),read.range=numeric(0))

shinyServer(function(input, output) {
  
  
  data.c<-reactive({
    if(input$coding==T){
      
      sum.df<-ddply(sum.df,~chr,coding.cut)
      return(sum.df)
    }else{
      return(sum.df)
    }
  })
  
  data<-reactive({
    
    regions.trim<-mutate(regions.bed,start=start+input$trim,stop=stop-input$trim)
    print(regions.trim)
    sum.df<-data.c()
    if(input$coding==F){
      sum.df<-ddply(sum.df,~chr,function(x){
        chr<-unique(x$chr)
        start<-regions.trim$start[match(x$chr,regions.trim$chr)]
        stop<-regions.trim$stop[match(x$chr,regions.trim$chr)]
        
        subset(x,pos>start & pos<stop)
      })
      #print(head(sum.df))
    }else{
      return(sum.df)
    }
  })
  
  output$myChart <- renderChart2({
    sum.df<-data()
    cut.df<-subset(x=sum.df,MapQ>input$MapQ & gc==input$gc_roc & freq.var>input$freq.var& Read_pos<input$pos[2] & Read_pos>input$pos[1]& exp.freq %in% input$exp.freq) ##### to subset of data plotted ########
    
    roc.ls<-dlply(cut.df,~Id,sum_roc)
    
    roc.df<-roc_df(roc.ls)
    
    roc.df.adj<-ddply(roc.df,~samp,adjust.coords,cut.df)
    r_plot<-subset(roc.df.adj,is.finite(threshold))
    r_plot<-mutate(r_plot,FDR=1-adj.specificity,threshold=format(threshold,scientific=T,digits=3),TP=adj.sensitivity*20)
    
    r_plot<-mutate(r_plot, gc = paste0("10^",sub(pat,"\\1",samp)),
                   exp.freq = sub(pat,"\\2",samp)
    )
    mutate(r_plot,exp.freq=ifelse(test = grepl("0",exp.freq),yes=as.numeric(exp.freq)/1000,no = as.numeric(exp.freq)/100))->r_plot
    
    
    rp <- nPlot(adj.sensitivity~FDR, group=c("exp.freq"), data = r_plot, type = "lineChart")
    rp$params$width = 450
    rp$params$height = 450
    rp$chart(forceX = c(0,0.003))
    rp$chart(forceY = c(0,1))
    rp$chart(tooltipContent = "#! function(key, x, y, e){
             return '<b>P=</b>: ' + e.point.threshold + 
             '<b> FP</b>: '+e.point.FP+ '<b> TP</b>: '+e.point.TP
  } !#")
    rp$yAxis( axisLabel = "Sensitivity" )
    rp$xAxis( axisLabel = "1-Sensitivity")
    #rp$save("roc.html")
    rp
})
  
  
  ## remaining variants
  data.remain<-reactive({
    sum.df<-data()
    remain.df<-subset(sum.df,MapQ>input$MapQ & freq.var>input$freq.var & p.val<input$p.val & gc %in% input$gc & Read_pos<input$pos[2] & Read_pos>input$pos[1] & exp.freq %in% input$exp.freq)
    #remain.reads<-subset(reads.df,Sample %in% remain.df$Id & Mutation %in% remain.df$mutation)
  })
  output$table<-renderTable({
    remain.df<-data.remain()
    table.out<-ddply(remain.df,~gc+exp.freq,summarize,wt_mut=length(which(category=="wt")),FP=length(which(category==F)),TP=length(which(category==T)),FDR=(possible_vars-20-length(wt_mut)-FP)/(possible_vars-20-length(wt_mut)),TDR=TP/20)
    table.out},digits=3)
  
  #output$cov<-renderPlot({})
  output$samp.dis<-renderPlot({
    remain.df<-data.remain()
    ddply(remain.df,~mutation+category,summarize,count=length(mutation))->remain.count
    
    ggplot(remain.count,aes(x=count,fill=category))+geom_histogram(position="dodge")+ggtitle("distribution of variants across samples")+xlab("number of samples")+ylab("number of varinats")+scale_y_log10()
  })
  output$freq<-renderPlot({
    remain.df<-data.remain()
    ggplot(remain.df,aes(x=freq.var,fill=category))+geom_histogram(position="dodge")+ggtitle("distribution of frequency")+xlab("Frequency")+ylab("number of varinats")+scale_x_log10()
  })
  output$mean.pos<-renderPlot({
    remain.df<-data.remain()
    ggplot(remain.df,aes(x=Read_pos,fill=category))+geom_histogram(position="dodge")+ggtitle("Average read position")+xlab("Read Position")+ylab("number of varinats")
  })
  output$mean.phred<-renderPlot({
    remain.df<-data.remain()
    ggplot(remain.df,aes(x=Phred,fill=category))+geom_histogram(position="dodge")+ggtitle("Average Phred")+xlab("Phred")+ylab("number of varinats")
  })
  output$mean.mapq<-renderPlot({
    remain.df<-data.remain()
    ggplot(remain.df,aes(x=MapQ,fill=category))+geom_histogram(position="dodge")+ggtitle("Average MapQ")+xlab("MapQ")+ylab("number of varinats")
  })
  output$position<-renderPlot({
    remain.df<-data.remain()
    prior.seg.length<-c()
    
    for(k in 1:length(regions.bed$chr)){ 
      prior.seg.length[k]<-sum(regions.bed$stop[1:k])  # the end positions of each segment relative to one sequence not including the trimming step
    }
    
    prior.seg.length<-c(0,prior.seg.length)
    #remain.df<-mutate(remain.df,chr=gsub("N_A","NA",chr))
    remain.df<-mutate(remain.df, # the concatenated position relative to the usual order
                      concat.pos=pos+prior.seg.length[match(chr,regions.bed$chr)])
    
    
    FP<-subset(remain.df,category==F)
    ggplot(FP,aes(x=concat.pos))+geom_histogram(binwidth=50)+ggtitle("FP locations")+geom_line(data=regions.l,aes(x=value,y=0,col=chr),size=2)
    
    
  })
  
  
  
  
  ##### Progress graph #####
  observeEvent(input$save,{
    remain.df<-data.remain()
    table.out<-ddply(remain.df,~gc+exp.freq,summarize,wt_mut=length(which(category=="wt")),FP=length(which(category==F)),TP=length(which(category==T)),FDR=(possible_vars-20-length(wt_mut)-FP)/(possible_vars-20-length(wt_mut)),TDR=TP/20)
    
    saved<-summarize(table.out,samp=length(gc),mean.fp=mean(FP),mean.tp=mean(TP),p.val=input$p.val,MapQ=input$MapQ,freq=input$freq.var,read.range=input$pos[2]-input$pos[1])
    
    table.save[input$save,]<<-saved
    output$saved.table<-renderTable({
      
      table.save},digits=4) 
    
  })  
  })
