## server.r
require(shiny)
require(rCharts)
require(plyr)
require(ggplot2)

source("~/Desktop/benchmarking_shiny/20_mut_data_set/roc_functions.R")
pat<-"([0-9]+)_([0-9]+)"


shinyServer(function(input, output) {

  
  dataInput<-reactive({  #### changes to data used ##########
    criteria<-paste0(input$dups,".*",input$method,".*",input$disp,".*","sum.csv$")
    file<-list.files(path = "~/Desktop/benchmarking_shiny/20_mut_data_set/processed_data",pattern = criteria,full.names = T)
    print(file)

    sum.df.raw<-read.csv(file,comment.char='#',stringsAsFactors = F)
    print(names(sum.df.raw))
    return(sum.df.raw)
  })
  
  data.wt<-reactive({
    criteria<-paste0(input$dups,".*",input$method,".*",input$disp,".*","wt.csv$")
    file<-list.files(path = "~/Desktop/benchmarking_shiny/20_mut_data_set/processed_data",pattern = criteria,full.names = T)
    print(file)

    wt<-read.csv(file,comment.char='#',stringsAsFactors = F)
    return(wt$wt)
  })
  
  data.c<-reactive({
    sum.df<-dataInput()
    if(input$coding==T){
      
      sum.df<-ddply(sum.df,~chr,coding.cut)
      return(sum.df)
    }else{
      return(sum.df)
    }
  })
  possible.vars<-reactive({
    wt_mut<-data.wt()
    poss_vars(mfe.df,wt_mut)
  })
  data<-reactive({
    regions.trim<-mutate(regions.bed,start=start+input$trim,stop=stop-input$trim)
    sum.df<-data.c()
    if(input$coding==F){
      sum.df<-ddply(sum.df,~chr,function(x){
        chr<-unique(x$chr)
        start<-regions.trim$start[match(x$chr,regions.trim$chr)]
        stop<-regions.trim$stop[match(x$chr,regions.trim$chr)]
        
        subset(x,pos>start & pos<stop)
      })
    }else{
      return(sum.df)
    }
  })
  
  output$myChart <- renderChart2({
    #print(total_possible_vars)
    sum.df<-data()
    
    possible_vars<-possible.vars()
    #print(possible_vars)
    cut.df<-subset(x=sum.df,MapQ>input$MapQ & Phred>input$Phred & gc==input$gc_roc & freq.var>input$freq.var& Read_pos<=input$pos[2] & Read_pos>=input$pos[1]& exp.freq %in% input$exp.freq) ##### to subset of data plotted ########
   #print(head(cut.df))
    h.roc.df<-hiseq.roc(cut.df,20,possible_vars,">")
    r_plot<-h.roc.df[order(h.roc.df$adj.specificity,decreasing = T),]

    
    r_plot<-mutate(r_plot,FDR=1-adj.specificity,threshold=format(threshold,scientific=T,digits=3))
    
    r_plot<-mutate(r_plot, gc = paste0("10^",sub(pat,"\\1",samp)))
    print(head(r_plot))
    
    
    #print(r_plot)
    rp <- nPlot(adj.sensitivity~FDR, group=c("exp.freq"), data = r_plot, type = "lineChart")
    rp$params$width = 450
    rp$params$height = 450
    rp$chart(forceX = c(0,0.005))
    rp$chart(forceY = c(0,1))
    rp$chart(tooltipContent = "#! function(key, x, y, e){
  return '<b>P=</b>: ' + e.point.threshold + 
  '<b> FP</b>: '+e.point.FP+ '<b> TP</b>: '+e.point.TP
} !#")
    rp$yAxis( axisLabel = "Sensitivity" )
    rp$xAxis( axisLabel = "1-Specificity")
    #rp$save("roc.html")
    rp
  })

  
  ## remaining variants
  data.remain<-reactive({
    sum.df<-data()
    
    remain.df<-subset(sum.df,MapQ>=input$MapQ & Phred>=input$Phred & freq.var>=input$freq.var & p.val<input$p.val & gc %in% input$gc & Read_pos<=input$pos[2] & Read_pos>=input$pos[1] & exp.freq %in% input$exp.freq & category %in% c(T,F))
    #remain.reads<-subset(reads.df,Sample %in% remain.df$Id & Mutation %in% remain.df$mutation)

  })
  output$table<-renderTable({
    remain.df<-data.remain()
    possible_vars<-possible.vars()
    print(possible_vars)
    remain.df$possible_vars<-possible_vars
    table.out<-ddply(remain.df,~gc+exp.freq,summarize,Sens =length(which(category==T))/20,TP=length(which(category==T)), Spec=(possible_vars[1]-length(which(category==F)))/(possible_vars[1]),FP=length(which(category==F)))
    print(table.out)
    table.out<-table.out[order(table.out$exp.freq,decreasing=T),]
    #names(table.out)[c(1,2)]<-c("Log(gc/ul)","Freq")
    #table.out$Freq<-paste0(as.character(table.out$Freq*100),"%")
    },digits=4)
  
  #output$cov<-renderPlot({})
  # output$samp.dis<-renderPlot({
  #   remain.df<-data.remain()
  #   ddply(remain.df,~mutation+category,summarize,count=length(mutation))->remain.count
  #   
  #   ggplot(remain.count,aes(x=count,fill=category))+geom_histogram(position="dodge")+ggtitle("distribution of variants across samples")+xlab("number of samples")+ylab("number of varinats")+scale_y_log10()
  # })
  output$freq<-renderPlot({
    remain.df<-data.remain()
    ggplot(remain.df,aes(x=freq.var,fill=category))+geom_histogram(position="dodge")+ggtitle("distribution of frequency")+xlab("Frequency")+ylab("number of varinats")+scale_x_log10()
  })
  output$mean.pos<-renderPlot({
    remain.df<-data.remain()
    ggplot(remain.df,aes(x=Read_pos,fill=category))+geom_histogram(position="dodge")+ggtitle("Average read position")+xlab("Read Position")+ylab("number of varinats")
  })
  output$mean.qual<-renderPlot({
    remain.df<-data.remain()
    ggplot(remain.df,aes(x=Phred,y=MapQ,color=category))+geom_point()+ggtitle("Quality Distributions")+xlab("Phred")+ylab("MapQ")
  })
  # output$mean.mapq<-renderPlot({
  #   remain.df<-data.remain()
  #   ggplot(remain.df,aes(x=MapQ,fill=category))+geom_histogram(position="dodge")+ggtitle("Average MapQ")+xlab("MapQ")+ylab("number of varinats")
  # })
  # output$position<-renderPlot({
  #   remain.df<-data.remain()
  #   prior.seg.length<-c()
  #   
  #   for(k in 1:length(regions.bed$chr)){ 
  #     prior.seg.length[k]<-sum(regions.bed$stop[1:k])  # the end positions of each segment relative to one sequence not including the trimming step
  #   }
  #   
  #   prior.seg.length<-c(0,prior.seg.length)
  #   #remain.df<-mutate(remain.df,chr=gsub("N_A","NA",chr))
  #   remain.df<-mutate(remain.df, # the concatenated position relative to the usual order
  #                     concat.pos=pos+prior.seg.length[match(chr,regions.bed$chr)])
  #   
  #   
  #   FP<-subset(remain.df,category==F)
  #   ggplot(FP,aes(x=concat.pos))+geom_histogram(binwidth=50)+ggtitle("FP locations")+geom_line(data=regions.l,aes(x=value,y=0,col=chr),size=2)
  #   
  #   
  # })
  
  
  
  
  ##### Progress graph #####
observeEvent(input$save,{
  remain.df<-data.remain()
  table.out<-ddply(remain.df,~gc+exp.freq,summarize,wt_mut=length(which(category=="wt")),FP=length(which(category==F)),TP=length(which(category==T)),FDR=(possible_vars-20-length(wt_mut)-FP)/(possible_vars-20-length(wt_mut)),TDR=TP/20)
  
  saved<-summarize(table.out,samp=length(gc),mean.fp=mean(FP),mean.tp=mean(TP),correction=paste(input$method),p.val=input$p.val,MapQ=input$MapQ,freq=input$freq.var,read.range=input$pos[2]-input$pos[1],dups=input$dups,distribution=input$disp)
  
  table.save[input$save,]<<-saved
  output$saved.table<-renderTable({

      table.save},digits=4) 
    
  })  
})

