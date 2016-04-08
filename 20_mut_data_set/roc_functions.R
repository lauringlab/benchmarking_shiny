require(pROC)
require(plyr)
require(Biostrings)
require(reshape2)
pat<-"([0-9]+)_([0-9]+)"
primer.cut<-function(x){ # a helper function to do this
  chr<-unique(x$chr)
  start<-regions.bed$start[match(x$chr,regions.bed$chr)]
  stop<-regions.bed$stop[match(x$chr,regions.bed$chr)]
  
  subset(x,pos>start & pos<stop)
}

mfe.df<-read.csv("./reference/mfe.csv")
reference.fasta<-"./reference/wsn33_wt_plasmid.fa" 
segments <- fasta.seqlengths(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=12, stop=segments-13, row.names=NULL) # the univeral primers are 12 and 13 bp long
regions.bed<-mutate(regions.bed,chr=as.character(chr))
regions.bed<-mutate(regions.bed,length=stop-start-1)
# The samples are named in the following format log(genome copy/ul)_expected frequency of variants We'll use that information here to get the gc# and exp.freq for each variant called
total_possible_vars<-sum(regions.bed$length)*3
poss_vars<-function(mfe.df,wt_mut){
  total_possible_vars<-sum(regions.bed$length)*3 # three possible variants/ position
  possible_vars<-total_possible_vars-wt_mut-dim(mfe.df)[1]-20
  }


find_freq.h<-function(Id){ # helper function to get the expected frequency from the sample name
  if(Id==0.05){
    x="5.0%"} else if ( Id==0.02){
      x="2.0%"} else if ( Id==0.01){
        x="1.0%"} else if ( Id==0.005){
          x="0.5%"} else if ( Id==0.002){
            x="0.2%"}
}

find_freq.v.h<-Vectorize(find_freq.h) # vectorized version of the helper function above





id<-function(x){
  x<-strsplit(x,".",fixed=T)[[1]][1]
}

seg<-function (x) {
  x <- strsplit(x, ".", fixed = T)[[1]][2]
}  

fill_in_plots<-function(x){
  spec<-min(x$adj.specificity[which(x$adj.specificity>0.995)])
  sense<-max(x$adj.sensitivity[x$adj.specificity==spec]) # incase there is a step right here
  
  y<-data.frame(threshold=x$threshold[1],Id=x$Id[1],adj.specificity=0.9950001,adj.sensitivity=sense,samp=unique(x$samp),exp.freq=x$exp.freq[1],FP=x$FP[1],TP=x$TP[1],sensitivity=x$sensitivity[1],specificity=x$specificity[1])
  
  if("gc" %in% names(x)){
    y$gc=x$gc[1]
  }
  return(y)
}
sum_roc<-function(x,direction){
  sum.df<-subset(x,category %in% c(TRUE,FALSE),select=c(category,p.val))# get the TRUE and false variant calls and p.vals
  if(length(which(sum.df$category==T))>0){ # filter out cases where there aren't any TP found
    if(length(which(sum.df$category==F))==0){
      sum.df<-rbind(sum.df,data.frame(category=F,p.val=1))
    }  
    roc(sum.df$category~sum.df$p.val,plot=F,CI=T,direction=direction)
  }
}

roc_df<-function(roc_analysis){ # get the coordinates and cut offs for all the points in the ROC object
  roc_analysis<-roc_analysis[unlist(lapply(roc_analysis,is.null))==F]
  all<-lapply(roc_analysis,coords,x="all")
  all.long<-lapply(all, melt) 
  all.long<-lapply(all.long, function(x){ 
    mutate(x, 
           Id=c(0,head(as.numeric(rownames(x))%/%3,-1))
    )
  }) # Set an id column to group the threshold,specificity, and sensitivity together with the same number by integer division. I adjust with c(0, head ..., -1) since 3%/%3=1 but it should be grouped with the 0's.
  roc.ls<-lapply(all.long,function(x) dcast(x,Id~Var1)) # dcast into columns of threshold, specificity and sensitivity
  roc.df<-do.call(rbind,roc.ls) # combine into data frame
  roc.df$samp<-unlist(lapply(rownames(roc.df),id)) # add sample id column
  
  mutate(roc.df, exp.freq=find_freq.v.h(samp))-> roc.df 
}  

adjust.coords<-function(roc.df,sum.df,possible_tp,possible_vars){ # adjust the sensitivity and specificity called by pROC  
  samp<-roc.df$samp[1]
  samp.df<-subset(sum.df,Id==samp)
  sense.factor<-length(which(samp.df$category==T))/possible_tp
  TN.samp<-length(which(samp.df$category==F))
  mutate(roc.df,adj.sensitivity=sensitivity*sense.factor,FP=(TN.samp-TN.samp*specificity),TP=adj.sensitivity*possible_tp,adj.specificity=(possible_vars-possible_tp-FP)/((possible_vars-possible_tp)))
}

roc_df.one<-function(roc_anal,thr){ # get the coordinates and cut offs for all the points in the ROC object
  roc_analysis<-roc_anal[unlist(lapply(roc_anal,is.null))==F]
  all<-lapply(roc_analysis,coords,x=thr,input="thr")
  
  #roc.ls<-lapply(all.long,function(x) dcast(x,Id~Var1)) # dcast into columns of threshold, specificity and sensitivity
  roc.df<-as.data.frame(do.call(rbind,all)) # combine into data frame
  roc.df$samp<-rownames(roc.df) # add sample id column
  
  mutate(roc.df, exp.freq=find_freq.v.h(samp))-> roc.df
}

hiseq.roc<-function(sum.df,possible_tp,possible_vars,direction){ # A function to run the roc calculations and adjustments
  roc.ls<-dlply(sum.df,~Id,sum_roc,direction)
  
  roc.df<-roc_df(roc.ls)
  
  
  roc.df.adj<-ddply(roc.df,~samp,adjust.coords,sum.df,possible_tp,possible_vars)
  roc.df.adj<-roc.df.adj[order(roc.df.adj$adj.sensitivity),]
  mutate(roc.df.adj, gc = as.numeric(sub(pat,"\\1",samp)), # This fills in the null column given in the roc_df function above. line 161
         exp.freq = sub(pat,"\\2",samp)
  )-> roc.df.adj
  fills<-ddply(roc.df.adj,~samp,fill_in_plots) # take out these lines to not fill in the plots
  roc.df.adj<-rbind(fills,roc.df.adj) # here too
  roc.df.adj$exp.freq<-as.character(roc.df.adj$exp.freq) # so that teh calculations below work since it is a factor here
  mutate(roc.df.adj,exp.freq=ifelse(test = grepl("0",exp.freq),yes=as.numeric(exp.freq)/1000,no = as.numeric(exp.freq)/100))->roc.df.adj
  roc.df.adj$exp.freq<-as.factor(roc.df.adj$exp.freq)
  roc.df.adj$exp.freq <- factor(roc.df.adj$exp.freq, levels = rev(levels(roc.df.adj$exp.freq)))
  return(roc.df.adj)
}

hiseq.roc.table<-function(sum.df,cut.off,possible_tp,possible_vars,direction){
  roc.ls<-dlply(sum.df,~Id,sum_roc,direction)
  roc.df<-roc_df.one(roc.ls,cut.off)
  roc.df.adj<-ddply(roc.df,~samp,adjust.coords,sum.df,possible_tp,possible_vars)
  mutate(roc.df.adj, gc = as.numeric(sub(pat,"\\1",samp)),
         exp.freq = sub(pat,"\\2",samp)
  )-> roc.df.adj
  
  mutate(roc.df.adj,exp.freq=ifelse(test = grepl("0",exp.freq),yes=as.numeric(exp.freq)/1000,no = as.numeric(exp.freq)/100))->roc.df.adj
  roc.df.adj<-roc.df.adj[order(roc.df.adj$exp.freq,decreasing = T),]
  roc.table<-subset(roc.df.adj,select=c(gc,exp.freq,adj.sensitivity,TP,adj.specificity,FP))
  return(roc.table)
}


### Limit to coding regions ###
coding<-read.csv("./data/wsn33.coding.bed.csv")

coding.cut<-function(x){
  chr<-unique(x$chr)
  start<-coding$start[match(x$chr,coding$chr)]
  stop<-coding$stop[match(x$chr,coding$chr)]
  
  subset(x,pos>start & pos<stop)
}
