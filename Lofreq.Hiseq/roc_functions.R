require(pROC)
require(plyr)
require(Biostrings)
require(reshape2)

reference.fasta<-"./wsn33_wt_plasmid.fa"
segments <- fasta.seqlengths(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=1, stop=segments, row.names=NULL)
regions.bed<-mutate(regions.bed,chr=as.character(chr))
#regions.bed<-mutate(regions.bed,chr=gsub("N_A","NA",chr))



possible_vars<-sum(regions.bed$stop)*3 # thre possible variants/ position
id<-function(x){
  x<-strsplit(x,".",fixed=T)[[1]][1]
}
prior.seg.length<-c()

for(k in 1:length(regions.bed$chr)){
  prior.seg.length[k]<-sum(regions.bed$stop[1:k])  # the end positions of each segment relative to one sequence not including the trimming step
}

prior.seg.length<-c(0,prior.seg.length)


mutate(regions.bed,concat.start=start+prior.seg.length[match(chr,regions.bed$chr)],concat.stop=concat.start+stop)->regions.bed

melt(subset(regions.bed,select=-c(stop,start)),id.vars="chr")->regions.l

seg<-function (x) {
  x <- strsplit(x, ".", fixed = T)[[1]][2]
}  
## These functions make an ROC by allowing pROC to calculate the positions based on only the variants that are provided, and then we adjust the senesitivity for with what we know.  For example if only 3 TP are found in one sample pROC will call that 100% Sensitivity but we will adjust it to 3/20 and so forth.
sum_roc<-function(x){
  sum.df<-subset(x,category %in% c(TRUE,FALSE),select=c(category,p.val))# get the TRUE and false variant calls and p.vals
  if(length(which(sum.df$category==T))>0){ # filter out cases where there aren't any TP found
    if(length(which(sum.df$category==F))==0){
      sum.df<-rbind(sum.df,data.frame(category=F,p.val=1))
    }  
    roc(sum.df$category~sum.df$p.val,plot=F,CI=T)
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
  
  pat<-"([0-9]+)_([0-9]+)"
  mutate(roc.df, gc_ul = sub(pattern = pat,"\\1",samp), # add columns to id the sample and the expected freq
         exp.freq = sub(pattern = pat,"\\2",samp),
  )-> roc.df
}  
adjust.coords<-function(roc.df,sum.df){ # adjsut the sensitivity and specificity called by pROC  
  samp<-roc.df$samp[1]
  samp.df<-subset(sum.df,Id==samp)
  sense.factor<-length(which(samp.df$category==T))/20
  TN.samp<-length(which(samp.df$category==F))
  mutate(roc.df,adj.sensitivity=sensitivity*sense.factor,FP=(TN.samp-TN.samp*specificity),adj.specificity=(possible_vars-20-FP)/((possible_vars-20)))
}

### Limit to coding regions ###
coding<-read.csv("./data/wsn33.coding.bed.csv")

coding.cut<-function(x){
  chr<-unique(x$chr)
  start<-coding$start[match(x$chr,coding$chr)]
  stop<-coding$stop[match(x$chr,coding$chr)]
  
  subset(x,pos>start & pos<stop)
}
