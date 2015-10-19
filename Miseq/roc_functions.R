require(pROC)
require(plyr)
require(Biostrings)
require(reshape2)

find_freq<-function(Id){
  if(Id=="Covaris_5"){
    x=0.05} else if ( Id=="Covaris_25"){
      x=0.025} else if ( Id=="Covaris_125"){
        x=0.0125} else if ( Id=="Covaris_063"){
          x=0.0063} else if ( Id=="Covaris_016"){
            x=0.0016}
}
find_freq.v<-Vectorize(find_freq)
reference.fasta<-"./data/pHW2000_PR8-N_A.fa"
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
## These functions make an ROC by allowing pROC to calculate the positions based on only the variants that are provided, and then we adjust the senesitivity for with what we know.  For example if only 3 TP are found in one sample pROC will call that 100% Sensitivity but we will adjust it to 3/possible_tp and so forth.
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
  
  mutate(roc.df, exp.freq=find_freq.v(samp))-> roc.df
}  

true_snv<-read.csv("./data/PR8_WSN33.csv",comment.char = "#")

true_snv<-subset(true_snv,Ref!="-" & Allele.1!="-")
mutate(true_snv,mutant=paste0(Name,"_",Ref,Ref.Pos,Allele.1))->true_snv
tp<-unique(true_snv$mutant)
possible_tp<-length(tp)
adjust.coords<-function(roc.df,sum.df){ # adjsut the sensitivity and specificity called by pROC  
  samp<-roc.df$samp[1]
  samp.df<-subset(sum.df,Id==samp)
  sense.factor<-length(which(samp.df$category==T))/possible_tp
  TN.samp<-length(which(samp.df$category==F))
  mutate(roc.df,adj.sensitivity=sensitivity*sense.factor,FP=(TN.samp-TN.samp*specificity),adj.specificity=(possible_vars-possible_tp-FP)/((possible_vars-possible_tp)))
}

### Limit to coding regions ###
coding<-read.csv("./data/pr8.coding.bed.csv")

coding.cut<-function(x){
  chr<-unique(x$chr)
  start<-coding$start[match(x$chr,coding$chr)]
  stop<-coding$stop[match(x$chr,coding$chr)]
  
  subset(x,pos>start & pos<stop)
}
