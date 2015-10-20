## Read in data
require(reshape2)
require(plyr)
lo.freqsum<-list.files(path = "./data/",pattern = "vcf_csv")
deep.sum<-list.files(path = "./data/",pattern = "sum.csv")
find_freq<-function(Id){
  if(Id=="Covaris_5"){
    x=0.05} else if ( Id=="Covaris_25"){
    x=0.025} else if ( Id=="Covaris_125"){
    x=0.0125} else if ( Id=="Covaris_063"){
    x=0.0063} else if ( Id=="Covaris_016"){
    x=0.0016}
}
find_freq.v<-Vectorize(find_freq)
## deepSNV ##

for (i in 1:length(deep.sum)){
  sum.data<-deep.sum[i]
  sum.df<-read.csv(paste0("./data/",sum.data),stringsAsFactors=F,comment.char = '#')
  
  true_snv<-read.csv("./data/PR8_WSN33.csv",comment.char = '#')
  true_snv<-subset(true_snv,Ref!="-" & Allele.1!="-")
  mutate(true_snv,mutant=paste0(Name,"_",Ref,Ref.Pos,Allele.1))->true_snv
  
  PR8_var<-read.csv("./data/PR8_variants.csv",comment.char = "#")
  
  mutate(PR8_var,mutant=paste0(Name,"_",Ref,Ref.Pos,Allele.1))->PR8_var
  sum.df<-mutate(sum.df,category=mutation %in% true_snv$mutant) # add Column for True and false variants
  sum.df$category[sum.df$mutation %in% PR8_var$mutant]<-"PR8"
  sum.df<-subset(sum.df,grepl("Covaris",Id))
  sum.df<-mutate(sum.df,exp.freq=find_freq.v(Id))
  
  write.csv(sum.df,paste0("./processed_data/",sum.data))
}


## lofreq ##

for(i in 1:length(lo.freqsum)){
  
  sum.data<-lo.freqsum[i]
  #read.data<-read_data_files[i]
  
  sum.df<-read.csv(paste0("./data/",sum.data),stringsAsFactors=F,comment.char = '#')
  
  sum.df<-rename(sum.df,c("CHROM"="chr",
                          "POS"="pos",
                          "ID"="Id",
                          "REF"="ref",
                          "ALT"="var",
                          "AF"="freq.var"
  ))
  
  sum.df<-mutate(sum.df,mutation=paste0(chr,"_",ref,pos,var))
  
  true_snv<-read.csv("./data/PR8_WSN33.csv",comment.char = '#')
  true_snv<-subset(true_snv,Ref!="-" & Allele.1!="-")
  mutate(true_snv,mutant=paste0(Name,"_",Ref,Ref.Pos,Allele.1))->true_snv
  sum.df<-mutate(sum.df,p.val=10^(QUAL/-10))
  PR8_var<-read.csv("./data/PR8_variants.csv",comment.char = "#")
  
  mutate(PR8_var,mutant=paste0(Name,"_",Ref,Ref.Pos,Allele.1))->PR8_var
  sum.df<-mutate(sum.df,category=mutation %in% true_snv$mutant) # add Column for True and false variants
  sum.df$category[sum.df$mutation %in% PR8_var$mutant]<-"PR8"
  sum.df<-subset(sum.df,grepl("Covaris",Id))
  sum.df<-mutate(sum.df,exp.freq=find_freq.v(Id))

  
  write.csv(sum.df,paste0("./processed_data/",sum.data))
  

}