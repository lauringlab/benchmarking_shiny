## Read in data

require(plyr)
require(reshape2)
require(Biostrings)
deep.sum<-list.files(path = "./data/",pattern = "sum.csv")
primer.cut<-function(x){ # a helper function to do this
  chr<-unique(x$chr)
  start<-regions.bed$start[match(x$chr,regions.bed$chr)]
  stop<-regions.bed$stop[match(x$chr,regions.bed$chr)]
  
  subset(x,pos>start & pos<stop)
}



for (i in 1:length(deep.sum)){
sum.data<-deep.sum[i]

sum.df<-read.csv(paste0("./data/",sum.data),stringsAsFactors=F)
sum.df$Id[sum.df$Id=="3_03"]<-"3_02" # correct an old naming error


### adapted from bencmarking_paper/scripts/figures.Rmd ######

mfe.df<-read.csv("./reference/mfe.csv") 
# Getting genome lenght and segment information from the wsn33 fasta file
reference.fasta<-"./reference/wsn33_wt_plasmid.fa" 
segments <- fasta.seqlengths(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=12, stop=segments-13, row.names=NULL) # the univeral primers are 12 and 13 bp long
regions.bed<-mutate(regions.bed,chr=as.character(chr))
# The samples are named in the following format log(genome copy/ul)_expected frequency of variants We'll use that information here to get the gc# and exp.freq for each variant called

pat<-"([0-9]+)_([0-9]+)" 
mutate(sum.df, gc = sub(pat,"\\1",Id),
       exp.freq = sub(pat,"\\2",Id)
)-> sum.df
mutate(sum.df,exp.freq=ifelse(test = grepl("0",exp.freq),yes=as.numeric(exp.freq)/1000,no = as.numeric(exp.freq)/100))->sum.df


sum.df<-ddply(sum.df,~chr,primer.cut) # removes any variants outside of priming sites (in this case it removes 0)
# Id the variants as T or F 
true_snv<-read.csv("./reference/mutant_id.csv",stringsAsFactor=F) # get the T variants
sum.df<-mutate(sum.df,category=mutation %in% true_snv$mutant) # add Column for True and false variants
wt1_mut<-subset(sum.df,Id=="WT1" & freq.var>0.01,select=mutation)
wt2_mut<-subset(sum.df,Id=="WT2" &freq.var>0.01 ,select=mutation)

wt_mut<-intersect(wt1_mut$mutation,wt2_mut$mutation) # there aren't any of these

wt.muts<-data.frame(wt=length(wt_mut)) #How many wt muts do we find.
sum.df$category[sum.df$mutation %in% wt_mut]<-"wt"
sum.df<-subset(sum.df,!(grepl("WT",Id))) # remove the WT sample

sum.df$category[sum.df$mutation %in% mfe.df$mutation]<-"mfe"


### read data

#read.df<-read.csv(paste0("./data/",read.data),stringsAsFactors=F)


#read.df<-mutate(read.df,category=sum.df$category[match(Mutation,sum.df$mutation)])

## write data

write.csv(sum.df,paste0("./processed_data/",sum.data))
write.csv(wt.muts,paste0("./processed_data/",sum.data,".wt.csv"))
#write.csv(read.df,paste0("./processed_data/",read.data))
}
