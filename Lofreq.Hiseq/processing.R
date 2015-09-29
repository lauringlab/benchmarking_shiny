## Read in data
require(reshape2)
require(plyr)
summaries<-c("all.removed.mapq_vcf.vcf_csv.csv")
#read_data_files<-c("BH.25.reads.csv",  "BH.30.reads.csv",   "bon.25.reads.csv", "bon.30.reads.csv")

for(i in length(summaries)){
  
sum.data<-summaries[i]
#read.data<-read_data_files[i]

sum.df<-read.csv(paste0("./data/",sum.data),stringsAsFactors=F,comment.char = '#')

sum.df<-rename(sum.df,c("CHROM"="chr",
                        "POS"="pos",
                        "ID"="Id",
                        "REF"="ref",
                        "ALT"="var",
                        "AF"="freq.var"
))

mutate(sum.df,mutation=paste0(chr,"_",ref,pos,var),category=mutation %in% tp$mutant)->sum.df

sum.df$Id[sum.df$Id=="3_03"]<-"3_02" # correct a naming error

WT<-subset(sum.df,grepl("WT",Id)) # Separate WT samples 
sum.df<-subset(sum.df,!(grepl("WT",Id))& Id!="Plasmid_control")



sum.df<-mutate(sum.df,p.val=10^(QUAL/-10))

pat<-"([0-9]+)_([0-9]+)"
mutate(sum.df, gc = sub(pat,"\\1",Id),
       exp.freq = sub(pat,"\\2",Id)
)-> sum.df
mutate(sum.df,exp.freq=ifelse(test = grepl("0",exp.freq),yes=as.numeric(exp.freq)/1000,no = as.numeric(exp.freq)/100))->sum.df

# TF column
true_snv<-read.csv("../Hiseq/data/mutant_id.csv",stringsAsFactor=F) # get the TP 
sum.df<-mutate(sum.df,category=mutation %in% true_snv$mutant) # add Column for True and false variants
wt1_mut<-subset(WT,Id=="WT1",select=mutation)
wt2_mut<-subset(WT,Id=="WT2",select=mutation)

wt_mut<-intersect(wt1_mut$mutation,wt2_mut$mutation)

length(wt_mut)
sum.df$category[sum.df$mutation %in% wt_mut]<-"wt"

### read data

# read.df<-read.csv(paste0("./data/",read.data),stringsAsFactors=F)
# 
# 
# read.df<-mutate(read.df,category=sum.df$category[match(Mutation,sum.df$mutation)])

## write data

write.csv(sum.df,paste0("./processed_data/",sum.data))
# write.csv(read.df,paste0("./processed_data/",read.data))
}