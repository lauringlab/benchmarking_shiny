## Read in data

require(plyr)
require(reshape2)
deep.sum<-list.files(path = "./data/",pattern = "sum.csv")

for (i in 1:length(deep.sum)){
sum.data<-deep.sum[i]

sum.df<-read.csv(paste0("./data/",sum.data),stringsAsFactors=F)
sum.df$Id[sum.df$Id=="3_03"]<-"3_02" # correct a naming error

WT<-subset(sum.df,grepl("WT",Id)) # Separate WT samples
sum.df<-subset(sum.df,!(grepl("WT",Id)))

pat<-"([0-9]+)_([0-9]+)"
mutate(sum.df, gc = sub(pat,"\\1",Id),
       exp.freq = sub(pat,"\\2",Id)
)-> sum.df
mutate(sum.df,exp.freq=ifelse(test = grepl("0",exp.freq),yes=as.numeric(exp.freq)/1000,no = as.numeric(exp.freq)/100))->sum.df

# TF column
true_snv<-read.csv("./data/mutant_id.csv",stringsAsFactor=F) # get the TP
sum.df<-mutate(sum.df,category=mutation %in% true_snv$mutant) # add Column for True and false variants
wt1_mut<-subset(WT,Id=="WT1",select=mutation)
wt2_mut<-subset(WT,Id=="WT2",select=mutation)

wt_mut<-intersect(wt1_mut$mutation,wt2_mut$mutation)

length(wt_mut)
sum.df$category[sum.df$mutation %in% wt_mut]<-"wt"

### read data

#read.df<-read.csv(paste0("./data/",read.data),stringsAsFactors=F)


#read.df<-mutate(read.df,category=sum.df$category[match(Mutation,sum.df$mutation)])

## write data

write.csv(sum.df,paste0("./processed_data/",sum.data))
#write.csv(read.df,paste0("./processed_data/",read.data))
}
