#Project title: Disturbance reduces recruitment 
#Authors: N. M. Kollars and J. J. Stachowicz 
#Script by: Nicole Kollars

#Nicole's working directory
setwd("C:/Users/nmkol/OneDrive/Desktop/Feedback_MS")

#########
#needed libraries and functions 
#########

library(RClone)
library(related)
library(reshape)
library(dplyr)

##########
#read in the data
##########

#genetic data
dat<-read.csv("Kollars_feedback_usats.csv")
dat<-dat[,1:29]###removes locus 11
dat<-dat[,-c(24,25)]###removes locus 8

#this gives the number of positions that were not sampled in the field 
dat.nocollection<-na.omit(dat[rowSums(dat[,10:27]) == 0, ])  
#this gives the samples that are missing data due to genotyping error
dat.NA<-dat[rowSums(is.na(dat)) > 0, ]  

##########
#Determining multi-locus genotypes
##########

#create a dataframe where individuals with missing alleles are removed
dat<-na.omit(dat)
#and remove samples that were not collected from the field
dat<-anti_join(dat,dat.nocollection)

#creates dataframes needed for RClone
#vector of populations
popvec<-dat[,1]
#vector of identifying info
idvec<-dat[,2:9]
#microsat data for RClone functions  
dat2<-dat[,10:ncol(dat)]

#psex, calculated using round robin method
p<-psex(dat2,RR=TRUE)
mlg<-as.numeric(p$genet)
mlg[c(which(is.na(mlg)))]<-c(which(is.na(mlg)))
#create a dataframe that binds the mlg designation to the dat
mlg<-cbind(dat,mlg)

#rarefaction curve of loci  
res<-sample_loci(dat2,nbrepeat=1000)
res.dat<-as.data.frame(do.call(cbind,res$raw_MLG))
boxplot(res$raw_MLG,
        xlab = "Number of loci sampled", ylab="Number of MLGs",
        ylim = c(0,140))

##########
#Using relatedness to group MLGs into MLLs
##########

#taking out repeated genotypes of assigned samples 
norep_mlg<-distinct(mlg,mlg,.keep_all=TRUE)
norep_mlg$id<-paste0(norep_mlg$block,norep_mlg$plot,"_",norep_mlg$sample,"T",norep_mlg$time)
norep<-norep_mlg[,10:29]
#create a file of the no repeats that is in the correct format 
#for the related package (no headers, tab-deliminated text file)
norep<- norep[,c(20,1:18)]
write.table(norep, file="norep.txt", sep="\t", 
            col.names=FALSE, row.names=FALSE, quote=FALSE)

#relatedness analysis
input<-readgenotypedata("norep.txt")
outfile<-coancestry("norep.txt",lynchli=1)
related<-outfile$relatedness

hist(related$lynchli,breaks=150)

#modify "related" file so that it is merged with the rest of the
#datset 
reldat<-related%>%
  rename(c("one" = "ind1.id","two"="ind2.id")) %>%
  select(pair.no,one,two,lynchli)
reldat$pair.no<-as.factor(reldat$pair.no)
reldat<-melt(reldat,id=c("pair.no","lynchli"))
reldat<-arrange(reldat,pair.no)
names(reldat)[3:4] <- c("indiv","id")

reldatall<-left_join(reldat,norep_mlg, by="id")
reldatall<-reldatall[,-c(5:13)] #this gets rid of the "X.1" column
reldatall<-arrange(reldatall,desc(lynchli))

#figure out number of repeats for each mlg 
num.repeats=data.frame()
for(i in unique(as.factor(mlg$mlg))){
  test<-mlg[mlg$mlg==i,]
  N<-length(test$mlg)
  df<-data.frame(cbind(i,N))
  num.repeats<-rbind(num.repeats,df)
}
names(num.repeats)[1] <- "mlg"

#combine dataset with relatedness, mlg, # repeats per mlg
reldatall$mlg<-as.character(reldatall$mlg)
relatedness<-left_join(reldatall,num.repeats,by="mlg")
relatedness <-relatedness[, c(1:4,23:24,5:22)]

#No. of allele differences between pairs for each pair combo
no.alleles=data.frame()
for(i in unique(as.factor(relatedness$pair.no))){
  test<-relatedness[relatedness$pair.no==i,]
  z<-sum((test[1,7:24]-test[2,7:24])!=0)
  df<-data.frame(cbind(i,z))
  no.alleles<-rbind(no.alleles,df)
}
names(no.alleles)[1] <- "pair.no"

#creates a dataframe with rxy for each pair # 
pairs<-select(relatedness,c("pair.no","lynchli"))
pairs<-distinct(pairs,pair.no, .keep_all= TRUE)
pairs.no<-left_join(pairs,no.alleles, by="pair.no")

####grouping of mlls####
#a priori rule: screen pairs in same/adjacent plots with 1 allele difference
oneallele<-pairs[pairs.no$z=="1",]
pairs.mll.test<-left_join(oneallele,relatedness,by="pair.no")

##removing pair 8532 because the two samples (which are both colonizers) 
#are spaitally distanced 
pairs.mll.test<-pairs.mll.test[pairs.mll.test$pair.no!="8532",]

#change alleles at distinct locus for each pair to NA 
mll.dat=data.frame()
for(i in unique(as.factor(pairs.mll.test$pair.no))){
  test<-pairs.mll.test[pairs.mll.test$pair.no==i,]
  diff<-test[1,8:25]-test[2,8:25]
  is.even <- function(x) x%%2 == 0
  diff[diff!=0]<-diff[ifelse(is.even(which(diff!=0)),which(diff!=0)-1,which(diff!=0)+1)]<-NA
  test[,which(is.na(diff))+7]<-NA
  mll.dat<-rbind(mll.dat,test)
}

#within each pair, replicate the number of rows for that pair to match
#number in the sample, remove the distinct locus and run psex 

mll.psex<-data.frame()
for( i in unique(mll.dat$pair.no)){
  pair1<-mll.dat[mll.dat$pair.no==i,]
  pair1<-as.data.frame(lapply(pair1, rep, pair1$N))
  idvec<-pair1[,1:7]
  pair.dat<-pair1[,8:25]
  pair.dat<-pair.dat %>% select_if(~ !any(is.na(.)))
  p<-psex(pair.dat,RR=TRUE)
  df<-data.frame(cbind(i,p))
  mll.psex<-rbind(mll.psex,df)
}

#psex confirms that these should be the same clone

#this provides the list of mlgs that need to be changed w/modified data  
to.change=data.frame()
for(i in unique(as.factor(pairs.mll.test$pair.no))){
  test<-pairs.mll.test[pairs.mll.test$pair.no==i,]
  high.N<-ifelse(as.numeric(test[test$indiv=="one",]$N)>
                   as.numeric(test[test$indiv=="two",]$N),"one","two")
  match.to<-test[test$indiv==high.N,]
  to.match<-test[test$indiv!=high.N,]
  to.match[,c(8:25)]<-match.to[,c(8:25)]
  to.match$real.mlg<-match.to$mlg
  to.change<-rbind(to.change,to.match)
}

to.change<-to.change[,c(6,26,8:25)]
to.change$mlg<-as.numeric(to.change$mlg)

#find all the need to change mlgs in the mlg dataset
samples2change<-mlg[mlg$mlg %in% to.change$mlg,]

#and swamp values in samples2 change with values in to.change 
test<-samples2change %>%
  left_join(to.change, by = "mlg")%>%
  mutate(mlg = real.mlg) %>%
  mutate(Zm01_1 = Zm01_1.y) %>%
  mutate(Zm01_2 = Zm01_2.y) %>%
  mutate(Zm02_1 = Zm02_1.y) %>%
  mutate(Zm02_2 = Zm02_2.y) %>% 
  mutate(Zm03_1 = Zm03_1.y) %>%
  mutate(Zm03_2 = Zm03_2.y) %>%
  mutate(Zm04_1 = Zm04_1.y) %>%
  mutate(Zm04_2 = Zm04_2.y) %>%
  mutate(Zm05_1 = Zm05_1.y) %>%
  mutate(Zm05_2 = Zm05_2.y) %>%
  mutate(Zm06_1 = Zm06_1.y) %>%
  mutate(Zm06_2 = Zm06_2.y) %>%
  mutate(Zm07_1 = Zm07_1.y) %>%
  mutate(Zm07_2 = Zm07_2.y) %>%
  mutate(Zm09_1 = Zm09_1.y) %>%
  mutate(Zm09_2 = Zm09_2.y) %>%
  mutate(Zm10_1 = Zm10_1.y) %>%
  mutate(Zm10_2 = Zm10_2.y)
test2<-test[,c(1:9,48:65,28)]
test2$mlg<-as.numeric(test2$mlg)

mlg.modified<-anti_join(mlg,samples2change)
mlg2<-full_join(mlg.modified,test2)

# Sort by vector name [z] then [x]
mlg2<-mlg2[
  with(mlg2, order(time, block,plot,sample)),
  ]

####################
###imputing########
##################
#confidentally (based on accumulation curve) need 7 loci
#selecting samples only missing data at 1 or 2 loci
dat.NA$sum<-rowSums(is.na(dat.NA[,10:27]))
to.impute<-dat.NA[dat.NA$sum<=4,]
to.impute<-to.impute[,-28]

#list of samples already assigned (without repeats)
#taking out repeated genotypes of assigned samples 
norep_mlg2<-distinct(mlg2,mlg,.keep_all=TRUE)

#merge to be assigned with already assigned 
impute.test<-full_join(norep_mlg2,to.impute)

#add id and prep for related package 
impute.test$id<-paste0(impute.test$block,impute.test$plot,"_",
                       impute.test$sample,"T",impute.test$time)
imputing<-impute.test[,10:29]

#create a file of the no repeats that is in the correct format 
#for the related package (no headers, tab-deliminated text file)
imputing<- imputing[,c(20,1:18)]
write.table(imputing, file="imputing.txt", sep="\t", 
            col.names=FALSE, row.names=FALSE, quote=FALSE)

#relatedness analysis
input<-readgenotypedata("imputing.txt")
outfile<-coancestry("imputing.txt",lynchli=1)
related<-outfile$relatedness

hist(related$lynchli,breaks=150)

#modify "related" file so that it is merged with the rest of the
#datset 
reldat<-related%>%
  rename(c("one" = "ind1.id","two"="ind2.id")) %>%
  select(pair.no,one,two,lynchli)
reldat$pair.no<-as.factor(reldat$pair.no)
reldat<-melt(reldat,id=c("pair.no","lynchli"))
reldat<-arrange(reldat,pair.no)
names(reldat)[3:4] <- c("indiv","id")

reldatall<-left_join(reldat,impute.test, by="id")
reldatall<-reldatall[,-c(5:13)] #this gets rid of the "X.1" column
reldatall<-arrange(reldatall,desc(lynchli))

#figure out number of repeats for each mlg 
num.repeats=data.frame()
for(i in unique(as.factor(mlg2$mlg))){
  test<-mlg2[mlg2$mlg==i,]
  N<-length(test$mlg)
  df<-data.frame(cbind(i,N))
  num.repeats<-rbind(num.repeats,df)
}
names(num.repeats)[1] <- "mlg"

#combine dataset with relatedness, mlg, # repeats per mlg
reldatall$mlg<-as.character(reldatall$mlg)
relatedness<-left_join(reldatall,num.repeats,by="mlg")
relatedness <-relatedness[, c(1:4,23:24,5:22)]

check.impute<-relatedness[relatedness$lynchli=="1",]

#this selects for all of the samples that are matching with a clone
check.impute<-check.impute[1:198,]
#imputing
dat.impute=data.frame()
for(i in unique(as.factor(check.impute$pair.no))){
  test<-check.impute[check.impute$pair.no==i,]
  test[2,5:24]<-test[1,5:24]
  dat.impute<-rbind(dat.impute,test)
}
imputed.dat<-dat.impute[dat.impute$indiv=="two",]
imputed.dat<-imputed.dat[,c(4:5,7:24)]  

##these are the samples that were assigned to more than one mlg
duplicates<-imputed.dat[duplicated(imputed.dat$id),]
dup.list<-imputed.dat[imputed.dat$id %in% duplicates$id,]

#need to remove dup.list from imputed.dat
imput.nodup<-anti_join(imputed.dat,dup.list)

##these are the samples without "1" matches 
to.impute$id<-paste0(to.impute$block,to.impute$plot,"_",to.impute$sample,"T",to.impute$time)
no.match<-anti_join(to.impute,imputed.dat,by="id")

test<-impute.test %>%
  left_join(imput.nodup, by = c("id")) %>%
  mutate(Zm01_1 = ifelse(is.na(Zm01_1.x), Zm01_1.y,Zm01_1.x )) %>%
  mutate(Zm01_2 = ifelse(is.na(Zm01_2.x), Zm01_2.y,Zm01_2.x )) %>%
  mutate(Zm02_1 = ifelse(is.na(Zm02_1.x), Zm02_1.y,Zm02_1.x )) %>%
  mutate(Zm02_2 = ifelse(is.na(Zm02_2.x), Zm02_2.y,Zm02_2.x )) %>% 
  mutate(Zm03_1 = ifelse(is.na(Zm03_1.x), Zm03_1.y,Zm03_1.x )) %>%
  mutate(Zm03_2 = ifelse(is.na(Zm03_2.x), Zm03_2.y,Zm03_2.x )) %>%
  mutate(Zm04_1 = ifelse(is.na(Zm04_1.x), Zm04_1.y,Zm04_1.x )) %>%
  mutate(Zm04_2 = ifelse(is.na(Zm04_2.x), Zm04_2.y,Zm04_2.x )) %>%
  mutate(Zm05_1 = ifelse(is.na(Zm05_1.x), Zm05_1.y,Zm05_1.x )) %>%
  mutate(Zm05_2 = ifelse(is.na(Zm05_2.x), Zm05_2.y,Zm05_2.x )) %>%
  mutate(Zm06_1 = ifelse(is.na(Zm06_1.x), Zm06_1.y,Zm06_1.x )) %>%
  mutate(Zm06_2 = ifelse(is.na(Zm06_2.x), Zm06_2.y,Zm06_2.x )) %>%
  mutate(Zm07_1 = ifelse(is.na(Zm07_1.x), Zm07_1.y,Zm07_1.x )) %>%
  mutate(Zm07_2 = ifelse(is.na(Zm07_2.x), Zm07_2.y,Zm07_2.x )) %>%
  mutate(Zm09_1 = ifelse(is.na(Zm09_1.x), Zm09_1.y,Zm09_1.x )) %>%
  mutate(Zm09_2 = ifelse(is.na(Zm09_2.x), Zm09_2.y,Zm09_2.x )) %>%
  mutate(Zm10_1 = ifelse(is.na(Zm10_1.x), Zm10_1.y,Zm10_1.x )) %>%
  mutate(Zm10_2 = ifelse(is.na(Zm10_2.x), Zm10_2.y,Zm10_2.x )) %>%
  mutate(mlg = ifelse(is.na(mlg.x), mlg.y,mlg.x )) 
test2<-test[,c(1:9,49:67)]
imputed<-na.omit(test2)
imputed$mlg<-as.numeric(imputed$mlg)

#samples to add to mlg2 
to.add<-anti_join(imputed,impute.test)

mlg3<-rbind(mlg2,to.add)
# Sort by vector name [z] then [x]
mlg3<-mlg3[
  with(mlg3, order(time, block,plot,sample)),
  ]

##########################################
###Continuing to impute duplicated samples 
##########################################

#make new list of samples not in MLG 
dat.NA2<-anti_join(dat.NA,to.add,by="name")

#selecting samples only missing data at 1 or 2loci
dat.NA2$sum<-rowSums(is.na(dat.NA2[,10:27]))

#list of samples to be imputed 
to.impute<-dat.NA2[dat.NA2$sum<=4,]
to.impute<-to.impute[,-28]

#list of samples already assigned (without repeats)
#taking out repeated genotypes of assigned samples 
norep_mlg3<-distinct(mlg3,mlg,.keep_all=TRUE)

#merge to be assigned with already assigned 
impute.test<-full_join(norep_mlg3,to.impute)

#add id and prep for related package 
impute.test$id<-paste0(impute.test$block,impute.test$plot,"_",
                       impute.test$sample,"T",impute.test$time)
imputing<-impute.test[,10:29]

#create a file of the no repeats that is in the correct format 
#for the related package (no headers, tab-deliminated text file)
imputing<- imputing[,c(20,1:18)]
write.table(imputing, file="imputing.txt", sep="\t", 
            col.names=FALSE, row.names=FALSE, quote=FALSE)

#relatedness analysis
input<-readgenotypedata("imputing.txt")
outfile<-coancestry("imputing.txt",lynchli=1)
related<-outfile$relatedness

hist(related$lynchli,breaks=150)

#modify "related" file so that it is merged with the rest of the
#datset 
detach("package:plyr")
reldat<-related%>%
  rename(c("one" = "ind1.id","two"="ind2.id")) %>%
  dplyr::select(pair.no,one,two,lynchli)
reldat$pair.no<-as.factor(reldat$pair.no)
reldat<-melt(reldat,id=c("pair.no","lynchli"))
reldat<-arrange(reldat,pair.no)
names(reldat)[3:4] <- c("indiv","id")

reldatall<-left_join(reldat,impute.test, by="id")
reldatall<-reldatall[,-c(5:13)] #this gets rid of the "X.1" column
reldatall<-arrange(reldatall,desc(lynchli))

#figure out number of repeats for each mlg 
num.repeats=data.frame()
for(i in unique(as.factor(mlg2$mlg))){
  test<-mlg2[mlg2$mlg==i,]
  N<-length(test$mlg)
  df<-data.frame(cbind(i,N))
  num.repeats<-rbind(num.repeats,df)
}
names(num.repeats)[1] <- "mlg"

#combine dataset with relatedness, mlg, # repeats per mlg
reldatall$mlg<-as.character(reldatall$mlg)
relatedness<-left_join(reldatall,num.repeats,by="mlg")
relatedness <-relatedness[, c(1:4,23:24,5:22)]

check.impute<-relatedness[relatedness$lynchli=="1",]

#this selects for all of the samples that are matching with a clone
check.impute<-check.impute[1:12,]
#sample to impute 4B_1T5 is assigned to 307 and 485
#sample to impute 4D_23T5 is assigned to 307 and 485
#sample to impute 10C_18T5 is assigned to 484 and 1070 

#imputing
dat.impute=data.frame()
for(i in unique(as.factor(check.impute$pair.no))){
  test<-check.impute[check.impute$pair.no==i,]
  test[2,5:24]<-test[1,5:24]
  dat.impute<-rbind(dat.impute,test)
}
imputed.dat<-dat.impute[dat.impute$indiv=="two",]
imputed.dat<-imputed.dat[,c(4:5,7:24)]  
##based on location, 4B_1T5 should be 307
#4D_23T5 should be 307
#10C_18T5 shoud be 1070
#remove 484 and 485
imputed.dat$mlg<-as.numeric(imputed.dat$mlg)
imputed.dat<-imputed.dat[imputed.dat$mlg!=484,]
imputed.dat<-imputed.dat[imputed.dat$mlg!=485,]

test<-impute.test %>%
  left_join(imputed.dat, by = c("id")) %>%
  mutate(Zm01_1 = ifelse(is.na(Zm01_1.x), Zm01_1.y,Zm01_1.x )) %>%
  mutate(Zm01_2 = ifelse(is.na(Zm01_2.x), Zm01_2.y,Zm01_2.x )) %>%
  mutate(Zm02_1 = ifelse(is.na(Zm02_1.x), Zm02_1.y,Zm02_1.x )) %>%
  mutate(Zm02_2 = ifelse(is.na(Zm02_2.x), Zm02_2.y,Zm02_2.x )) %>% 
  mutate(Zm03_1 = ifelse(is.na(Zm03_1.x), Zm03_1.y,Zm03_1.x )) %>%
  mutate(Zm03_2 = ifelse(is.na(Zm03_2.x), Zm03_2.y,Zm03_2.x )) %>%
  mutate(Zm04_1 = ifelse(is.na(Zm04_1.x), Zm04_1.y,Zm04_1.x )) %>%
  mutate(Zm04_2 = ifelse(is.na(Zm04_2.x), Zm04_2.y,Zm04_2.x )) %>%
  mutate(Zm05_1 = ifelse(is.na(Zm05_1.x), Zm05_1.y,Zm05_1.x )) %>%
  mutate(Zm05_2 = ifelse(is.na(Zm05_2.x), Zm05_2.y,Zm05_2.x )) %>%
  mutate(Zm06_1 = ifelse(is.na(Zm06_1.x), Zm06_1.y,Zm06_1.x )) %>%
  mutate(Zm06_2 = ifelse(is.na(Zm06_2.x), Zm06_2.y,Zm06_2.x )) %>%
  mutate(Zm07_1 = ifelse(is.na(Zm07_1.x), Zm07_1.y,Zm07_1.x )) %>%
  mutate(Zm07_2 = ifelse(is.na(Zm07_2.x), Zm07_2.y,Zm07_2.x )) %>%
  mutate(Zm09_1 = ifelse(is.na(Zm09_1.x), Zm09_1.y,Zm09_1.x )) %>%
  mutate(Zm09_2 = ifelse(is.na(Zm09_2.x), Zm09_2.y,Zm09_2.x )) %>%
  mutate(Zm10_1 = ifelse(is.na(Zm10_1.x), Zm10_1.y,Zm10_1.x )) %>%
  mutate(Zm10_2 = ifelse(is.na(Zm10_2.x), Zm10_2.y,Zm10_2.x )) %>%
  mutate(mlg = ifelse(is.na(mlg.x), mlg.y,mlg.x )) 
test2<-test[,c(1:9,49:67)]
test2$mlg<-as.numeric(test2$mlg)

#samples to add to mlg2 
to.add<-anti_join(test2,impute.test)

mlg4<-rbind(mlg3,to.add)
# Sort by vector name [z] then [x]
mlg4<-mlg4[
  with(mlg4, order(time, block,plot,sample)),
  ]

##########################################
###Continuing to impute: unmatched samples 
##########################################

#make new list of samples not in MLG 
dat.NA3<-anti_join(dat.NA2,to.add,by="name")

#selecting samples only missing data at 1 or 2 loci
dat.NA3$sum<-rowSums(is.na(dat.NA3[,10:27]))

#list of samples to be imputed 
to.impute<-dat.NA3[dat.NA3$sum<=4,]
to.impute<-to.impute[,-28]

#list of samples already assigned (without repeats)
#taking out repeated genotypes of assigned samples 
norep_mlg4<-distinct(mlg4,mlg,.keep_all=TRUE)

#merge to be assigned with already assigned 
impute.test<-full_join(norep_mlg4,to.impute)

#add id and prep for related package 
impute.test$id<-paste0(impute.test$block,impute.test$plot,"_",
                       impute.test$sample,"T",impute.test$time)
imputing<-impute.test[,10:29]

#create a file of the no repeats that is in the correct format 
#for the related package (no headers, tab-deliminated text file)
imputing<- imputing[,c(20,1:18)]
write.table(imputing, file="imputing.txt", sep="\t", 
            col.names=FALSE, row.names=FALSE, quote=FALSE)

#relatedness analysis
input<-readgenotypedata("imputing.txt")
outfile<-coancestry("imputing.txt",lynchli=1)
related<-outfile$relatedness

hist(related$lynchli,breaks=150)

#modify "related" file so that it is merged with the rest of the
#datset 
detach("package:plyr")
reldat<-related%>%
  rename(c("one" = "ind1.id","two"="ind2.id")) %>%
  dplyr::select(pair.no,one,two,lynchli)
reldat$pair.no<-as.factor(reldat$pair.no)
reldat<-melt(reldat,id=c("pair.no","lynchli"))
reldat<-arrange(reldat,pair.no)
names(reldat)[3:4] <- c("indiv","id")

reldatall<-left_join(reldat,impute.test, by="id")
reldatall<-reldatall[,-c(5:13)] #this gets rid of the "X.1" column
reldatall<-arrange(reldatall,desc(lynchli))

#figure out number of repeats for each mlg 
num.repeats=data.frame()
for(i in unique(as.factor(mlg2$mlg))){
  test<-mlg2[mlg2$mlg==i,]
  N<-length(test$mlg)
  df<-data.frame(cbind(i,N))
  num.repeats<-rbind(num.repeats,df)
}
names(num.repeats)[1] <- "mlg"

#combine dataset with relatedness, mlg, # repeats per mlg
reldatall$mlg<-as.character(reldatall$mlg)
relatedness<-left_join(reldatall,num.repeats,by="mlg")
relatedness <-relatedness[, c(1:4,23:24,5:22)]

#this sample is unique in that it is missing data at 1 locus and 
#different at 3/8 remaining loci; closest relative is far away  
T18A12<-impute.test[impute.test$id=="8A_12T1",]
T18A12$mlg<-2400
T18A12[,c(16:17)]<-999

#these samples should all be the same clone, but unique as they are 
#missing data at 1 loci but differ at 2 of the 8 remaining loci.
impute.test$id<-as.factor(impute.test$id)
dat7A<-impute.test%>%
  filter(id %in% c("7A_19T1","7A_26T1","7A_8T5","7A_21T5"))
dat7A$mlg<-2500
dat7A[,c(10:11)]<-999

#this sample is missing data at 2 loci and mismatches 1/7 loci, 
#with a nearby clone 
pairid<-arrange(relatedness[relatedness$id=="5D_34T1",],desc(lynchli))
test<-relatedness[relatedness$pair.no==5256,]
test[2,5:24]<-test[1,5:24]
T15D34<-test[test$indiv=="two",]

#this sample is missing data at 1 loci and mismatches 1/8 loci, 
#with a nearby clone 
pairid<-arrange(relatedness[relatedness$id=="9A_34T5",],desc(lynchli))
test<-relatedness[relatedness$pair.no==7289,]
test[2,5:24]<-test[1,5:24]
T59A34<-test[test$indiv=="two",]

imputed.dat3<-rbind(T15D34,T59A34)
imputed.dat3<-imputed.dat3[,c(4:5,7:24)]

imputed.dat4<-rbind(T18A12,dat7A) 
imputed.dat4<-imputed.dat4[,-29]

test<-impute.test %>%
  left_join(imputed.dat3, by = c("id")) %>%
  mutate(Zm01_1 = ifelse(is.na(Zm01_1.x), Zm01_1.y,Zm01_1.x )) %>%
  mutate(Zm01_2 = ifelse(is.na(Zm01_2.x), Zm01_2.y,Zm01_2.x )) %>%
  mutate(Zm02_1 = ifelse(is.na(Zm02_1.x), Zm02_1.y,Zm02_1.x )) %>%
  mutate(Zm02_2 = ifelse(is.na(Zm02_2.x), Zm02_2.y,Zm02_2.x )) %>% 
  mutate(Zm03_1 = ifelse(is.na(Zm03_1.x), Zm03_1.y,Zm03_1.x )) %>%
  mutate(Zm03_2 = ifelse(is.na(Zm03_2.x), Zm03_2.y,Zm03_2.x )) %>%
  mutate(Zm04_1 = ifelse(is.na(Zm04_1.x), Zm04_1.y,Zm04_1.x )) %>%
  mutate(Zm04_2 = ifelse(is.na(Zm04_2.x), Zm04_2.y,Zm04_2.x )) %>%
  mutate(Zm05_1 = ifelse(is.na(Zm05_1.x), Zm05_1.y,Zm05_1.x )) %>%
  mutate(Zm05_2 = ifelse(is.na(Zm05_2.x), Zm05_2.y,Zm05_2.x )) %>%
  mutate(Zm06_1 = ifelse(is.na(Zm06_1.x), Zm06_1.y,Zm06_1.x )) %>%
  mutate(Zm06_2 = ifelse(is.na(Zm06_2.x), Zm06_2.y,Zm06_2.x )) %>%
  mutate(Zm07_1 = ifelse(is.na(Zm07_1.x), Zm07_1.y,Zm07_1.x )) %>%
  mutate(Zm07_2 = ifelse(is.na(Zm07_2.x), Zm07_2.y,Zm07_2.x )) %>%
  mutate(Zm09_1 = ifelse(is.na(Zm09_1.x), Zm09_1.y,Zm09_1.x )) %>%
  mutate(Zm09_2 = ifelse(is.na(Zm09_2.x), Zm09_2.y,Zm09_2.x )) %>%
  mutate(Zm10_1 = ifelse(is.na(Zm10_1.x), Zm10_1.y,Zm10_1.x )) %>%
  mutate(Zm10_2 = ifelse(is.na(Zm10_2.x), Zm10_2.y,Zm10_2.x )) %>%
  mutate(mlg = ifelse(is.na(mlg.x), mlg.y,mlg.x )) 
test2<-test[,c(1:9,49:67)]
test2$mlg<-as.numeric(test2$mlg)

#samples to add to mlg2 
to.add<-anti_join(test2,impute.test)
to.add2<-rbind(to.add,imputed.dat4)

mlg5<-rbind(mlg4,to.add2)
# Sort by vector name [z] then [x]
mlg5<-mlg5[
  with(mlg5, order(time, block,plot,sample)),
  ]

require(plyr)
timereps<-ddply(mlg5, c("time"), summarise, 
                N=length(trt))
timereps

write.csv(mlg5,"Kollars_feedback_mlgs.csv")


