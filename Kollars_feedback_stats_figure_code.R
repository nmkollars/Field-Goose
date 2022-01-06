#Project title: Disturbance reduces recruitment 
#Authors: N. M. Kollars and J. J. Stachowicz 
#Script by: Nicole Kollars

#Nicole's working directory
setwd("C:/Users/nmkol/OneDrive - UC Davis/Field_Goose/drafts/Feedback_MS")

#########
#needed libraries and functions 
#########

library(RClone)
library(tidyverse)
library(reshape2)
library(scales)
library(AER)
library(multcomp)

#function to add minior tick marks 
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

##########
#"Kollars_feedback_mlgs.csv" has the imputed and adjusted MLGS. It is the output file 
#of the code in Appendix_S4. 
##########
dat<-read.csv("Kollars_feedback_mlgs.csv")
dat<-dat[,-1]

#creates dataframes to use in RClone
#vector of populations
popvec<-dat[,1]
#vector of identifying info
idvec<-dat[,1:9]
#microsat data for RClone functions  
dat2<-dat[,10:27]

##########
#clonal summary stats
##########
clones<-clonal_index(dat2,vecpop=popvec)
#writing csv file so that it is converted to a dataframe
write.csv(clones,"clones.csv")
#in excel, added the identifying information (block, trt, etc)
#for each plot. T. 
#The modified file is "Kollars_feedback_clones.csv" 
clones<-read.csv("Kollars_feedback_clones.csv")
#selecting only variables of interest
clones<-clones[,1:7]

#adjust dataframe from long to wide
clones.T1 <- clones %>%
  filter(time == "T1") %>%
  rename(T1.G=G,T1.N=N) %>%
  select_at(vars(block,trt,T1.N,T1.G))

clones.T5 <- clones %>%
  filter(time == "T5") %>%
  rename(T5.N=N,T5.G=G) %>%
  select_at(vars(block,trt,T5.N,T5.G))

clones2<-full_join(clones.T1, clones.T5)
#create net change in # of genotypes variable 
clones2$geno<-clones2$T5.G-clones2$T1.G

##now add the # of colonization and extinction events
#to do this, start by making dataset of relative abundance numbers for each
#mlg at T1 and T5

dat3<- dat %>%
  select_at(vars(block,time,trt,mlg))

require(plyr)
df_total=data.frame()
for(i in unique(dat3$trt)){
  treatment<-dat3[dat3$trt==i,]
  
  for(plotnum in as.character(unique(treatment$block))){
    ten<-treatment[treatment$block==plotnum,]
    
    #timepoint 1
    ten1<-count(ten[ten$time=="1",],vars="mlg")
    ten1$ra<-ten1$freq/sum(ten1$freq)
    ten1$plot<-rep(plotnum,length(ten1$ra))
    ten1$time<-rep(1,length(ten1$ra))
    
    #timepoint 5
    ten5<-count(ten[ten$time=="5",],vars="mlg")
    ten5$ra<-ten5$freq/sum(ten5$freq)
    ten5$plot<-rep(plotnum,length(ten5$ra))
    ten5$time<-rep(5,length(ten5$ra))
    
    #combining all of the timepoints together for the plot
    tenall<-rbind(ten1,ten5)               
    
    #inputting zeros 
    
    #gives the mlg names that need to be repsented at each timepoint
    ten_genos<-unique(tenall$mlg)
    richness<-length(ten_genos)
    
    tenra<-data.frame(select_at(tenall,vars("mlg","ra","plot","time")))
    
    toexpand<-tenra
    
    expandto<-as.data.frame(matrix(data=0,nrow=richness*2,ncol=3))
    names(expandto)<-c("mlg","plot","time")
    
    times<-c(1,5)
    
    #Copy over the labeling for the categorical columns
    expandto$mlg<-rep(ten_genos)
    expandto$plot<-rep(plotnum)
    expandto$time<-rep(times,each=length(ten_genos))
    
    allw0s<-merge(toexpand,expandto,all=TRUE)
    allw0s$ra<-replace(allw0s$ra,is.na(allw0s$ra),0)
    df <- data.frame(allw0s)
    df$i<-i
    df_total <- rbind(df_total,df)
  }
}  

detach("package:plyr",unload=TRUE)

rel<- df_total %>%
  select_at(vars(plot,i,time,mlg,ra)) %>%
  rename("block"="plot","trt"="i")
rel2 <- spread(rel, time, ra) 
colnames(rel2)<-c("block","trt","mlg","T1","T5")  

#########
#Colonization/Extinction events 
#########

#raw number of colonization and extinction events 
#for each treatment, and for each unique mlg
#colonization = # of times that ra for T1 =0 but T5>0
#extinction = # of times that ra for T1>0 and T5 = 0 

###colonization events 
col=data.frame()
for(j in unique(rel$trt)){
  trt<-rel[rel$trt==j,]
  colonization<-trt[trt$time=="1",][trt[trt$time=="1",]$ra==0,]
  df <- data.frame(colonization)
  col <- rbind(col,df)
}
require(plyr)
sumcol<- ddply(col, c("trt","block"), summarise, 
               N=length(ra))
names(sumcol)[3] <- "col"

###extinction events
ext=data.frame()
for(j in unique(rel$trt)){
  trt<-rel[rel$trt==j,]
  extinction<-trt[trt$time=="5",][trt[trt$time=="5",]$ra==0,]
  df <- data.frame(extinction)
  ext <- rbind(ext,df)
}

sumext<- ddply(ext, c("trt","block"), summarise, 
               N=length(ra))
names(sumext)[3] <- "ext"

colext<-full_join(sumext,sumcol)
colext[is.na(colext)] <- 0

#inputting zeros 
toexpand<-colext
expandto<-as.data.frame(matrix(data=0,nrow=36,ncol=4))
names(expandto)<-c("block","trt","col","ext")

#Copy over the labeling for the categorical columns
expandto$block<-c(10,10,10,10,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9)
expandto$block<-as.character(expandto$block)
expandto<-arrange(expandto,block)
expandto$trt<-rep(c(c("0X","1X","2X","4X"),#10
                    c("0X","1X","2X","4X"),#2
                    c("0X","1X","2X","4X"),#3
                    c("0X","1X","2X","4X"),#4
                    c("0X","1X", "2X","4X"),#5
                    c("0X","1X","2X","4X"),#6
                    c("0X","1X","2X","4X"),#7
                    c("0X","1X","2X","4X"),#8
                    c("0X","1X","2X","4X")))#9
dat4<-left_join(expandto,toexpand,by=c("block","trt"))
detach("package:plyr",unload=TRUE)

dat4<-dat4%>%
  rename(col=col.y,ext=ext.y)%>%
  select_at(vars(block,trt,col,ext))%>%
  mutate(col = ifelse(is.na(col), 0, col))%>%
  mutate(ext = ifelse(is.na(ext), 0, ext))

#merge #col and #ext with the rest of the clones data 
all.dat<-merge(dat4,clones2)

require(plyr)
#total number of colonization and extinction events
colSums(all.dat[,3:4])

##add in the shoot count data
shoots<-read.csv("Kollars_feedback_shoots.csv",check.names=FALSE)
##and merge with the rest of the data  
all.dat<-merge(all.dat,shoots)

################
#SHOOT DENSITIES
#################
#remove the flowering shoots here
all<-all.dat[,-(31:34)]
all<-melt(all,id.vars=c(1:10))
all$variable<-as.Date(all$variable,"%m/%d/%Y")
all<-na.omit(all)
require(plyr)

##multiplying by 9 gives the estimated shoot densities of entire genetics sampling plot
sumall<- ddply(all, c("variable","trt"), summarise, 
               N=length(value),
               mean=mean(value*9),
               sd=sd(value*9),
               se=sd/sqrt(N))
sumall
pd<-position_dodge(0.1)
custom_breaks<-seq(0,900,25)
clip<-as.Date(c("01/14/2017","02/10/2017","03/09/2017","04/29/2017","02/28/2018","03/27/2018"),"%m/%d/%Y")
clip2<-as.Date(c("01/07/2017","02/03/2017","03/02/2017","04/22/2017","02/21/2018","03/20/2018"),"%m/%d/%Y")

P1<-ggplot(sumall,aes(variable,y=mean, color=trt))+
  geom_point(position=pd, size=2)+
  geom_ribbon(aes(x=variable, ymin=(mean-se), ymax=(mean+se),
                  fill=trt, color=trt),alpha=0.2) +
  scale_y_continuous(limits=c(50,825),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,2, inverse = TRUE))+
  scale_x_date(labels = date_format("%b"),date_breaks = "1 month",
               limits=c(as.Date("10/15/2016","%m/%d/%Y"),
                        as.Date("10/31/2018","%m/%d/%Y")))+
  scale_color_manual("Disturbance",labels = c("0X", "1X", "2X", "4X"),
                     values=c("#000000","#999999","#E69F00","#56B4E9"))+ 
  scale_fill_manual("Disturbance",labels = c("0X", "1X", "2X", "4X"),
                    values=c("#000000","#999999","#E69F00","#56B4E9"))+
  geom_segment(x=as.Date("04/01/2018","%m/%d/%Y"), y=0, 
           xend=as.Date("06/01/2018","%m/%d/%Y"), yend=0,
           col="red",lwd=10) +
  ylab("Shoot density (# shoots per 0.36 m^2)\n") +xlab("\nMonth (2016 - 2018)")+
  theme_classic()+
  theme(legend.position= c(0.84,0.175),
        axis.text.x = element_text(angle = 45,hjust=1),
        text=element_text(size=8))+
  geom_vline(xintercept=clip,lwd=0.5,lty=2,color="gray80")+
  annotate("text", x =clip2, y = 140, 
           label = c("4X (1 of 4)", "1X, 2X (1 of 2), 4X (2 of 2)",
                     "2X (2 of 2), 4X (3 of 4)",
                    "4X (4 of 4)", "1X;  2X (1 of 2)","2X (2 of 2)"),angle=90,size=2)

#Figure 1 inset 
all.dat$flo.dens.apr<-all.dat$no.flowers.april.2018/all.dat[,26]
all.dat$flo.dens.may<-all.dat$no.flowers.may.2018/all.dat[,27]
all.dat$flo.dens.jun<-all.dat$no.flowers.june.2018/all.dat[,28]
all.flo<-all.dat[,c(1:10,35:37)]
dat10<-melt(all.flo,id.vars=c(1:10))
names(dat10)[12] <- "dens"
sumshoots<- ddply(dat10, c("variable","trt"), summarise, 
                  N=length(dens),
                  mean=mean(dens,na.rm=TRUE),
                  sd=sd(dens,na.rm=TRUE),
                  se=sd/sqrt(N))
sumshoots
p2<-ggplot(sumshoots,aes(x=as.numeric(variable),y=mean,col=trt,fill=trt))+
  geom_point(cex=2)+
  geom_ribbon(aes(x=as.numeric(variable), ymin=(mean-se), ymax=(mean+se),
                  fill=trt, color=trt),alpha=0.2)+
  scale_color_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                     values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_fill_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                    values=c("#000000","#999999","#E69F00","#56B4E9"))+
  theme_classic()+
  xlab("Spring 2018")+ylab("Proportion of shoots flowering")+
  scale_x_continuous(limits=c(0.8,3.2),breaks=c(1,2,3),
                     labels=c("April","May","June"))+
  theme(text=element_text(size=8),
        legend.position="none")
P1 + annotation_custom(ggplotGrob(p2), xmin = as.Date("05/01/2017","%m/%d/%Y"), 
                       xmax = as.Date("02/01/2018","%m/%d/%Y"), 
                       ymin = 530, ymax = 850)
  

ggsave("Kollars_Stach_Ecol_Figure1_revised.tiff",dpi=600,height=5,width=6,units="in",compression="lzw")


##########
#Initial Richness vs Final Richness
##########

#draw lines on plot with m1 so that curves fit 
m1<-glm(T5.G~T1.G*trt,family="poisson",data=all.dat)
plot(m1)
dispersiontest(m1)#not overdispersed 
summary(m1)
ndat<-data.frame()
for (trt in unique(all.dat$trt)){
  trtr<-all.dat[all.dat$trt==trt,]
  rand<-with(trtr, data_frame(T1.G = seq(min(T1.G), max(T1.G),
                                         length = length(T1.G))))
  df<-data.frame(cbind(trt,rand))
  ndat<-rbind(ndat,df)
}

ndat$predict<- predict(m1,newdata=ndat, type="response")
ndat<- ndat[with(ndat, order(trt,T1.G)), ]

ilink <- family(m1)$linkinv
ndat<- bind_cols(ndat, setNames(as_tibble(predict(m1,newdata=ndat,se.fit = TRUE)[1:2]),
                                c('fit_link','se_link')))

## create the interval and backtransform
ndat <- mutate(ndat,
               fit_resp  = ilink(fit_link),
               right_upr = ilink(fit_link + (2*se_link)),
               right_lwr = ilink(fit_link - (2*se_link)))
#Plot
custom_breaks <- seq(-20, 20, 1)
set.seed(20)
ggplot(ndat,aes(x=T1.G,y=predict,color=trt,fill=trt))+
  geom_abline(intercept = 0, slope = 1,lwd=0.8,lty=2,color="grey80")+
  geom_ribbon(aes(x=T1.G,ymin = right_lwr, ymax = right_upr),
              alpha = 0.1)+
  geom_jitter(data=all.dat, aes(x=T1.G, y=T5.G,color=trt),cex=2)+
  geom_line(size = 1.5)+
  theme_bw()+
  scale_color_manual("Clipping frequency",labels = c("0X", "1X", "2X", "4X (recovery)"),
                     values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_fill_manual("Clipping frequency",labels = c("0X", "1X", "2X", "4X (recovery)"),
                    values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_y_continuous(limits=c(0,20),breaks = custom_breaks,
                     labels = every_nth(custom_breaks, 2, inverse = TRUE)) +
  scale_x_continuous(limits=c(0.6,10.5),breaks = custom_breaks,
                     labels = every_nth(custom_breaks, 1, inverse = TRUE))+
  facet_wrap(~trt)+
  xlab("Initial # of genotypes")+ylab("Final # of genotypes ")+
  theme(text=element_text(size=8),
        legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("feedback_richness.tiff",dpi=600,height=5,width=5,units="in",compression="lzw")

#############
##Stats######
#############
library(car)
##"clipping" is data for 0X, 1X, and 2X only
clipping<-all.dat[all.dat$trt!="4X",]
clipping$trt<-as.factor(clipping$trt)

##"recovery" is data comparing 0X and 4X only 
recovery<-all.dat[all.dat$trt!="1X" & all.dat$trt!="2X", ]
recovery$trt<-as.numeric(as.factor(recovery$trt))
recovery$trt[recovery$trt=="4"] <- "2"
recovery$trt<-as.factor(recovery$trt)

#CLIPPING: Raw change in the # genos
m1<-glm(geno~T1.G+trt+T1.G*trt,family="gaussian",data=clipping)
summary(m1)
m2<-glm(geno~T1.G+trt,family="gaussian",data=clipping)
summary(m2)
m3<-glm(geno~trt,family="gaussian",data=clipping)
summary(m3)
m4<-glm(geno~T1.G,family="gaussian",data=clipping)
summary(m4)
Anova(m1,test="F")
#m3 is best fit 
#model fit
plot(m3)
#coefficent estimates
summary(m3)
#post-hoc
summary(glht(m3, mcp(trt="Tukey")))

#RECOVERY:raw change in #genos  
m1<-glm(geno~T1.G*trt,family="gaussian",data=recovery)
m2<-glm(geno~T1.G+trt,family="gaussian",data=recovery)
m3<-glm(geno~trt,family="gaussian",data=recovery)
m4<-glm(geno~T1.G,family="gaussian",data=recovery)
Anova(m1,test="F")
#model fit
plot(m4)
#coefficent estimates
summary(m4)

#CLIPPING: # Colonization
m1<-glm(col~T1.G*trt,family="poisson",data=clipping)
dispersiontest(m1)#not overdispersed#not overdispersed
summary(m1)
m2<-glm(col~T1.G+trt,family="poisson",data=clipping)
summary(m2)
m3<-glm(col~trt,family="poisson",data=clipping)
summary(m3)
m4<-glm(col~T1.G,family="poisson",data=clipping)
#interaction test 
Anova(m1)
#m2 is best fit
#model fit
plot(m2)
#coefficent estimates
summary(m2)
#posthoc
summary(glht(m2, mcp(trt="Tukey")))

#RECOVERY: Colonization
#recovery
m1.1<-glm(col~T1.G*trt,family="poisson",data=recovery)
dispersiontest(m1.1)#notoverdispersed 
m2.1<-glm(col~T1.G+trt,family="poisson",data=recovery)
m3.1<-glm(col~trt,family="poisson",data=recovery)
m4.1<-glm(col~T1.G,family="poisson",data=recovery)
#interaction test 
Anova(m1.1)
#model fit
plot(m2.1)
#coefficent estimates
summary(m2.1)

#CLIPPING: # Extinctions
m1<-glm(ext~T1.G*trt,family="poisson",data=clipping)
dispersiontest(m1)#not overdispersed 
summary(m1)
m2<-glm(ext~T1.G+trt,family="poisson",data=clipping)
summary(m2)
m3<-glm(ext~trt,family="poisson",data=clipping)
summary(m3)
m4<-glm(ext~T1.G,family="poisson",data=clipping)
summary(m4)
Anova(m1)
#m4 is best fit 
#model fit
plot(m4)
#coefficent estimates
summary(m4)

#RECOVERY: Extinctions
#recovery 
m1.1<-glm(ext~T1.G*trt,family="poisson",data=recovery)
dispersiontest(m1.1)#not overdispersed
m2.1<-glm(ext~T1.G+trt,family="poisson",data=recovery)
m3.1<-glm(ext~trt,family="poisson",data=recovery)
m4.1<-glm(ext~T1.G,family="poisson",data=recovery)
Anova(m1.1)
#model fit
plot(m4.1)
#coefficent estimates
summary(m4.1)

#CLIPPING: Flowers 
clipping$flo.dens.apr<-clipping$no.flowers.april.2018/clipping[,26]
clipping$flo.dens.may<-clipping$no.flowers.may.2018/clipping[,27]
clipping$flo.dens.jun<-clipping$no.flowers.june.2018/clipping[,28]
clip.flo<-clipping[,c(1:10,35:37)]
dat10<-melt(clip.flo,id.vars=c(1:10))
names(dat10)[12] <- "dens"
#(permutation analysis)
numperm=1000
print(anova(lm(dens ~ variable + trt*T1.G, data = dat10)))

retval <- matrix(0, ncol = 4, nrow = numperm + 1)
retval[1, ] <- anova(lm(dens ~ variable + trt*T1.G, data = dat10))[1:4, 4]

for (i in 1:numperm) {
  dat10$permtime <- sample(dat10$variable, length(dat10$variable), replace = F)
  dat10$permtrt <- sample(dat10$trt, length(dat10$trt), replace = F)
  dat10$permG<-sample(dat10$T1.G,length(dat10$T1.G),replace=F)
  retval[i + 1, ] <- anova(lm(dens ~ permtime + permtrt*permG, data = dat10))[1:4,4]
}

par(mfrow = c(2, 2))
probALL = NULL
###
for (i in 1:4) {
  prob <- sum(retval[-1, i] >= retval[1, i])/(dim(retval)[1] - 1)
  probALL <- rbind(probALL, prob)
  hist(retval[-1, i],main = paste("Observed F=", round(retval[1, i],2),";", "p =",prob),
       xlab = "Simulated F-statistic")
  points(c(retval[1, i], retval[1, i]), c(0, numperm), type = "l", lwd = 2, 
         col = "gray",lty=2)
}

#RECOVERY: Flowering
recovery$flo.dens.apr<-recovery$no.flowers.april.2018/recovery[,26]
recovery$flo.dens.may<-recovery$no.flowers.may.2018/recovery[,27]
recovery$flo.dens.jun<-recovery$no.flowers.june.2018/recovery[,28]
reco.flo<-recovery[,c(1:10,35:37)]
dat10<-melt(reco.flo,id.vars=c(1:10))
names(dat10)[12] <- "dens"

#stats (permutation analysis)
numperm=1000
print(anova(lm(dens ~ variable + trt*T1.G, data = dat10)))

retval <- matrix(0, ncol = 4, nrow = numperm + 1)
retval[1, ] <- anova(lm(dens ~ variable + trt*T1.G, data = dat10))[1:4, 4]

for (i in 1:numperm) {
  dat10$permtime <- sample(dat10$variable, length(dat10$variable), replace = F)
  dat10$permtrt <- sample(dat10$trt, length(dat10$trt), replace = F)
  dat10$permG<-sample(dat10$T1.G,length(dat10$T1.G),replace=F)
  retval[i + 1, ] <- anova(lm(dens ~ permtime + permtrt*permG, data = dat10))[1:4,4]
}

par(mfrow = c(2, 2))
probALL = NULL
###
for (i in 1:4) {
  prob <- sum(retval[-1, i] >= retval[1, i])/(dim(retval)[1] - 1)
  probALL <- rbind(probALL, prob)
  hist(retval[-1, i],main = paste("Observed F=", round(retval[1, i],2),";", "p =",prob),
       xlab = "Simulated F-statistic")
  points(c(retval[1, i], retval[1, i]), c(0, numperm), type = "l", lwd = 2, 
         col = "gray",lty=2)
}

###############
####FIGURES####
###############

###FIGURE 2: NET CHANGE IN GENOS###
require(plyr)
sumgeno<- ddply(all.dat, c("trt"), summarise, 
                N=length(geno),
                mean=mean(geno,na.rm=TRUE),
                sd=sd(geno,na.rm=TRUE),
                se=sd/sqrt(N))
sumgeno
custom_breaks<-seq(-5,5,1)
ggplot(sumgeno,aes(x=as.numeric(as.factor(trt)),y=mean,fill=trt,color=trt))+
  geom_hline(yintercept = 0,lty=2)+
  geom_vline(xintercept=3.5,lty=1)+
  geom_jitter(data=all.dat,aes(x=as.numeric(as.factor(trt)),
    y=geno,color=trt,fill=trt,alpha=0.9),
    cex=2,width=0.2,height=0.1)+
  geom_point(cex=3,color="black")+
  geom_errorbar(aes(ymin=mean+se,ymax=mean-se),width=0.2,
                position=position_dodge(0.9),lwd=1,color="black")+
  theme_classic()+ xlab(NULL)+ ylab("Net change in # of genotypes")+
  scale_color_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                     values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_fill_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                    values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_x_continuous(limits=c(0.8,4.2),breaks=c(1,2,3,4),labels=c("0X","1X","2X","4X"))+
  scale_y_continuous(limits=c(-4.25,3.25),breaks = custom_breaks,
                     labels = every_nth(custom_breaks, 1, inverse = TRUE))+ 
  theme(text=element_text(size=8),legend.position="none")+
  geom_text(aes(x=1,y=3,label="A"),color="black",size=3)+
  geom_text(aes(x=2,y=3,label="B"),color="black",size=3)+
  geom_text(aes(x=3,y=3,label="B"),color="black",size=3)
ggsave("Kollars_Stach_Ecol_Figure2_Revised.tiff",dpi=600,height=3,width=3,units="in",compression="lzw")

###Figure3: # COLONIZATIONS AND EXTINCTIONS###
sum.col<- ddply(all.dat, c("trt"), summarise, 
                N=length(col),
                mean=mean(col,na.rm=TRUE),
                sd=sd(col,na.rm=TRUE),
                se=sd/sqrt(N))
sum.col
custom_breaks<-seq(-10,10,1)
clip_gain<-ggplot(sum.col,aes(x=as.numeric(as.factor(trt)),y=mean,color=trt))+
  geom_vline(xintercept=3.5,lty=1)+
  geom_jitter(data=all.dat,aes(x=as.numeric(as.factor(trt)),y=col,color=trt,fill=trt,alpha=0.9),
              cex=2,width=0.2,height=0.1)+
  geom_errorbar(aes(ymin=mean+se,ymax=mean-se),width=0.2,
                position=position_dodge(0.9),lwd=1,color="black")+
  geom_point(cex=3,color="black")+
  geom_text(aes(x=1,y=6,label="A"),color="black",size=3)+
  geom_text(aes(x=2,y=6,label="AB"),color="black",size=3)+
  geom_text(aes(x=3,y=6,label="B"),color="black",size=3)+
  theme_classic()+ xlab(NULL)+ ylab("# of genotypes gained")+
  scale_color_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                     values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_fill_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                    values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_x_continuous(limits=c(0.8,4.2),breaks=c(1,2,3,4),labels=c("0X","1X","2X","4X"))+
  scale_y_continuous(limits=c(-0.2,6),breaks = custom_breaks,
                     labels = every_nth(custom_breaks, 1, inverse = TRUE))+ 
  theme(text=element_text(size=8),legend.position="none")
clip_gain
sum.ext<- ddply(all.dat, c("trt"), summarise, 
                N=length(ext),
                mean=mean(ext,na.rm=TRUE),
                sd=sd(ext,na.rm=TRUE),
                se=sd/sqrt(N))
sum.ext
custom_breaks<-seq(-10,0,1)
clip_loss<-ggplot(sum.ext,aes(x=as.numeric(as.factor(trt)),y=-mean,color=trt))+
  geom_hline(yintercept=0,lty=2)+
  geom_vline(xintercept=3.5,lty=1)+
  geom_jitter(data=all.dat,aes(x=as.numeric(as.factor(trt)),y=-ext,
                                color=trt,fill=trt,alpha=0.9),
              cex=2,width=0.2,height=0.1)+
  geom_errorbar(aes(ymin=-mean+se,ymax=-mean-se),width=0.2,
                position=position_dodge(0.9),lwd=1, color="black")+
  geom_point(cex=3,color="black")+
  theme_classic()+ xlab(NULL)+ ylab("# of genotypes lost")+
  scale_color_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                     values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_fill_manual(NULL,labels = c("0X", "1X", "2X","4X"),
                    values=c("#000000","#999999","#E69F00","#56B4E9"))+
  scale_x_continuous(limits=c(0.8,4.2),breaks=c(1,2,3,4),labels=c("0X","1X","2X","4X"))+
  scale_y_continuous(limits=c(-7.1,0.3),breaks = custom_breaks,
                     labels = every_nth(custom_breaks, 1, inverse = TRUE))+ 
  theme(text=element_text(size=8),legend.position="none")
clip_loss

library(cowplot)
plot_grid(clip_gain, clip_loss, labels = c("A)","B)"),label_size=8)
ggsave("Kollars_Stach_Ecol_Figure3_Revised.tiff",dpi=600,height=3,width=6,units="in",compression="lzw")



