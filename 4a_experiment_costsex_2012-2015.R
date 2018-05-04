# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Experimental test: parasites can increase the cost of sex; A
# Comparison of asexual frequencies in offspring following controlled exposures, 2012-2015

######
# Experiment, 2012-2015
rm(list=ls())
library(aod)
library(reshape)

df<-read.csv("Bins_babies.csv")
df<-subset(df,df$ploidy>0)
summary(df)
attach(df)

# logistic regression with sexual v asexual babies in a tank
# by treatment and year

# subset by ploidy and treatment
aC=subset(df,df$treatment==0 & df$ploidy==3)
aE=subset(df,df$treatment==1 & df$ploidy==3)
sC=subset(df,df$treatment==0 & df$ploidy==2)
sE=subset(df,df$treatment==1 & df$ploidy==2)

# control triploids
aCT=tapply(aC$ploidy,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","triploid")

# exposed triploids
aET=tapply(aE$ploidy,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","triploid")

# control diploids
sCT=tapply(sC$ploidy,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","diploid")

# exposed diploids
sET=tapply(sE$ploidy,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","diploid")

#merge
treatment=rep("exposed",length(LaET$triploid))
ET=cbind(LsET,LaET$triploid,treatment)
names(ET)=c("year","replicate","diploid","triploid","treatment")

treatment=rep("control",length(LaCT$triploid))
CT=cbind(LsCT,LaCT$triploid,treatment)
names(CT)=c("year","replicate","diploid","triploid","treatment")

T=rbind(ET,CT)

lrfit<-glm( cbind(triploid,diploid) ~ treatment*factor(year), 
            family = binomial,data=T)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/40 

# quasi
lrfitq<-glm( cbind(triploid,diploid) ~ treatment*factor(year), 
             family = quasibinomial,data=T)

lrfitq2<-glm( cbind(triploid,diploid) ~ treatment+factor(year), 
              family = quasibinomial,data=T)
anova(lrfitq,lrfitq2,test="F")

lrfitq3<-glm( cbind(triploid,diploid) ~ treatment, 
              family = quasibinomial,data=T)
anova(lrfitq2,lrfitq3,test="F")

lrfitq4<-glm( cbind(triploid,diploid) ~ factor(year), 
              family = quasibinomial,data=T)
anova(lrfitq2,lrfitq4,test="F")

#####
# Repeat: exclude 2012
# Experiment: 2013-2015
Tlat<-subset(T,T$year!=2012)
lrfit<-glm( cbind(triploid,diploid) ~ treatment*factor(year), 
            family = binomial,data=Tlat)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/30

# quasibinomial
lrfitq<-glm( cbind(triploid,diploid) ~ treatment*factor(year), 
             family = quasibinomial,data=Tlat)
lrfitq2<-glm( cbind(triploid,diploid) ~ treatment+factor(year), 
              family = quasibinomial,data=Tlat)
anova(lrfitq,lrfitq2,test="F")

lrfitq3<-glm( cbind(triploid,diploid) ~ treatment, 
              family = quasibinomial,data=Tlat)
anova(lrfitq2,lrfitq3,test="F")

lrfitq4<-glm( cbind(triploid,diploid) ~ factor(year), 
              family = quasibinomial,data=Tlat)
anova(lrfitq2,lrfitq4,test="F")

detach(df)

###########################################
# Calculations

#########################################
#babies in control v exposed
rm(list=ls())
library(aod)
library(reshape)
df<-read.csv("Bins_babies.csv")
df<-subset(df,df$ploidy>0)
summary(df)
attach(df)

# subset by ploidy and treatment
aC=subset(df,df$treatment==0 & df$ploidy==3)
aE=subset(df,df$treatment==1 & df$ploidy==3)
sC=subset(df,df$treatment==0 & df$ploidy==2)
sE=subset(df,df$treatment==1 & df$ploidy==2)

# control triploids first
aCT=tapply(aC$ploidy,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","triploid")

# exposed triploids
aET=tapply(aE$ploidy,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","triploid")

# control diiploids
sCT=tapply(sC$ploidy,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","diploid")

# exposed diploids
sET=tapply(sE$ploidy,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","diploid")

#merge them
treatment=rep("exposed",length(LaET$triploid))
ET=cbind(LsET,LaET$triploid,treatment)
names(ET)=c("year","replicate","diploid","triploid","treatment")

treatment=rep("control",length(LaCT$triploid))
CT=cbind(LsCT,LaCT$triploid,treatment)
names(CT)=c("year","replicate","diploid","triploid","treatment")

T=rbind(ET,CT)
freq=T$triploid/(T$triploid+T$diploid)
sum=T$triploid+T$diploid
Tf=cbind(T,freq,sum)

# sample sizes
mean(Tf$sum)
sd(Tf$sum)/sqrt(length(Tf$sum))

###########2012-2015####################
mean<-tapply(Tf$freq,list(Tf$treatment),mean)
sd<-tapply(Tf$freq,list(Tf$treatment),sd)
se<-sd/sqrt(length(Tf$treatment)/2)

# sample sizes
mean<-tapply(Tf$sum,list(Tf$treatment),mean)
sd<-tapply(Tf$sum,list(Tf$treatment),sd)
se<-sd/sqrt(length(Tf$treatment)/2)

###########2013-2015####################
Tf2<-subset(Tf,Tf$year>2012)

mean<-tapply(Tf2$freq,list(Tf2$treatment),mean)
sd<-tapply(Tf2$freq,list(Tf2$treatment),sd)
se<-sd/sqrt(length(Tf2$treatment)/2)

# sample sizes
mean<-tapply(Tf2$sum,list(Tf2$treatment),mean)
sd<-tapply(Tf2$sum,list(Tf2$treatment),sd)
se<-sd/sqrt(length(Tf2$treatment)/2)

detach(df)
#######################################################
#frequencies of triploid and diploid adults in tanks
#includes males (see methods)

rm(list=ls())
df<-read.csv("Bins_adults_withmales.csv")
df<-subset(df,df$ploidy>0)
summary(df)
attach(df)

# subset by ploidy and treatment
aC=subset(df,df$treatment==0 & df$ploidy==3)
aE=subset(df,df$treatment==1 & df$ploidy==3)
sC=subset(df,df$treatment==0 & df$ploidy==2)
sE=subset(df,df$treatment==1 & df$ploidy==2)

# control triploids first
aCT=tapply(aC$ploidy,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","triploid")

# exposed triploids
aET=tapply(aE$ploidy,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","triploid")

# control diiploids
sCT=tapply(sC$ploidy,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","diploid")

# exposed diploids
sET=tapply(sE$ploidy,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","diploid")

#merge them
treatment=rep("exposed",length(LaET$triploid))
ET=cbind(LsET,LaET$triploid,treatment)
names(ET)=c("year","replicate","diploid","triploid","treatment")

treatment=rep("control",length(LaCT$triploid))
CT=cbind(LsCT,LaCT$triploid,treatment)
names(CT)=c("year","replicate","diploid","triploid","treatment")

T=rbind(ET,CT)
freq=T$triploid/(T$triploid+T$diploid)
sum=T$triploid+T$diploid
Tf=cbind(T,freq,sum)

mean(Tf$freq)
sd(Tf$freq)/sqrt(length(Tf$freq))

# sample sizes
mean(Tf$sum)
sd(Tf$sum)/sqrt(length(Tf$sum))

###########2012-2015####################
# by treatment
mean<-tapply(Tf$freq,list(Tf$treatment),mean)
sd<-tapply(Tf$freq,list(Tf$treatment),sd)
se<-sd/sqrt(length(Tf$treatment)/2)

# sample sizes
mean<-tapply(Tf$sum,list(Tf$treatment),mean)
sd<-tapply(Tf$sum,list(Tf$treatment),sd)
se<-sd/sqrt(length(Tf$treatment)/2)

###########2013-2015####################
Tf2<-subset(Tf,Tf$year>2012)
mean(Tf2$freq)
sd(Tf2$freq)/sqrt(36)

mean<-tapply(Tf2$freq,list(Tf2$treatment),mean)
sd<-tapply(Tf2$freq,list(Tf2$treatment),sd)
se<-sd/sqrt(length(Tf2$treatment)/2)

# sample sizes
mean<-tapply(Tf2$sum,list(Tf2$treatment),mean)
sd<-tapply(Tf2$sum,list(Tf2$treatment),sd)
se<-sd/sqrt(length(Tf2$treatment)/2)

detach(df)
##############################################################
#fold-changes

rm(list=ls())
df<-read.csv("Bins_adults+babies.csv")
df<-subset(df,df$ploidy>0)
summary(df)
attach(df)

# babies
dfB<-subset(df,df$cohort=="baby")
# subset by ploidy and treatment
aC=subset(dfB,dfB$treatment==0 & dfB$ploidy==3)
aE=subset(dfB,dfB$treatment==1 & dfB$ploidy==3)
sC=subset(dfB,dfB$treatment==0 & dfB$ploidy==2)
sE=subset(dfB,dfB$treatment==1 & dfB$ploidy==2)

# control triploids first
aCT=tapply(aC$ploidy,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","triploid")

# exposed triploids
aET=tapply(aE$ploidy,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","triploid")

# control diiploids
sCT=tapply(sC$ploidy,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","diploid")

# exposed diploids
sET=tapply(sE$ploidy,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","diploid")

#merge
treatment=rep("exposed",length(LaET$triploid))
ET=cbind(LsET,LaET$triploid,treatment)
names(ET)=c("year","replicate","diploid","triploid","treatment")

treatment=rep("control",length(LaCT$triploid))
CT=cbind(LsCT,LaCT$triploid,treatment)
names(CT)=c("year","replicate","diploid","triploid","treatment")

T=rbind(ET,CT)
freq=T$triploid/(T$triploid+T$diploid)
sum=T$triploid+T$diploid
Tfbabies=cbind(T,freq,sum)
cohort=rep("baby",length(Tfbabies$treatment))
Tfbabies=cbind(Tfbabies,cohort)

# adults
dfA<-subset(df,df$cohort=="adult")
# subset by ploidy and treatment
aC=subset(dfA,dfA$treatment==0 & dfA$ploidy==3)
aE=subset(dfA,dfA$treatment==1 & dfA$ploidy==3)
sC=subset(dfA,dfA$treatment==0 & dfA$ploidy==2)
sE=subset(dfA,dfA$treatment==1 & dfA$ploidy==2)

# control triploids first
aCT=tapply(aC$ploidy,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","triploid")

# exposed triploids
aET=tapply(aE$ploidy,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","triploid")

# control diiploids
sCT=tapply(sC$ploidy,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","diploid")

# exposed diploids
sET=tapply(sE$ploidy,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","diploid")

#merge them
treatment=rep("exposed",length(LaET$triploid))
ET=cbind(LsET,LaET$triploid,treatment)
names(ET)=c("year","replicate","diploid","triploid","treatment")

treatment=rep("control",length(LaCT$triploid))
CT=cbind(LsCT,LaCT$triploid,treatment)
names(CT)=c("year","replicate","diploid","triploid","treatment")

T=rbind(ET,CT)
freq=T$triploid/(T$triploid+T$diploid)
sum=T$triploid+T$diploid
Tfadults=cbind(T,freq,sum)
cohort=rep("adult",length(Tfadults$treatment))
Tfadults=cbind(Tfadults,cohort)

Tf=rbind(Tfadults,Tfbabies)
Tf$treatment<-factor(Tf$treatment, levels=c("control","exposed"))
Tfadults$treatment<-factor(Tfadults$treatment, levels=c("control","exposed"))
Tfbabies$treatment<-factor(Tfbabies$treatment, levels=c("control","exposed"))

#means
Tfa=mean(Tf$freq[1:48])
Tfa_se=sd(Tf$freq[1:48])/sqrt(48)

TfbE=mean(Tf$freq[49:72])
TfbE_se=sd(Tf$freq[49:72])/sqrt(24)

TfbC=mean(Tf$freq[73:96])
TfbC_se=sd(Tf$freq[73:96])/sqrt(24)

#fold changes
foldC=TfbC/Tfa
se_foldC=foldC*sqrt(((Tfa_se/Tfa)^2)+((TfbC_se/TfbC)^2))

foldE=TfbE/Tfa
se_foldE=foldE*sqrt(((Tfa_se/Tfa)^2)+((TfbE_se/TfbE)^2))


#### exclude 2012
Tf2<-subset(Tf,Tf$year>2012)

#means
Tfa=mean(Tf2$freq[1:36])
Tfa_se=sd(Tf2$freq[1:36])/sqrt(36)

TfbE=mean(Tf2$freq[37:54])
TfbE_se=sd(Tf2$freq[37:72])/sqrt(18)

TfbC=mean(Tf2$freq[55:72])
TfbC_se=sd(Tf2$freq[55:72])/sqrt(18)

#fold changes
foldC=TfbC/Tfa
se_foldC=foldC*sqrt(((Tfa_se/Tfa)^2)+((TfbC_se/TfbC)^2))

foldE=TfbE/Tfa
se_foldE=foldE*sqrt(((Tfa_se/Tfa)^2)+((TfbE_se/TfbE)^2))

detach(df)

###################################################
#ambiguous DNA content
#babies
rm(list=ls())
df<-read.csv("Bins_babies.csv")
summary(df)
attach(df)

# total
d<-subset(df,df$ploidy!="NA")
l<-tapply(d$ploidy,list(d$year,d$tank),length)

#zero
d2<-subset(df,df$ploidy==0)
l2<-tapply(d2$ploidy,list(d2$year,d2$tank),length)
l2[is.na(l2)] <- 0

mean(l2/l)
sd(l2/l)/sqrt(48)

detach(df)

# adults
rm(list=ls())
df<-read.csv("Bins_adults.csv")
summary(df)
attach(df)

# total
d<-subset(df,df$ploidy!="NA")
l<-tapply(d$ploidy,list(d$year,d$tank),length)

#zero
d2<-subset(df,df$ploidy==0)
l2<-tapply(d2$ploidy,list(d2$year,d2$tank),length)
l2[is.na(l2)] <- 0

mean(l2/l)
sd(l2/l)/sqrt(48)

detach(df)
