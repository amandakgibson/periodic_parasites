# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Temporal shifts in relative infection prevalence, part B
# Comparison of Microphallus infection rates in the field in 2001-2005 and 2012-2016

######
# Field, all sites, females
rm(list=ls())
library(aod)
library(lme4)
library(geepack)
library(lmtest)

df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & length>0 & ploidy>0)
summary(df)
attach(df)

geefit1<-geeglm(mic~site*factor(ploidy)+factor(period)*factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")

geefit2<-geeglm(mic~site+factor(period)*factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit1,geefit2)  

geefit3<-geeglm(mic~site+factor(period)+factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit2,geefit3)

#results
anova(geefit2,test="LRT")
summary(geefit2)

detach(df)

######
# Variance in infection rates
# Field, all sites, females
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0  & ploidy>0)
summary(df)

attach(df)
library(reshape)

# asexual infection rates
asex<-subset(df,df$ploidy==3)
aT<-t(tapply(asex$mic,list(asex$year,asex$site),length))
aI<-t(tapply(asex$mic,list(asex$year,asex$site),sum))
aF<-aI/aT
longaF<-melt(aF,id=c("year","frequency"))

# sexual
sex<-subset(df,df$ploidy==2)
sT<-t(tapply(sex$mic,list(sex$year,sex$site),length))
sI<-t(tapply(sex$mic,list(sex$year,sex$site),sum))
sF<-sI/sT
longsF<-melt(sF,id=c("year","frequency"))

longsF<-na.omit(longsF)
longaF<-na.omit(longaF)

#means
mean = c(mean(longsF$value),mean(longaF$value))
se = c((sd(longsF$value)/sqrt(length(longsF$value))),sd(longaF$value)/sqrt(length(longaF$value)))

# Coefficient of variation
cvS<-sd(longsF$value)/mean(longsF$value)
cvA<-sd(longaF$value)/mean(longaF$value)

# bootstrap
# sexual
m = 10000 
cvBootS=rep(0,m)

for (i in 1:m) { 
  dboot<-longsF[sample(1:50,50, replace=TRUE), ]
  cvBootS[i]<-sd(dboot$value)/mean(dboot$value)
}

cvBootS = sort(cvBootS,decreasing=FALSE) 
lowS<-cvBootS[m*0.025]
highS<-cvBootS[m*0.975] 

# asexual
cvBootA=rep(0,m)

for (i in 1:m) { 
  dboot<-longaF[sample(1:46,46, replace=TRUE), ]
  cvBootA[i]<-sd(dboot$value)/mean(dboot$value)
}

cvBootA = sort(cvBootA,decreasing=FALSE) 
lowA<-cvBootA[m*0.025]
highA<-cvBootA[m*0.975] 

######
# Repeat GEE: four shared sites
# Field, four sites, females
rm(list=ls())
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df2<-subset(df2,df2$site!="swend" & df2$site!="jms")
df2$site <- factor(df2$site)
df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & length>0 & ploidy>0)
summary(df)
attach(df)

geefit1<-geeglm(mic~site*factor(period) + site*factor(ploidy)+factor(period)*factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")

summary(geefit1)
anova(geefit1,test="LRT")

geefit2<-geeglm(mic~site*factor(period) + factor(period)*factor(ploidy)+ length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit2, test="LRT")

geefit3<-geeglm(mic~site +factor(period)*factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit2,geefit3) 

geefit4<-geeglm(mic~site +factor(period)+factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit3,geefit4) 

anova(geefit3,test="LRT")

detach(df)

######
# Repeat GEE: just sexual females
# Field, all sites, sexual females
rm(list=ls())
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & length>0 & ploidy==2)
summary(df)
attach(df)

geefit1<-geeglm(mic~site+ factor(period) + length, data=df, family=binomial,id=year,corstr="exchangeable")

anova(geefit1,test="LRT")
summary(geefit1)

detach(df)

######
# Repeat GEE: just asexual females
# Field, all sites, asexual females
rm(list=ls())
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & length>0 & ploidy==3)
summary(df)

attach(df)

geefit1<-geeglm(mic~site+ factor(period) + length, data=df, family=binomial,id=year,corstr="exchangeable")

anova(geefit1,test="LRT")
summary(geefit1)

detach(df)

######
# Repeat GEE: males and females
# Field, all sites, males + females
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016_withmales.csv")
df<-rbind(df1,df2)
df<-subset(df,length>0 & ploidy>0)
summary(df)

attach(df)

geefit1<-geeglm(mic~site*factor(ploidy)+factor(period)*factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
geefit2<-geeglm(mic~site+factor(period) * factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit1,geefit2) 

anova(geefit2,test="LRT")

detach(df)

######
# Repeat GEE: males and females, four shared sites
# Field, four sites, males + females
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016_withmales.csv")
df2<-subset(df2,df2$site!="swend" & df2$site!="jms")
df2$site <- factor(df2$site)
df<-rbind(df1,df2)
df<-subset(df,length>0 & ploidy>0)
summary(df)
attach(df)

geefit1<-geeglm(mic~site+factor(period)*factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
geefit2<-geeglm(mic~site+factor(period) + factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit1,geefit2) 
geefit3<-geeglm(mic~site*factor(period)+factor(period)*factor(ploidy) + length, data=df, family=binomial,id=year,corstr="exchangeable")
anova(geefit3)

detach(df)

#####################################################
# Calculations

#######################################################
# Microphallus prevalence in sexual and asexual females
# 6 sites
rm(list=ls())

df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0  & ploidy>0)
summary(df)
attach(df)

# infected and uninfected by year/site
# asexual
asex<-subset(df,df$ploidy==3)
aT<-t(tapply(asex$mic,list(asex$year,asex$site),length))
aI<-t(tapply(asex$mic,list(asex$year,asex$site),sum))
aF<-aI/aT
longaF<-melt(aF,id=c("year","frequency"))

# infected and uninfected by year/site
# sexual
sex<-subset(df,df$ploidy==2)
sT<-t(tapply(sex$mic,list(sex$year,sex$site),length))
sI<-t(tapply(sex$mic,list(sex$year,sex$site),sum))
sF<-sI/sT
longsF<-melt(sF,id=c("year","frequency"))

longsF<-na.omit(longsF)
longaF<-na.omit(longaF)

#means
mean = c(mean(longsF$value),mean(longaF$value))
se = c((sd(longsF$value)/sqrt(length(longsF$value))),sd(longaF$value)/sqrt(length(longaF$value)))


# sample sizes
# per replicate 
# sexual
longsT<-melt(sT)
longsT<-na.omit(longsT)
s<-mean(longsT$value)
s_se<-sd(longsT$value)/sqrt(length(longsT$value))

#asexual
longaT<-melt(aT)
longaT<-na.omit(longaT)
a<-mean(longaT$value)
a_se<-sd(longaT$value)/sqrt(length(longaT$value))

# just 2001 to 2005
dfp1<-subset(df,df$period==0)

# infected and uninfected by year/site
# asexual
asex<-subset(dfp1,dfp1$ploidy==3)
aT<-t(tapply(asex$mic,list(asex$year,asex$site),length))
aI<-t(tapply(asex$mic,list(asex$year,asex$site),sum))
aF<-aI/aT
longaF<-melt(aF,id=c("year","frequency"))

# infected and uninfected by year/site
# sexual
sex<-subset(dfp1,dfp1$ploidy==2)
sT<-t(tapply(sex$mic,list(sex$year,sex$site),length))
sI<-t(tapply(sex$mic,list(sex$year,sex$site),sum))
sF<-sI/sT
longsF<-melt(sF,id=c("year","frequency"))

longsF<-na.omit(longsF)
longaF<-na.omit(longaF)

#means
mean = c(mean(longsF$value),mean(longaF$value))
se = c((sd(longsF$value)/sqrt(length(longsF$value))),sd(longaF$value)/sqrt(length(longaF$value)))

# sample sizes
# per replicate 
# sexual
longsT<-melt(sT)
longsT<-na.omit(longsT)
s<-mean(longsT$value)
s_se<-sd(longsT$value)/sqrt(length(longsT$value))

#asexual
longaT<-melt(aT)
longaT<-na.omit(longaT)
a<-mean(longaT$value)
a_se<-sd(longaT$value)/sqrt(length(longaT$value))

# mean 2001-2005
earlyA<-subset(longaF,longaF$X2<2005)
earlyS<-subset(longsF,longsF$X2<2005)
mean = c(mean(earlyS$value),mean(earlyA$value))
se = c((sd(earlyS$value)/sqrt(length(earlyS$value))),sd(earlyA$value)/sqrt(length(earlyA$value)))

# mean 2001-2005
A05<-subset(longaF,longaF$X2==2005)
S05<-subset(longsF,longsF$X2==2005)
mean = c(mean(S05$value),mean(A05$value))
se = c((sd(S05$value)/sqrt(length(S05$value))),sd(A05$value)/sqrt(length(A05$value)))


###########
# 4 sites
rm(list=ls())
detach(df)
df1<-read.csv("../Field-2001-2005.csv")
df2<-read.csv("../Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0  & ploidy>0 & df$site!="swend" & df$site!="jms")
summary(df)

attach(df)

# infected and uninfected by year/site
# asexual
asex<-subset(df,df$ploidy==3)
aT<-t(tapply(asex$mic,list(asex$year,asex$site),length))
aI<-t(tapply(asex$mic,list(asex$year,asex$site),sum))
aF<-aI/aT
longaF<-melt(aF,id=c("year","frequency"))

# infected and uninfected by year/site
# sexual
sex<-subset(df,df$ploidy==2)
sT<-t(tapply(sex$mic,list(sex$year,sex$site),length))
sI<-t(tapply(sex$mic,list(sex$year,sex$site),sum))
sF<-sI/sT
longsF<-melt(sF,id=c("year","frequency"))

longsF<-na.omit(longsF)
longaF<-na.omit(longaF)

#means
mean = c(mean(longsF$value),mean(longaF$value))
se = c((sd(longsF$value)/sqrt(length(longsF$value))),sd(longaF$value)/sqrt(length(longaF$value)))

# sample sizes
# per replicate 
# sexual
longsT<-melt(sT)
longsT<-na.omit(longsT)
s<-mean(longsT$value)
s_se<-sd(longsT$value)/sqrt(length(longsT$value))

#asexual
longaT<-melt(aT)
longaT<-na.omit(longaT)
a<-mean(longaT$value)
a_se<-sd(longaT$value)/sqrt(length(longaT$value))

detach(df)

##################################
# overall prevalence across time periods
# 6 sites
rm(list=ls())
detach(df)
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0  & ploidy>0)
summary(df)
attach(df)

T<-t(tapply(df$mic,list(df$year,df$site),length))
I<-t(tapply(df$mic,list(df$year,df$site),sum))
F<-I/T
longF<-melt(F,id=c("year","frequency"))
longF<-na.omit(longF)

p1<-subset(longF,longF$X2<2012)
p2<-subset(longF,longF$X2>2005)

mean = c(mean(p1$value),mean(p2$value))
se = c((sd(p1$value)/sqrt(length(p1$value))),sd(p2$value)/sqrt(length(p2$value)))

# sample sizes
# per replicate 
# sexual
longT<-melt(T)
longT<-na.omit(longT)
p1<-subset(longT,longT$X2<2012)
mp1<-mean(p1$value)
p1_se<-sd(p1$value)/sqrt(length(p1$value))
length(p1$value)

p2<-subset(longT,longT$X2>2005)
mp2<-mean(p2$value)
p2_se<-sd(p2$value)/sqrt(length(p2$value))
length(p2$value)

# 4 sites
df<-subset(df,df$site!="jms" & df$site!="swend")

T<-t(tapply(df$mic,list(df$year,df$site),length))
I<-t(tapply(df$mic,list(df$year,df$site),sum))
F<-I/T
longF<-melt(F,id=c("year","frequency"))

longF<-na.omit(longF)

p1<-subset(longF,longF$X2<2012)
p2<-subset(longF,longF$X2>2005)

mean = c(mean(p1$value),mean(p2$value))
se = c((sd(p1$value)/sqrt(length(p1$value))),sd(p2$value)/sqrt(length(p2$value)))
