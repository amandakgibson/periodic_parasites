# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Asexual frequency declined between sampling periods
# Comparison of frequency of asexuals in the field in 2001-2005 vs. 2012-2016

######
# Field, all sites, females
rm(list=ls())
detach(df)
library(aod)
library(lme4)
library(geepack)
library(lmtest)
library(reshape)

df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & ploidy>0)
summary(df)
attach(df)

# triploids
dfA<-subset(df,df$ploidy==3)
matrixA<-tapply(dfA$ploidy,list(dfA$year,dfA$site),length)
longA<-melt(matrixA,id=c("year","site"))
colnames(longA)=c("year","site","triploid")

# diploids
dfS<-subset(df,df$ploidy==2)
matrixS<-tapply(dfS$ploidy,list(dfS$year,dfS$site),length)
longS<-melt(matrixS,id=c("year","site"))
colnames(longS)=c("year","site","diploid")

long<-cbind(longA,longS$diploid)
colnames(long)=c("year","site","triploid","diploid")

period<-rep(c(0,0,0,0,0,1,1,1,1,1),6)

long<-cbind(period,long)

#delete non-existent rows where we don't have data
long_final<-long[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]

long_final[is.na(long_final)] <- 0
df.manip<-data.frame(long_final)

detach(df)
attach(df.manip)

# now the model
lrfit <- glm(cbind(triploid,diploid)~ site+factor(period), family=binomial,data=df.manip)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/43 

#quasi
lrfitq1<-glm(cbind(triploid,diploid)~ site+factor(period), family=quasibinomial,data=df.manip)
lrfitq2<-glm(cbind(triploid,diploid)~ site, family=quasibinomial,df=df.manip)
anova(lrfitq1,lrfitq2,test="F")
lrfitq3<-glm(cbind(triploid,diploid)~ period, family=quasibinomial,df=df.manip)
anova(lrfitq1,lrfitq3,test="F")
lrfitq4<-glm(cbind(triploid,diploid)~+1,family=quasibinomial, df=df.manip)
anova(lrfitq1,lrfitq4,test="F")

anova(lrfitq1, test="F")

detach(df.manip)

######
# Repeat: four shared sites
# Field, four sites, females

rm(list=ls())
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df2<-subset(df2,df2$site!="swend" & df2$site!="jms")
df2$site <- factor(df2$site)

df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & ploidy>0)
summary(df)
attach(df)

# triploids
dfA<-subset(df,df$ploidy==3)
matrixA<-tapply(dfA$ploidy,list(dfA$year,dfA$site),length)
longA<-melt(matrixA,id=c("year","site"))
colnames(longA)=c("year","site","triploid")

# diploids
dfS<-subset(df,df$ploidy==2)
matrixS<-tapply(dfS$ploidy,list(dfS$year,dfS$site),length)
longS<-melt(matrixS,id=c("year","site"))
colnames(longS)=c("year","site","diploid")

long<-cbind(longA,longS$diploid)
colnames(long)=c("year","site","triploid","diploid")

period<-rep(c(0,0,0,0,0,1,1,1,1,1),4)

long<-cbind(period,long)

long[is.na(long)] <- 0
df.manip<-data.frame(long)

detach(df)
attach(df.manip)

# now the model
lrfit <- glm(cbind(triploid,diploid)~ site*factor(period), family=binomial,data=df.manip)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/32 

#quasi
lrfitq1<-glm(cbind(triploid,diploid)~ site*factor(period), family=quasibinomial,data=df.manip)
summary(lrfitq1)
plot(lrfitq1)
anova(lrfitq1,test="F")

detach(df.manip)

######
# # Repeat: males+females
# Field, all sites, females+males
rm(list=ls())
detach(df)
library(aod)
library(lme4)
library(geepack)
library(lmtest)
library(reshape)

df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016_withmales.csv")
df<-rbind(df1,df2)
df<-subset(df,ploidy>0)
summary(df)
attach(df)

# triploids
dfA<-subset(df,df$ploidy==3)
matrixA<-tapply(dfA$ploidy,list(dfA$year,dfA$site),length)
longA<-melt(matrixA,id=c("year","site"))
colnames(longA)=c("year","site","triploid")

# diploids
dfS<-subset(df,df$ploidy==2)
matrixS<-tapply(dfS$ploidy,list(dfS$year,dfS$site),length)
longS<-melt(matrixS,id=c("year","site"))
colnames(longS)=c("year","site","diploid")

long<-cbind(longA,longS$diploid)
colnames(long)=c("year","site","triploid","diploid")
period<-rep(c(0,0,0,0,0,1,1,1,1,1),6)
long<-cbind(period,long)

#delete non-existent rows where we don't have data
long_final<-long[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-41,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]
long_final<-long_final[-46,]

long_final[is.na(long_final)] <- 0
df.manip<-data.frame(long_final)

detach(df)
attach(df.manip)

lrfit <- glm(cbind(triploid,diploid)~ site+factor(period), family=binomial,data=df.manip)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/43 

#quasi
lrfitq1<-glm(cbind(triploid,diploid)~ site+factor(period), family=quasibinomial,data=df.manip)
summary(lrfitq1)
plot(lrfitq1)
anova(lrfitq1,test="F")

detach(df.manip)

######
# Repeat: four shared sites, males and females
# Field, for sites, females+males
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016_withmales.csv")
df2<-subset(df2,df2$site!="swend" & df2$site!="jms")
df2$site <- factor(df2$site)

df<-rbind(df1,df2)
df<-subset(df,ploidy>0)
summary(df)
attach(df)

# triploids
dfA<-subset(df,df$ploidy==3)
matrixA<-tapply(dfA$ploidy,list(dfA$year,dfA$site),length)
longA<-melt(matrixA,id=c("year","site"))
colnames(longA)=c("year","site","triploid")

# diploids
dfS<-subset(df,df$ploidy==2)
matrixS<-tapply(dfS$ploidy,list(dfS$year,dfS$site),length)
longS<-melt(matrixS,id=c("year","site"))
colnames(longS)=c("year","site","diploid")

long<-cbind(longA,longS$diploid)
colnames(long)=c("year","site","triploid","diploid")

period<-rep(c(0,0,0,0,0,1,1,1,1,1),4)
long<-cbind(period,long)
long[is.na(long)] <- 0
df.manip<-data.frame(long)

detach(df)
attach(df.manip)

lrfit <- glm(cbind(triploid,diploid)~ site+factor(period), family=binomial,data=df.manip)
par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/35 

#quasi
lrfitq1<-glm(cbind(triploid,diploid)~ site+factor(period), family=quasibinomial,data=df.manip)
summary(lrfitq1)
plot(lrfitq1)
anova(lrfitq1,test="F")

detach(df.manip)

############################
#Calculations

########################################
# females, asexual vs. sexual
# mean prevalence of asexuals
# mean over each site*year 
# matches Fig. 1e
rm(list=ls())
df1<-read.csv("Field-2001-2005.csv")
df2<-read.csv("Field-2012-2016.csv")
df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & ploidy>0)
summary(df)
attach(df)

df$ploidy[df$ploidy==2]<-0
df$ploidy[df$ploidy==3]<-1

# sexuals and asexuals by site and year
T<-t(tapply(df$ploidy,list(df$year,df$site),length))
A<-t(tapply(df$ploidy,list(df$year,df$site),sum))
F<-A/T
longF<-melt(F,id=c("site","year","frequency"))

# means
mean = c(mean(longF$value[1:30],na.rm=T),mean(longF$value[31:60]))
se = c((sd(longF$value[1:30],na.rm=T)/sqrt(20)),
       sd(longF$value[31:60])/sqrt(length(longF$value[31:60])))

# sample sizes
# per replicate 
longT<-melt(T)
s<-c(mean(longT$value[1:30],na.rm=T),mean(longT$value[31:60]))
s_se = c((sd(longT$value[1:30],na.rm=T)/sqrt(20)),
         sd(longT$value[31:60])/sqrt(length(longT$value[31:60])))

detach(df)
###############repeat for 4 restricted sites
rm(list=ls())
df1<-read.csv("../Field-2001-2005.csv")
df2<-read.csv("../Field-2012-2016.csv")
df2<-subset(df2,df2$site!="swend" & df2$site!="jms")
df2$site <- factor(df2$site)

df<-rbind(df1,df2)
df<-subset(df,df$sex==0 & ploidy>0)
summary(df)
attach(df)

df$ploidy[df$ploidy==2]<-0
df$ploidy[df$ploidy==3]<-1

# sexuals and asexuals by site and year
T<-t(tapply(df$ploidy,list(df$year,df$site),length))
A<-t(tapply(df$ploidy,list(df$year,df$site),sum))
F<-A/T
longF<-melt(F,id=c("site","year","frequency"))

# means
mean = c(mean(longF$value[1:20]),mean(longF$value[21:40]))
se = c((sd(longF$value[1:20])/sqrt(length(longF$value[1:20]))),
       sd(longF$value[21:40])/sqrt(length(longF$value[21:40])))

# sample sizes
# per replicate 
longT<-melt(T)
s<-c(mean(longT$value[1:20]),mean(longT$value[21:40]))
s_se = c((sd(longT$value[1:20])/sqrt(length(longT$value[1:20]))),
         sd(longT$value[21:40])/sqrt(length(longT$value[21:40])))

