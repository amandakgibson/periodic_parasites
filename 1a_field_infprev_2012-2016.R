# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Temporal shifts in relative infection prevalence, part A
# Comparison of Microphallus infection rates in the field in 2012-2016

######
# Field, 2012-2016, all sites, females
rm(list=ls())
library(aod)
library(lmtest)

df<-read.csv("Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)

attach(df)

##########
# year and site and ploidy as predictors of sex
# covariate length

lrfit <- glm(mic~site*factor(ploidy) + factor(year)*factor(ploidy) + length, family=binomial)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/1702

lrfit2 <- glm(mic~site*factor(ploidy) + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit2,lrfit,test="LRT")

lrfit3 <- glm(mic~site + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit3,lrfit2,test="LRT")

# results
summary(lrfit2)
anova(lrfit2,test="Chisq")

# pseudo-R2
lr <- glm( mic ~ 
             +  1, family = binomial)
1-logLik(lrfit2)/logLik(lr)

# site*ploidy interaction
# halfway against jms
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#half-swamp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#half-swend
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#half-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#jms-swamp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#jms-swend
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#jms-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#swamp-swend
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#swamp-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#swend-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 


l <- cbind(0,0,0,0,0,0,1,0,0,0,0,0,-1,0,0,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

detach(df)

######
# Repeat: four shared sites only
# Field, 2012-2016, four sites, females
df<-subset(df,df$sex==0 & df$ploidy>0 & site!="jms" & site!="swend")
summary(df)
attach(df)

lrfit <- glm(mic~site*factor(ploidy) + factor(year)*factor(ploidy) + length, family=binomial)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/1119

# do we need the interaction term?
lrfit2 <- glm(mic~site*factor(ploidy) + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit2,lrfit,test="LRT")

lrfit3 <- glm(mic~site + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit3,lrfit2,test="LRT") 

# results
summary(lrfit2)
anova(lrfit2,test="Chisq")

# evaluate pseudo-R2
#intercept only model
lr <- glm( mic ~ 
             +  1, family = binomial)
1-logLik(lrfit2)/logLik(lr)

detach(df)

######
# Repeat: with males included
# Field, 2012-2016, all sites, females+males
rm(list=ls())

df<-read.csv("Field-2012-2016_withmales.csv")
df<-subset(df,df$ploidy>0)
summary(df)
attach(df)

lrfit <- glm(mic~site*factor(ploidy) + factor(year)*factor(ploidy) + length, family=binomial)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/2136

lrfit2 <- glm(mic~site*factor(ploidy) + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit2,lrfit,test="LRT")

lrfit3 <- glm(mic~site + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit3,lrfit2,test="LRT")

summary(lrfit2)
anova(lrfit2,test="Chisq")

# pseudo-R2
lr <- glm( mic ~ 
             +  1, family = binomial)
1-logLik(lrfit2)/logLik(lr)

# ####
# site*ploidy interaction

l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#half-swamp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#half-swend
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#half-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#jms-swamp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#jms-swend
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#jms-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#swamp-swend
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l)

#swamp-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#swend-wp
l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

detach(df)

######
# Repeat: four shared sites, with males
# Field, 2012-2016, four sites, females+males
rm(list=ls())
df<-read.csv("Field-2012-2016_withmales.csv")
df<-subset(df,df$ploidy>0 & site!="jms" & site!="swend")
summary(df)
attach(df)

lrfit <- glm(mic~site*factor(ploidy) + factor(year)*factor(ploidy) + length, family=binomial)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/1470

lrfit2 <- glm(mic~site*factor(ploidy) + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit2,lrfit,test="LRT")

lrfit3 <- glm(mic~site + factor(year)+factor(ploidy) + length, family=binomial)
anova(lrfit3,lrfit2,test="LRT")

summary(lrfit2)
anova(lrfit2,test="Chisq")

lr <- glm( mic ~ 
             +  1, family = binomial)
1-logLik(lrfit2)/logLik(lr)

detach(df)

#####################################
# Calculations

# mean MICROPHALLUS infection prevalence
# mean over each site*year 
# matches boxplot insets of Figure 1
rm(list=ls())
df<-read.csv("../Field-2012-2016.csv") 
df<-subset(df,df$sex==0 & df$ploidy>0)
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
se = c((sd(longsF$value)/sqrt(length(longsF$value))),
       sd(longaF$value)/sqrt(length(longaF$value)))

subaF<-subset(longaF,longaF$X1!="swamp" & longaF$X1!="jms")
mean(subaF$value)
sd(subaF$value)/sqrt(length(subaF$value))
max(subaF$value)

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

# fold difference by site
names(longsF)=c("site","year","mic")
meanS_site<-tapply(longsF$mic,c(longsF$site),mean)

names(longaF)=c("site","year","mic")
meanA_site<-tapply(longaF$mic,c(longaF$site),mean)

fold=meanS_site/meanA_site
mean1=tapply(longsF$mic,list(longsF$site),mean)
mean2=tapply(longaF$mic,list(longaF$site),mean)
se1=tapply(longsF$mic,list(longsF$site),sd)/sqrt(5)
se2=tapply(longaF$mic,list(longaF$site),sd)/sqrt(5)

# se of this relationship
se_fold=fold*sqrt(((se1/mean1)^2)+((se2/mean2)^2))

detach(df)

########################################
# repeat for 4 restricted sites
rm(list=ls())
df<-read.csv("../Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0 & df$site!="jms" & df$site!="swend")
summary(df)
attach(df)

year<-factor(df$year)
ploidy<-factor(df$factor)
mic<-factor(df$mic)

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
se = c((sd(longsF$value)/sqrt(length(longsF$value))),
       sd(longaF$value)/sqrt(length(longaF$value)))

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

#########################################
# Site means
rm(list=ls())
df<-read.csv("../Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

year<-factor(df$year)
ploidy<-factor(df$factor)
mic<-factor(df$mic)

# infected and uninfected by year/site
T<-t(tapply(df$mic,list(df$year,df$site),length))
I<-t(tapply(df$mic,list(df$year,df$site),sum))
F<-I/T
longF<-melt(F,id=c("year","frequency"))
colnames(longF)=c("site","year","mic")

# by site
site_mean<-tapply(longF$mic,list(longF$site),mean)
site_sd<-tapply(longF$mic,list(longF$site),sd)
site_se<-site_sd/sqrt(5)

# sample size
longT<-melt(T)
longT<-na.omit(longT)
T_mean<-mean(longT$value)
T_se<-sd(longT$value)/sqrt(length(longT$value))

####################################
# Year means

T<-t(tapply(df$mic,list(df$year,df$site),length))
I<-t(tapply(df$mic,list(df$year,df$site),sum))
F<-I/T
longF<-melt(F,id=c("year","frequency"))
colnames(longF)=c("site","year","mic")

# by year
year_mean<-tapply(longF$mic,list(longF$year),mean)
year_sd<-tapply(longF$mic,list(longF$year),sd)
year_se<-year_sd/sqrt(6)

# sample size
longT<-melt(T)
longT<-na.omit(longT)
T_mean<-mean(longT$value)
T_se<-sd(longT$value)/sqrt(length(longT$value))

detach(df)

###############################################
##ambiguous DNA content######################
rm(list=ls())
df<-read.csv("../Field-2012-2016.csv") 
df<-subset(df,df$sex==0)
summary(df)

attach(df)

d<-subset(df,df$ploidy!="NA")
l<-tapply(d$ploidy,list(d$site,d$year),length)

d2<-subset(df,df$ploidy==0)
l2<-tapply(d2$ploidy,list(d2$site,d2$year),length)
l2[is.na(l2)] <- 0

mean(l2/l)
sd(l2/l)/sqrt(30)

detach(df)


