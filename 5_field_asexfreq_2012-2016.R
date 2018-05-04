# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Field: Asexual females increase in frequency
# Change in the frequency of asexuals from 2012 to 2016

#####
# Field, females, all sites, 2012-2016
rm(list=ls())
library(aod)
df<-read.csv("Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

df2=df
df2$ploidy[df2$ploidy==2]<-0
df2$ploidy[df2$ploidy==3]<-1

lrfit <- glm(ploidy~site*factor(year), data=df2, family=binomial)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/1693

lrfit2 <- glm(ploidy~site+factor(year), data=df2,family=binomial)
anova(lrfit2,lrfit,test="LRT")

summary(lrfit)
anova(lrfit,test="Chisq")

# odds ratios
exp(coef(lrfit))
exp(cbind(OR = coef(lrfit), confint(lrfit)))

# pseudo-R2
lr <- glm( ploidy~ 
             +  1, df=df2, family = binomial)
1-logLik(lrfit)/logLik(lr)

# reduced model
summary(lrfit2)
anova(lrfit2,test="Chisq")

# odds ratios
exp(coef(lrfit2))
exp(cbind(OR = coef(lrfit2), confint(lrfit2)))

# yearly changes
# 2016 against 
# 2015
l <- cbind(0,0,0,0,0,0,0,0,1,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

# 2015 against 
# 2014
l <- cbind(0,0,0,0,0,0,0,1,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

# 2014 against 
# 2013
l <- cbind(0,0,0,0,0,0,1,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#2013 against
# 2012
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), Terms = 7)

detach(df)

#####
# Repeat with four shared sites
# Field, females, four shared sites, 2012-2016
rm(list=ls())
df<-read.csv("Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0 & df$site!="swend" & df$site!="jms")
summary(df)
attach(df)

df2=df
df2$ploidy[df2$ploidy==2]<-0
df2$ploidy[df2$ploidy==3]<-1

lrfit <- glm(ploidy~site*factor(year), data=df2, family=binomial)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/1693

lrfit2 <- glm(ploidy~site+factor(year), data=df2,family=binomial)
anova(lrfit2,lrfit,test="LRT")

summary(lrfit)
anova(lrfit,test="Chisq")

# odds ratios
exp(coef(lrfit))
exp(cbind(OR = coef(lrfit), confint(lrfit)))

# pseudo-R2
lr <- glm( ploidy~ 
             +  1, df=df2, family = binomial)
1-logLik(lrfit)/logLik(lr)

# reduced model
summary(lrfit2)
anova(lrfit2,test="Chisq")

# odds ratios
exp(coef(lrfit2))
exp(cbind(OR = coef(lrfit2), confint(lrfit2)))

# 2016 against 
# 2015
l <- cbind(0,0,0,0,0,0,1,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 
# 0.018 # only really an increase in 2016...

# 2015 against 
# 2014
l <- cbind(0,0,0,0,0,1,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 
# same - 0.91

# 2014 against 
# 2013
l <- cbind(0,0,0,0,1,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 
# more 0.81

#2013 against
# 2012
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), Terms = 5)
# 0.33

detach(df)

#####
# Repeat with males
# Field, males+females, all sites, 2012-2016

rm(list=ls())
df<-read.csv("Field-2012-2016_withmales.csv")
df<-subset(df,df$ploidy>0)
summary(df)
attach(df)

df2=df
df2$ploidy[df2$ploidy==2]<-0
df2$ploidy[df2$ploidy==3]<-1

lrfit <- glm(ploidy~site*factor(year), data=df2, family=binomial)
par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/2127

lrfit2 <- glm(ploidy~site+factor(year), data=df2,family=binomial)
anova(lrfit2,lrfit,test="LRT")

# pseudo-R2
lr <- glm( ploidy~ 
             +  1, df=df2, family = binomial)
1-logLik(lrfit)/logLik(lr)

summary(lrfit)
anova(lrfit,test="Chisq")
summary(lrfit2)
anova(lrfit2,test="Chisq")

# odds ratios
exp(coef(lrfit2))
exp(cbind(OR = coef(lrfit2), confint(lrfit2)))

# yearly changes
# 2016 against 
# 2015
l <- cbind(0,0,0,0,0,0,0,0,1,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

# 2015 against 
# 2014
l <- cbind(0,0,0,0,0,0,0,1,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

# 2014 against 
# 2013
l <- cbind(0,0,0,0,0,0,1,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#2013 against
# 2012
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), Terms = 7)

detach(df)
#####
# Repeat with males and four shared sites
# Field, males+females, four shared sites, 2012-2016

rm(list=ls())
df<-read.csv("Field-2012-2016_withmales.csv")
df<-subset(dfdf$ploidy>0 & df$site!="swend" & data$site!="jms")
summary(df)
attach(df)

df2=df
df2$ploidy[df2$ploidy==2]<-0
df2$ploidy[df2$ploidy==3]<-1

lrfit <- glm(ploidy~site*factor(year), data=df2, family=binomial)
par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/2127

lrfit2 <- glm(ploidy~site+factor(year), data=df2,family=binomial)
anova(lrfit2,lrfit,test="LRT")

# pseudo-R2
lr <- glm( ploidy~ 
             +  1, df=df2, family = binomial)
1-logLik(lrfit)/logLik(lr)

summary(lrfit)
anova(lrfit,test="Chisq")
summary(lrfit2)
anova(lrfit2,test="Chisq")

# odds ratios
exp(coef(lrfit2))
exp(cbind(OR = coef(lrfit2), confint(lrfit2)))

# 2016 against 
# 2015
l <- cbind(0,0,0,0,0,0,1,-1)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l)

# 2015 against 
# 2014
l <- cbind(0,0,0,0,0,1,-1,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

# 2014 against 
# 2013
l <- cbind(0,0,0,0,1,-1,0,0)
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), L = l) 

#2013 against
# 2012
wald.test(b = coef(lrfit2), Sigma = vcov(lrfit2), Terms = 5) 

detach(df)
############################################
# Calculations

########################################
# females, asexual vs. sexual
# mean frequency asexuals
# mean of each year, across sites
# matches Fig. 2

rm(list=ls())
df<-read.csv("Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
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
mean = c(mean(longF$value[1:6]),mean(longF$value[7:12]),mean(longF$value[13:18]),
         mean(longF$value[19:24]),mean(longF$value[25:30]))
se = c(sd(longF$value[1:6])/sqrt(6),sd(longF$value[7:12])/sqrt(6),sd(longF$value[13:18])/sqrt(6),
       sd(longF$value[19:24])/sqrt(6),sd(longF$value[25:30])/sqrt(6))

# sample sizes
# per replicate 
longT<-melt(T)
s = mean(longT$value)
s_se = sd(longT$value)/sqrt(length(longT$value))

#####################################
# use mean to calculate cost of sex for each year
qt = mean[1:4]
qt1 = mean[2:5]

c = (qt1*(1-qt))/(qt*(1-qt1))

detach(df)
###################restrict to 4 sites
rm(list=ls())
df<-read.csv("Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0 & df$site!="swend" & df$site!="jms")
df$site <- factor(df$site)
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
mean = c(mean(longF$value[1:4]),mean(longF$value[5:8]),mean(longF$value[9:12]),
         mean(longF$value[13:16]),mean(longF$value[17:20]))
se = c(sd(longF$value[1:4])/sqrt(4),sd(longF$value[5:8])/sqrt(4),sd(longF$value[9:12])/sqrt(4),
       sd(longF$value[13:16])/sqrt(4),sd(longF$value[17:20])/sqrt(4))


# sample sizes
# per replicate 
longT<-melt(T)
s = mean(longT$value)
s_se = sd(longT$value)/sqrt(length(longT$value))

detach(df)


########################################
# probability of sampling an asexual female, by site

rm(list=ls())
df<-read.csv("Field-2012-2016.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

df$ploidy[df$ploidy==2]<-0
df$ploidy[df$ploidy==3]<-1

# sexuals and asexuals by site and year
T<-t(tapply(df$ploidy,list(df$year,df$site),length))
A<-t(tapply(df$ploidy,list(df$year,df$site),sum))
F<-A/T
longF<-melt(F,id=c("site","year","frequency"))
names(longF)=c("site","year","frequency")

# summarize by site
means<-tapply(longF$frequency,list(longF$site),mean)
se<-tapply(longF$frequency,list(longF$site),sd)/sqrt(5)

#sample sizes
longT<-melt(T,id=c("site","year","frequency"))
names(longT)=c("site","year","frequency")

mean(longT$frequency)
sd(longT$frequency)/sqrt(length(longT$frequency))

detach(df)
########################################
# probability of mic infection by site
rm(list=ls())
df<-read.csv("Field-2012-2016.csv") 
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

# infected and uninfected by year/site
# asexual
T<-t(tapply(df$mic,list(df$year,df$site),length))
I<-t(tapply(df$mic,list(df$year,df$site),sum))
F<-I/T
longF<-melt(F,id=c("year","frequency"))
names(longF)=c("site","year","frequency")

# summarize by site
means<-tapply(longF$frequency,list(longF$site),mean)
se<-tapply(longF$frequency,list(longF$site),sd)/sqrt(5)

#sample sizes
longT<-melt(T,id=c("site","year","frequency"))
names(longT)=c("site","year","frequency")

mean(longT$frequency)
sd(longT$frequency)/sqrt(length(longT$frequency))

detach(df)

########################################
# probability of mic infection by site
# for asexuals only
rm(list=ls())
df<-read.csv("Field-2012-2016.csv") 
df<-subset(df,df$sex==0 & df$ploidy==3)
summary(df)
attach(df)

# infected and uninfected by year/site
# asexual
T<-t(tapply(df$mic,list(df$year,df$site),length))
I<-t(tapply(df$mic,list(df$year,df$site),sum))
F<-I/T
longF<-melt(F,id=c("year","frequency"))
names(longF)=c("site","year","frequency")

# summarize by site
# infected and uninfected by year/site
# asexual
asex<-subset(df,df$ploidy==3)
aT<-t(tapply(asex$mic,list(asex$year,asex$site),length))
aI<-t(tapply(asex$mic,list(asex$year,asex$site),sum))
aF<-aI/aT
longaF<-melt(aF,id=c("year","frequency"))
longaF<-na.omit(longaF)

#means
m=mean(longaF$value)
se = sd(longaF$value)/sqrt(length(longaF$value))

subaF<-subset(longaF,longaF$X1!="swamp" & longaF$X1!="jms")
mean(subaF$value)
sd(subaF$value)/sqrt(length(subaF$value))
max(subaF$value)

# sample sizes
# per replicate 

#asexual
longaT<-melt(aT)
longaT<-na.omit(longaT)
a<-mean(longaT$value)
a_se<-sd(longaT$value)/sqrt(length(longaT$value))
