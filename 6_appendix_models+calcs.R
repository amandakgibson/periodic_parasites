# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Appendix models and calculations

#####
# Field: male v diploid female infection prevalence
rm(list=ls())
library(aod)
df<-read.csv("Field-2012-2016_withmales.csv")
df<-subset(df,df$ploidy==2)
summary(df)
attach(df)

lrfit <- glm(mic~site*factor(year)*factor(sex) + length, family=binomial)

par(mfrow=c(2,2))
plot(lrfit)
# overdispersion test
sumsquare = sum(residuals(lrfit, type = "pearson")^2)
sumsquare/1712

lrfit2 <- glm(mic~site*factor(year)+site*factor(sex)+factor(year)*factor(sex) + length, family=binomial)
anova(lrfit2,lrfit,test="LRT")
lrfit3 <- glm(mic~site*factor(year)+factor(year)*factor(sex) + length, family=binomial)
anova(lrfit3,lrfit2,test="LRT")

summary(lrfit3)
anova(lrfit3,test="Chisq")

# odds ratios
exp(coef(lrfit3))
exp(cbind(OR = coef(lrfit3), confint(lrfit3)))

lr <- glm( mic~ 
             +  1, family = binomial)
1-logLik(lrfit3)/logLik(lr)

detach(df)

########################################
# Field calculations

########################################
#  Microphallus in males and females
rm(list=ls())
df<-read.csv("Field-2012-2016_withmales.csv")
df<-subset(df,df$ploidy==2)
summary(df)
attach(df)

year<-factor(df$year)
ploidy<-factor(df$factor)
mic<-factor(df$mic)
sex<-factor(df$sex)

# infected and uninfected by year/site
# female
Fem<-subset(df,df$sex==0)
FT<-t(tapply(Fem$mic,list(Fem$year,Fem$site),length))
FI<-t(tapply(Fem$mic,list(Fem$year,Fem$site),sum))
FF<-FI/FT
longFF<-melt(FF,id=c("year","frequency"))

# infected and uninfected by year/site
# male
Male<-subset(df,df$sex==1)
MT<-t(tapply(Male$mic,list(Male$year,Male$site),length))
MI<-t(tapply(Male$mic,list(Male$year,Male$site),sum))
MF<-MI/MT
longMF<-melt(MF,id=c("year","frequency"))

longFF<-na.omit(longFF)
longMF<-na.omit(longMF)

#means
mean = c(mean(longFF$value),mean(longMF$value))
se = c((sd(longFF$value)/sqrt(length(longFF$value))),
       sd(longMF$value)/sqrt(length(longMF$value)))

# sample sizes
# per replicate 
# female
longFT<-melt(FT)
longFT<-na.omit(longFT)
F<-mean(longFT$value)
F_se<-sd(longFT$value)/sqrt(length(longFT$value))

#male
longMT<-melt(MT)
longMT<-na.omit(longMT)
M<-mean(longMT$value)
M_se<-sd(longMT$value)/sqrt(length(longMT$value))

detach(df)

#############################################
# Experiment: male v diploid female infection rate, Microphallus
rm(list=ls())
df<-read.csv("Bins_adults_withmales.csv")
df<-subset(df,df$ploidy==2)
summary(df)
attach(df)

# subset by sex and treatment
fC=subset(df,df$treatment==0 & df$sex==0)
fE=subset(df,df$treatment==1 & df$sex==0)
mC=subset(df,df$treatment==0 & df$sex==1)
mE=subset(df,df$treatment==1 & df$sex==1)

# control females
fCT=tapply(fC$mic,list(fC$year,fC$replicate),length)
LfCT=melt(fCT)
names(LfCT)=c("year","replicate","total")

fCM=tapply(fC$mic,list(fC$year,fC$replicate),sum)
LfCM=melt(fCM)
names(LfCM)=c("year","replicate","mic")

Cf<-cbind(LfCT,LfCM$mic)
names(Cf)=c("year","replicate","total","mic")
healthy=Cf$total-Cf$mic
sex=rep(0,length(Cf$year))
treatment=rep("control",length(Cf$year))
tank=rep(1:24)
Cf=cbind(Cf,healthy,sex,treatment,tank)

# exposed females
fET=tapply(fE$mic,list(fE$year,fE$replicate),length)
LfET=melt(fET)
names(LfET)=c("year","replicate","total")

fEM=tapply(fE$mic,list(fE$year,fE$replicate),sum)
LfEM=melt(fEM)
names(LfEM)=c("year","replicate","mic")

Ef<-cbind(LfET,LfEM$mic)
names(Ef)=c("year","replicate","total","mic")
healthy=Ef$total-Ef$mic
sex=rep(0,length(Ef$year))
treatment=rep("exposed",length(Ef$year))
tank=rep(25:48)
Ef=cbind(Ef,healthy,sex,treatment,tank)

# control males
mCT=tapply(mC$mic,list(mC$year,mC$replicate),length)
LmCT=melt(mCT)
names(LmCT)=c("year","replicate","total")

mCM=tapply(mC$mic,list(mC$year,mC$replicate),sum)
LmCM=melt(mCM)
names(LmCM)=c("year","replicate","mic")

Cm<-cbind(LmCT,LmCM$mic)
names(Cm)=c("year","replicate","total","mic")
healthy=Cm$total-Cm$mic
sex=rep(1,length(Cm$year))
treatment=rep("control",length(Cm$year))
tank=rep(1:24)
Cm=cbind(Cm,healthy,sex,treatment,tank)

# exposed males
mET=tapply(mE$mic,list(mE$year,mE$replicate),length)
LmET=melt(mET)
names(LmET)=c("year","replicate","total")

mEM=tapply(mE$mic,list(mE$year,mE$replicate),sum)
LmEM=melt(mEM)
names(LmEM)=c("year","replicate","mic")

Em<-cbind(LmET,LmEM$mic)
names(Em)=c("year","replicate","total","mic")
healthy=Em$total-Em$mic
sex=rep(1,length(Em$year))
treatment=rep("exposed",length(Em$year))
tank=rep(25:48)
Em=cbind(Em,healthy,sex,treatment,tank)

inf<-rbind(Cf,Cm,Ef,Em)

library(geepack)
geefit1<-geeglm( cbind(mic,healthy) ~ treatment+factor(year)+factor(sex)+
                   treatment*factor(year)+treatment*factor(sex)+factor(year)*factor(sex)+
                   treatment*factor(year)*factor(sex), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
summary(geefit1)
anova(geefit1,test="Wald")

detach(df)
######
# Experiment: male v diploid female infection rate; all castrating parasites
rm(list=ls())
df<-read.csv("Bins_adults_withmales.csv")
df<-subset(data,data$ploidy==2)
summary(df)
attach(df)

# subset by sex and treatment
fC=subset(df,df$treatment==0 & df$sex==0)
fE=subset(df,df$treatment==1 & df$sex==0)
mC=subset(df,df$treatment==0 & df$sex==1)
mE=subset(df,df$treatment==1 & df$sex==1)

# control females
fCT=tapply(fC$castrated,list(fC$year,fC$replicate),length)
LfCT=melt(fCT)
names(LfCT)=c("year","replicate","total")

fCM=tapply(fC$castrated,list(fC$year,fC$replicate),sum)
LfCM=melt(fCM)
names(LfCM)=c("year","replicate","cas")

Cf<-cbind(LfCT,LfCM$cas)
names(Cf)=c("year","replicate","total","cas")
healthy=Cf$total-Cf$cas
sex=rep(0,length(Cf$year))
treatment=rep("control",length(Cf$year))
tank=rep(1:24)
Cf=cbind(Cf,healthy,sex,treatment,tank)

# exposed females
fET=tapply(fE$castrated,list(fE$year,fE$replicate),length)
LfET=melt(fET)
names(LfET)=c("year","replicate","total")

fEM=tapply(fE$castrated,list(fE$year,fE$replicate),sum)
LfEM=melt(fEM)
names(LfEM)=c("year","replicate","cas")

Ef<-cbind(LfET,LfEM$cas)
names(Ef)=c("year","replicate","total","cas")
healthy=Ef$total-Ef$cas
sex=rep(0,length(Ef$year))
treatment=rep("exposed",length(Ef$year))
tank=rep(25:48)
Ef=cbind(Ef,healthy,sex,treatment,tank)

# control males
mCT=tapply(mC$castrated,list(mC$year,mC$replicate),length)
LmCT=melt(mCT)
names(LmCT)=c("year","replicate","total")

mCM=tapply(mC$castrated,list(mC$year,mC$replicate),sum)
LmCM=melt(mCM)
names(LmCM)=c("year","replicate","cas")

Cm<-cbind(LmCT,LmCM$cas)
names(Cm)=c("year","replicate","total","cas")
healthy=Cm$total-Cm$cas
sex=rep(1,length(Cm$year))
treatment=rep("control",length(Cm$year))
tank=rep(1:24)
Cm=cbind(Cm,healthy,sex,treatment,tank)

# exposed males
mET=tapply(mE$castrated,list(mE$year,mE$replicate),length)
LmET=melt(mET)
names(LmET)=c("year","replicate","total")

mEM=tapply(mE$cas,list(mE$year,mE$replicate),sum)
LmEM=melt(mEM)
names(LmEM)=c("year","replicate","cas")

Em<-cbind(LmET,LmEM$cas)
names(Em)=c("year","replicate","total","cas")
healthy=Em$total-Em$cas
sex=rep(1,length(Em$year))
treatment=rep("exposed",length(Em$year))
tank=rep(25:48)
Em=cbind(Em,healthy,sex,treatment,tank)

inf<-rbind(Cf,Cm,Ef,Em)

geefit1<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(sex)+
                   treatment*factor(year)+treatment*factor(sex)+factor(year)*factor(sex)+
                   treatment*factor(year)*factor(sex), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
summary(geefit1)
anova(geefit1,test="Wald")
wald.test(b = coef(geefit1), Sigma = vcov(geefit1), Terms = 14:16) 

geefit2<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(sex)+
                   treatment*factor(year)+treatment*factor(sex)+factor(year)*factor(sex), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
anova(geefit2,test="Wald")
summary(geefit2)
wald.test(b = coef(geefit2), Sigma = vcov(geefit2), Terms = 10)

geefit3<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(sex)+
                   treatment*factor(year)+factor(year)*factor(sex),
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
anova(geefit3,test="LRT")

detach(df)

######
# Experiment: frequency of males in control v exposed mesocosms
rm(list=ls())
df<-read.csv("Bins_adults_withmales.csv")
df<-subset(df,df$ploidy==2)
summary(df)
attach(df)

# subset by treatment
C=subset(df,df$treatment==0)
E=subset(df,df$treatment==1)

# control tanks
CT=tapply(C$sex,list(C$year,C$replicate),length)
LCT=melt(CT)
names(LCT)=c("year","replicate","total")

CM=tapply(C$sex,list(C$year,C$replicate),sum)
LCM=melt(CM)
names(LCM)=c("year","replicate","male")

control<-cbind(LCT,LCM$male)
names(control)=c("year","replicate","total","male")
female=control$total-control$male
treatment=rep("control",length(control$year))
control=cbind(control,female,treatment)

# exposed tanks
ET=tapply(E$sex,list(E$year,E$replicate),length)
LET=melt(ET)
names(LET)=c("year","replicate","total")

EM=tapply(E$sex,list(E$year,E$replicate),sum)
LEM=melt(EM)
names(LEM)=c("year","replicate","male")

exposed<-cbind(LET,LEM$male)
names(exposed)=c("year","replicate","total","male")
female=exposed$total-exposed$male
treatment=rep("exposed",length(exposed$year))
exposed=cbind(exposed,female,treatment)

# frequencies
cfreq=control$male/control$total
efreq=exposed$male/exposed$total

mean(cfreq)
sd(cfreq)/sqrt(24)

mean(efreq)
sd(efreq)/sqrt(24)

# paired t-test
t.test(cfreq,efreq,paired=TRUE) 

library(nortest)
ad.test(cfreq) 
ad.test(efreq) 

# non-parametric
wilcox.test(cfreq,efreq,paired=TRUE)

detach(df)
######
# Experiment: frequency of males in control v exposed mesocosms
rm(list=ls())
df<-read.csv("Bins_adults.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

# subset by ploidy and treatment
aC=subset(df,df$treatment==0 & df$ploidy==3)
aE=subset(df,df$treatment==1 & df$ploidy==3)
sC=subset(df,df$treatment==0 & df$ploidy==2)
sE=subset(df,df$treatment==1 & df$ploidy==2)

# control asexuals
aCT=tapply(aC$ploidy,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","triploid")

# exposed asexuals
aET=tapply(aE$ploidy,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","triploid")

# control sexuals
sCT=tapply(sC$ploidy,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","diploid")

# exposed sexuals
sET=tapply(sE$ploidy,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","diploid")

LCT=cbind(LaCT,LsCT$diploid)
treatment=rep("control",length(LCT$triploid))
LCT=cbind(LCT,treatment)
names(LCT)=c("year","replicate","triploid","diploid","treatment")

LET=cbind(LaET,LsET$diploid)
treatment=rep("exposed",length(LET$triploid))
LET=cbind(LET,treatment)
names(LET)=c("year","replicate","triploid","diploid","treatment")

cfreq3=LCT$triploid/(LCT$triploid+LCT$diploid)
efreq3=LET$triploid/(LET$triploid+LET$diploid)

mean(cfreq3)
mean(efreq3)
mean(cfreq3-efreq3)

t.test(cfreq3,efreq3,paired=T) 
wilcox.test(cfreq3,efreq3,paired=TRUE) 

ad.test(cfreq3)
ad.test(efreq3)

detach(df)
################################
# Experiment Calculations
rm(list=ls())
df<-read.csv("Bins_adults_withmales.csv")
df<-subset(df,df$ploidy==2)
summary(df)
attach(df)

###################################################
# males vs females in control tanks, microphallus

Cf<-subset(df,df$treatment==0 &df$sex==0)
tCf=tapply(Cf$mic,list(Cf$replicate,Cf$year),length)
mCf=tapply(Cf$mic,list(Cf$replicate,Cf$year),sum)
fCf=mCf/tCf

Cm<-subset(df,df$treatment==0 & df$sex==1)
tCm=tapply(Cm$mic,list(Cm$replicate,Cm$year),length)
mCm=tapply(Cm$mic,list(Cm$replicate,Cm$year),sum)
fCm=mCm/tCm

mean=c(mean(fCf),mean(fCm))
se=c(sd(fCf)/sqrt(24),sd(fCm)/sqrt(24))

# sample sizes
mean=c(mean(tCf),mean(tCm))
se=c(sd(tCf)/sqrt(24),sd(tCm)/sqrt(24))

###################################################
# males vs females in exposed tanks, microphallus

Ef<-subset(df,df$treatment==1 &df$sex==0)
tEf=tapply(Ef$mic,list(Ef$replicate,Ef$year),length)
mEf=tapply(Ef$mic,list(Ef$replicate,Ef$year),sum)
fEf=mEf/tEf

Em<-subset(df,df$treatment==1 & df$sex==1)
tEm=tapply(Em$mic,list(Em$replicate,Em$year),length)
mEm=tapply(Em$mic,list(Em$replicate,Em$year),sum)
fEm=mEm/tEm

mean=c(mean(fEf),mean(fEm))
se=c(sd(fEf)/sqrt(24),sd(fEm)/sqrt(24))

# sample sizes
mean=c(mean(tEf),mean(tEm))
se=c(sd(tEf)/sqrt(24),sd(tEm)/sqrt(24))

###################################################
# males vs females in control tanks, all parasites

Cf<-subset(df,df$treatment==0 &df$sex==0)
tCf=tapply(Cf$castrated,list(Cf$replicate,Cf$year),length)
mCf=tapply(Cf$castrated,list(Cf$replicate,Cf$year),sum)
fCf=mCf/tCf

Cm<-subset(df,df$treatment==0 & df$sex==1)
tCm=tapply(Cm$castrated,list(Cm$replicate,Cm$year),length)
mCm=tapply(Cm$castrated,list(Cm$replicate,Cm$year),sum)
fCm=mCm/tCm

mean=c(mean(fCf),mean(fCm))
se=c(sd(fCf)/sqrt(24),sd(fCm)/sqrt(24))

# sample sizes
mean=c(mean(tCf),mean(tCm))
se=c(sd(tCf)/sqrt(24),sd(tCm)/sqrt(24))

###################################################
# males vs females in exposed tanks, all parasites

Ef<-subset(df,df$treatment==1 &df$sex==0)
tEf=tapply(Ef$castrated,list(Ef$replicate,Ef$year),length)
mEf=tapply(Ef$castrated,list(Ef$replicate,Ef$year),sum)
fEf=mEf/tEf

Em<-subset(df,df$treatment==1 & df$sex==1)
tEm=tapply(Em$castrated,list(Em$replicate,Em$year),length)
mEm=tapply(Em$castrated,list(Em$replicate,Em$year),sum)
fEm=mEm/tEm

mean=c(mean(fEf),mean(fEm))
se=c(sd(fEf)/sqrt(24),sd(fEm)/sqrt(24))

# sample sizes
mean=c(mean(tEf),mean(tEm))
se=c(sd(tEf)/sqrt(24),sd(tEm)/sqrt(24))

detach(df)
#####################################################
# sex ratio in control v exposed tanks
rm(list=ls())
df<-read.csv("Bins_adults_withmales.csv")
df<-subset(df,df$ploidy==2)
summary(df)
attach(df)

# subset by treatment
C=subset(df,df$treatment==0)
E=subset(df,df$treatment==1)

# now divide each subset by tank and year
# control tanks
CT=tapply(C$sex,list(C$year,C$replicate),length)
LCT=melt(CT)
names(LCT)=c("year","replicate","total")

CM=tapply(C$sex,list(C$year,C$replicate),sum)
LCM=melt(CM)
names(LCM)=c("year","replicate","male")

control<-cbind(LCT,LCM$male)
names(control)=c("year","replicate","total","male")
female=control$total-control$male
treatment=rep("control",length(control$year))
control=cbind(control,female,treatment)

# exposed tanks
ET=tapply(E$sex,list(E$year,E$replicate),length)
LET=melt(ET)
names(LET)=c("year","replicate","total")

EM=tapply(E$sex,list(E$year,E$replicate),sum)
LEM=melt(EM)
names(LEM)=c("year","replicate","male")

exposed<-cbind(LET,LEM$male)
names(exposed)=c("year","replicate","total","male")
female=exposed$total-exposed$male
treatment=rep("exposed",length(exposed$year))
exposed=cbind(exposed,female,treatment)

cfreq=control$male/control$total
efreq=exposed$male/exposed$total

mean=c(mean(cfreq),mean(efreq))
se=c(sd(cfreq)/sqrt(length(cfreq)),sd(efreq/sqrt(length(efreq))))

# sample sizes
mean=c(mean(control$total),mean(exposed$total))
se=c(sd(control$total)/sqrt(length(control$total)),sd(exposed$total)/sqrt(length(exposed$total)))

detach(df)
#############################################################
# ratio diploids to triploids in control v exposed tanks
rm(list=ls())
df<-read.csv("Bins_adults.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

# subset by ploidy and treatment
aC=subset(df,df$treatment==0 & df$ploidy==3)
aE=subset(df,df$treatment==1 & df$ploidy==3)
sC=subset(df,df$treatment==0 & df$ploidy==2)
sE=subset(df,df$treatment==1 & df$ploidy==2)

# get lengths by treatment and replicate
# control asexuals
aCT=tapply(aC$ploidy,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","triploid")

# exposed asexuals
aET=tapply(aE$ploidy,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","triploid")

# control sexuals
sCT=tapply(sC$ploidy,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","diploid")

# exposed sexuals
sET=tapply(sE$ploidy,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","diploid")

LCT=cbind(LaCT,LsCT$diploid)
treatment=rep("control",length(LCT$triploid))
LCT=cbind(LCT,treatment)
names(LCT)=c("year","replicate","triploid","diploid","treatment")

LET=cbind(LaET,LsET$diploid)
treatment=rep("exposed",length(LET$triploid))
LET=cbind(LET,treatment)
names(LET)=c("year","replicate","triploid","diploid","treatment")

cfreq3=LCT$triploid/(LCT$triploid+LCT$diploid)
efreq3=LET$triploid/(LET$triploid+LET$diploid)

mean=c(mean(cfreq3),mean(efreq3))
se=c(sd(cfreq3)/sqrt(length(cfreq3)),sd(efreq3/sqrt(length(efreq3))))

# sample sizes
mean(LCT$triploid+LCT$diploid)
mean(LET$triploid+LET$diploid)

sd(LCT$triploid+LCT$diploid)/sqrt(length(LCT$triploid))
sd(LET$triploid+LET$diploid)/sqrt(length(LET$triploid))

