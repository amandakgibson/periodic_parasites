# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Experimental test: coevolving parasites have reduced infectivity to asexual hosts
# Comparison of Microphallus infection rates following controlled exposures, 2012-2015

######
# Experiment, females, Microphallus
detach(df)

rm(list=ls())
detach(df)
library(aod)
library(reshape)
library(geepack)

df<-read.csv("Bins_adults.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

# subset by ploidy and treatment
aC=subset(df,df$treatment==0 & df$ploidy==3)
aE=subset(df,df$treatment==1 & df$ploidy==3)
sC=subset(df,df$treatment==0 & df$ploidy==2)
sE=subset(df,df$treatment==1 & df$ploidy==2)

# group by tank and year
aCT=tapply(aC$mic,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","total")

aCM=tapply(aC$mic,list(aC$year,aC$replicate),sum)
LaCM=melt(aCM)
names(LaCM)=c("year","replicate","mic")

Ca<-cbind(LaCT,LaCM$mic)
names(Ca)=c("year","replicate","total","mic")
healthy=Ca$total-Ca$mic
ploidy=rep(3,length(Ca$year))
treatment=rep("control",length(Ca$year))
tank=rep(1:24)
Ca=cbind(Ca,healthy,ploidy,treatment,tank)

# exposed asexuals
aET=tapply(aE$mic,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","total")

aEM=tapply(aE$mic,list(aE$year,aE$replicate),sum)
LaEM=melt(aEM)
names(LaEM)=c("year","replicate","mic")

Ea<-cbind(LaET,LaEM$mic)
names(Ea)=c("year","replicate","total","mic")
healthy=Ea$total-Ea$mic
ploidy=rep(3,length(Ea$year))
treatment=rep("exposed",length(Ea$year))
tank=rep(25:48)
Ea=cbind(Ea,healthy,ploidy,treatment,tank)

# control sexuals
sCT=tapply(sC$mic,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","total")

sCM=tapply(sC$mic,list(sC$year,sC$replicate),sum)
LsCM=melt(sCM)
names(LsCM)=c("year","replicate","mic")

Cs<-cbind(LsCT,LsCM$mic)
names(Cs)=c("year","replicate","total","mic")
healthy=Cs$total-Cs$mic
ploidy=rep(2,length(Cs$year))
treatment=rep("control",length(Cs$year))
tank=rep(1:24)
Cs=cbind(Cs,healthy,ploidy,treatment,tank)

# exposed sexuals
sET=tapply(sE$mic,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","total")

sEM=tapply(sE$mic,list(sE$year,sE$replicate),sum)
LsEM=melt(sEM)
names(LsEM)=c("year","replicate","mic")

Es<-cbind(LsET,LsEM$mic)
names(Es)=c("year","replicate","total","mic")
healthy=Es$total-Es$mic
ploidy=rep(2,length(Es$year))
treatment=rep("exposed",length(Es$year))
tank=rep(25:48)
Es=cbind(Es,healthy,ploidy,treatment,tank)

inf<-rbind(Ca,Cs,Ea,Es)

geefit1<-geeglm( cbind(mic,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy)+
                   treatment*factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
wald.test(b = coef(geefit1), Sigma = vcov(geefit1), Terms = 14:16)

geefit2<-geeglm( cbind(mic,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
wald.test(b = coef(geefit2), Sigma = vcov(geefit2), Terms = 10)


geefit3<-geeglm( cbind(mic,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")

anova(geefit3,test="LRT")

detach(df)

######
# Repeat: all parasites
# Experiment, females, all castrating trematodes

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
aCT=tapply(aC$castrated,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","total")

aCM=tapply(aC$castrated,list(aC$year,aC$replicate),sum)
LaCM=melt(aCM)
names(LaCM)=c("year","replicate","cas")

Ca<-cbind(LaCT,LaCM$cas)
names(Ca)=c("year","replicate","total","cas")
healthy=Ca$total-Ca$cas
ploidy=rep(3,length(Ca$year))
treatment=rep("control",length(Ca$year))
tank=rep(1:24)
Ca=cbind(Ca,healthy,ploidy,treatment,tank)

# exposed asexuals
aET=tapply(aE$castrated,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","total")

aEM=tapply(aE$castrated,list(aE$year,aE$replicate),sum)
LaEM=melt(aEM)
names(LaEM)=c("year","replicate","cas")

Ea<-cbind(LaET,LaEM$cas)
names(Ea)=c("year","replicate","total","cas")
healthy=Ea$total-Ea$cas
ploidy=rep(3,length(Ea$year))
treatment=rep("exposed",length(Ea$year))
tank=rep(25:48)
Ea=cbind(Ea,healthy,ploidy,treatment,tank)

# control sexuals
sCT=tapply(sC$castrated,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","total")

sCM=tapply(sC$castrated,list(sC$year,sC$replicate),sum)
LsCM=melt(sCM)
names(LsCM)=c("year","replicate","cas")

Cs<-cbind(LsCT,LsCM$cas)
names(Cs)=c("year","replicate","total","cas")
healthy=Cs$total-Cs$cas
ploidy=rep(2,length(Cs$year))
treatment=rep("control",length(Cs$year))
tank=rep(1:24)
Cs=cbind(Cs,healthy,ploidy,treatment,tank)

# exposed sexuals
sET=tapply(sE$castrated,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","total")

sEM=tapply(sE$castrated,list(sE$year,sE$replicate),sum)
LsEM=melt(sEM)
names(LsEM)=c("year","replicate","cas")

Es<-cbind(LsET,LsEM$cas)
names(Es)=c("year","replicate","total","cas")
healthy=Es$total-Es$cas
ploidy=rep(2,length(Es$year))
treatment=rep("exposed",length(Es$year))
tank=rep(25:48)
Es=cbind(Es,healthy,ploidy,treatment,tank)

inf<-rbind(Ca,Cs,Ea,Es)

geefit1<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy)+
                   treatment*factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
anova(geefit1,test="LRT")
wald.test(b = coef(geefit1), Sigma = vcov(geefit1), Terms = 14:16)

geefit2<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy), 
                 family = binomial,df=inf,
                 id=tank,corstr="exchangeable")
anova(geefit2,test="LRT")
summary(geefit2)

detach(df)

######
# Repeat: males and females
# Experiment, males+females, Microphallus
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

# control asexuals
aCT=tapply(aC$mic,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","total")

aCM=tapply(aC$mic,list(aC$year,aC$replicate),sum)
LaCM=melt(aCM)
names(LaCM)=c("year","replicate","mic")

Ca<-cbind(LaCT,LaCM$mic)
names(Ca)=c("year","replicate","total","mic")
healthy=Ca$total-Ca$mic
ploidy=rep(3,length(Ca$year))
treatment=rep("control",length(Ca$year))
tank=rep(1:24)
Ca=cbind(Ca,healthy,ploidy,treatment,tank)

# exposed asexuals
aET=tapply(aE$mic,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","total")

aEM=tapply(aE$mic,list(aE$year,aE$replicate),sum)
LaEM=melt(aEM)
names(LaEM)=c("year","replicate","mic")

Ea<-cbind(LaET,LaEM$mic)
names(Ea)=c("year","replicate","total","mic")
healthy=Ea$total-Ea$mic
ploidy=rep(3,length(Ea$year))
treatment=rep("exposed",length(Ea$year))
tank=rep(25:48)
Ea=cbind(Ea,healthy,ploidy,treatment,tank)

# control sexuals
sCT=tapply(sC$mic,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","total")

sCM=tapply(sC$mic,list(sC$year,sC$replicate),sum)
LsCM=melt(sCM)
names(LsCM)=c("year","replicate","mic")

Cs<-cbind(LsCT,LsCM$mic)
names(Cs)=c("year","replicate","total","mic")
healthy=Cs$total-Cs$mic
ploidy=rep(2,length(Cs$year))
treatment=rep("control",length(Cs$year))
tank=rep(1:24)
Cs=cbind(Cs,healthy,ploidy,treatment,tank)

# exposed sexuals
sET=tapply(sE$mic,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","total")

sEM=tapply(sE$mic,list(sE$year,sE$replicate),sum)
LsEM=melt(sEM)
names(LsEM)=c("year","replicate","mic")

Es<-cbind(LsET,LsEM$mic)
names(Es)=c("year","replicate","total","mic")
healthy=Es$total-Es$mic
ploidy=rep(2,length(Es$year))
treatment=rep("exposed",length(Es$year))
tank=rep(25:48)
Es=cbind(Es,healthy,ploidy,treatment,tank)

inf<-rbind(Ca,Cs,Ea,Es)

geefit1<-geeglm( cbind(mic,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy)+
                   treatment*factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
summary(geefit1)
wald.test(b = coef(geefit1), Sigma = vcov(geefit1), Terms = 14:16)

geefit2<-geeglm( cbind(mic,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
anova(geefit2,test="LRT")
summary(geefit2)
wald.test(b = coef(geefit2), Sigma = vcov(geefit2), Terms = 10)

geefit3<-geeglm( cbind(mic,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
anova(geefit3,test="LRT")

detach(df)

######
#Repeat: males and females; all parasites
# Experiment, males+females, all castrating parasites
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

# control asexuals
aCT=tapply(aC$castrated,list(aC$year,aC$replicate),length)
LaCT=melt(aCT)
names(LaCT)=c("year","replicate","total")

aCM=tapply(aC$castrated,list(aC$year,aC$replicate),sum)
LaCM=melt(aCM)
names(LaCM)=c("year","replicate","cas")

Ca<-cbind(LaCT,LaCM$cas)
names(Ca)=c("year","replicate","total","cas")
healthy=Ca$total-Ca$cas
ploidy=rep(3,length(Ca$year))
treatment=rep("control",length(Ca$year))
tank=rep(1:24)
Ca=cbind(Ca,healthy,ploidy,treatment,tank)

# exposed asexuals
aET=tapply(aE$castrated,list(aE$year,aE$replicate),length)
LaET=melt(aET)
names(LaET)=c("year","replicate","total")

aEM=tapply(aE$castrated,list(aE$year,aE$replicate),sum)
LaEM=melt(aEM)
names(LaEM)=c("year","replicate","cas")

Ea<-cbind(LaET,LaEM$cas)
names(Ea)=c("year","replicate","total","cas")
healthy=Ea$total-Ea$cas
ploidy=rep(3,length(Ea$year))
treatment=rep("exposed",length(Ea$year))
tank=rep(25:48)
Ea=cbind(Ea,healthy,ploidy,treatment,tank)

# control sexuals
sCT=tapply(sC$castrated,list(sC$year,sC$replicate),length)
LsCT=melt(sCT)
names(LsCT)=c("year","replicate","total")

sCM=tapply(sC$castrated,list(sC$year,sC$replicate),sum)
LsCM=melt(sCM)
names(LsCM)=c("year","replicate","cas")

Cs<-cbind(LsCT,LsCM$cas)
names(Cs)=c("year","replicate","total","cas")
healthy=Cs$total-Cs$cas
ploidy=rep(2,length(Cs$year))
treatment=rep("control",length(Cs$year))
tank=rep(1:24)
Cs=cbind(Cs,healthy,ploidy,treatment,tank)

# exposed sexuals
sET=tapply(sE$castrated,list(sE$year,sE$replicate),length)
LsET=melt(sET)
names(LsET)=c("year","replicate","total")

sEM=tapply(sE$castrated,list(sE$year,sE$replicate),sum)
LsEM=melt(sEM)
names(LsEM)=c("year","replicate","cas")

Es<-cbind(LsET,LsEM$cas)
names(Es)=c("year","replicate","total","cas")
healthy=Es$total-Es$cas
ploidy=rep(2,length(Es$year))
treatment=rep("exposed",length(Es$year))
tank=rep(25:48)
Es=cbind(Es,healthy,ploidy,treatment,tank)

inf<-rbind(Ca,Cs,Ea,Es)

geefit1<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy)+
                   treatment*factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
summary(geefit1)
wald.test(b = coef(geefit1), Sigma = vcov(geefit1), Terms = 14:16)

geefit2<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+treatment*factor(ploidy)+factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
anova(geefit2,test="LRT")
summary(geefit2)
wald.test(b = coef(geefit2), Sigma = vcov(geefit2), Terms = 10)

geefit3<-geeglm( cbind(cas,healthy) ~ treatment+factor(year)+factor(ploidy)+
                   treatment*factor(year)+factor(year)*factor(ploidy), 
                 family = binomial,data=inf,
                 id=tank,corstr="exchangeable")
anova(geefit3,test="LRT")

detach(df)

##################################
# Calculations
rm(list=ls())
df<-read.csv("Bins_adults.csv")
df<-subset(df,df$sex==0 & df$ploidy>0)
summary(df)
attach(df)

###################################################
# increase in infection from control to exposed tanks#
# control tanks
# all females, yes/no mic by tank, year
C<-subset(df,df$treatment==0)
tC=tapply(C$mic,list(C$replicate,C$year),length)
mC=tapply(C$mic,list(C$replicate,C$year),sum)
fC=mC/tC

# exposed tanks
# all females, yes/no mic by tank, year
E<-subset(df,df$treatment==1)
tE=tapply(E$mic,list(E$replicate,E$year),length)
mE=tapply(E$mic,list(E$replicate,E$year),sum)
fE=mE/tE

mean=c(mean(fC),mean(fE))
mean(fE)/mean(fC) # fold increase
se=c(sd(fC)/sqrt(24),sd(fE)/sqrt(24))

# sample sizes
mean=c(mean(tC),mean(tE))
se=c(sd(tC)/sqrt(24),sd(tE)/sqrt(24))

###################################################
# differential infection in control tanks#
Ca<-subset(df,df$treatment==0 &df$ploidy==3)
tCa=tapply(Ca$mic,list(Ca$replicate,Ca$year),length)
mCa=tapply(Ca$mic,list(Ca$replicate,Ca$year),sum)
fCa=mCa/tCa

Cs<-subset(df,df$treatment==0 & df$ploidy==2)
tCs=tapply(Cs$mic,list(Cs$replicate,Cs$year),length)
mCs=tapply(Cs$mic,list(Cs$replicate,Cs$year),sum)
fCs=mCs/tCs

m1=mean(fCa)
m2=mean(fCs)
se1=sd(fCa)/sqrt(24)
se2=sd(fCs)/sqrt(24)

fold<-mean(fCa)/mean(fCs) # fold increase
se_fold=fold*sqrt(((se1/m1)^2)+((se2/m2)^2))

# sample sizes
mean=c(mean(tCa),mean(tCs))
se=c(sd(tCa)/sqrt(24),sd(tCs)/sqrt(24))

###################################################
# differential infection in exposed tanks#
Ea<-subset(df,df$treatment==1 &df$ploidy==3)
tEa=tapply(Ea$mic,list(Ea$replicate,Ea$year),length)
mEa=tapply(Ea$mic,list(Ea$replicate,Ea$year),sum)
fEa=mEa/tEa

Es<-subset(df,df$treatment==1 & df$ploidy==2)
tEs=tapply(Es$mic,list(Es$replicate,Es$year),length)
mEs=tapply(Es$mic,list(Es$replicate,Es$year),sum)
fEs=mEs/tEs

m1=mean(fEa)
m2=mean(fEs)
se1=sd(fEa)/sqrt(24)
se2=sd(fEs)/sqrt(24)

fold<-mean(fEa)/mean(fEs) # fold increase
se_fold=fold*sqrt(((se1/m1)^2)+((se2/m2)^2))

# sample sizes
mean=c(mean(tEa),mean(tEs))
se=c(sd(tEa)/sqrt(24),sd(tEs)/sqrt(24))

# overall sample sizes
ta<-cbind(tCa,tEa)
mean(ta)
sd(ta)/sqrt(48)

ts<-cbind(tCs,tEs)
mean(ts)
sd(ts)/sqrt(48)

######## repeat for all parasites##########################

###################################################
# increase in infection from control to exposed tanks#
# control tanks
# all females, yes/no mic by tank, year
C<-subset(df,df$treatment==0)
tC=tapply(C$castrated,list(C$replicate,C$year),length)
mC=tapply(C$castrated,list(C$replicate,C$year),sum)
fC=mC/tC

# exposed tanks
# all females, yes/no mic by tank, year
E<-subset(df,df$treatment==1)
tE=tapply(E$castrated,list(E$replicate,E$year),length)
mE=tapply(E$castrated,list(E$replicate,E$year),sum)
fE=mE/tE

mean=c(mean(fC),mean(fE))
mean(fE)/mean(fC) # fold increase
se=c(sd(fC)/sqrt(24),sd(fE)/sqrt(24))

# sample sizes
mean=c(mean(tC),mean(tE))
se=c(sd(tC)/sqrt(24),sd(tE)/sqrt(24))

###################################################
# differential infection in control tanks#
Ca<-subset(df,df$treatment==0 &df$ploidy==3)
tCa=tapply(Ca$castrated,list(Ca$replicate,Ca$year),length)
mCa=tapply(Ca$castrated,list(Ca$replicate,Ca$year),sum)
fCa=mCa/tCa

Cs<-subset(df,df$treatment==0 & df$ploidy==2)
tCs=tapply(Cs$castrated,list(Cs$replicate,Cs$year),length)
mCs=tapply(Cs$castrated,list(Cs$replicate,Cs$year),sum)
fCs=mCs/tCs

mean=c(mean(fCa),mean(fCs))
mean(fCs)/mean(fCa) # fold increase
se=c(sd(fCa)/sqrt(24),sd(fCs)/sqrt(24))

m1=mean(fCa)
m2=mean(fCs)
se1=sd(fCa)/sqrt(24)
se2=sd(fCs)/sqrt(24)

fold<-mean(fCa)/mean(fCs) # fold increase
se_fold=fold*sqrt(((se1/m1)^2)+((se2/m2)^2))

# sample sizes
mean=c(mean(tCa),mean(tCs))
se=c(sd(tCa)/sqrt(24),sd(tCs)/sqrt(24))

###################################################
# differential infection in exposed tanks#
Ea<-subset(df,df$treatment==1 &df$ploidy==3)
tEa=tapply(Ea$castrated,list(Ea$replicate,Ea$year),length)
mEa=tapply(Ea$castrated,list(Ea$replicate,Ea$year),sum)
fEa=mEa/tEa

Es<-subset(df,df$treatment==1 & df$ploidy==2)
tEs=tapply(Es$castrated,list(Es$replicate,Es$year),length)
mEs=tapply(Es$castrated,list(Es$replicate,Es$year),sum)
fEs=mEs/tEs

mean=c(mean(fEa),mean(fEs))
mean(fEs)/mean(fEa) # fold increase
se=c(sd(fEa)/sqrt(24),sd(fEs)/sqrt(24))

m1=mean(fEa)
m2=mean(fEs)
se1=sd(fEa)/sqrt(24)
se2=sd(fEs)/sqrt(24)

fold<-mean(fEa)/mean(fEs) # fold increase
se_fold=fold*sqrt(((se1/m1)^2)+((se2/m2)^2))

# sample sizes
mean=c(mean(tEa),mean(tEs))
se=c(sd(tEa)/sqrt(24),sd(tEs)/sqrt(24))

detach(df)
