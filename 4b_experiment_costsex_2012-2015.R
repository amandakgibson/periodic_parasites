# Periodic, parasite-mediated selection for and against sex
# Gibson, Delph, Vergara, Lively

# Experimental test: parasites can increase the cost of sex; B
# Estimating the cost of sex in presence and absence of parasites, 2012-2015

######
# Experiment; 2012-2015; model fitting

rm(list=ls())
library(reshape)

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
treatment=rep(2,length(LaET$triploid))
ET=cbind(LsET,LaET$triploid,treatment)
names(ET)=c("year","replicate","diploid","triploid","treatment")

treatment=rep(1,length(LaCT$triploid))
CT=cbind(LsCT,LaCT$triploid,treatment)
names(CT)=c("year","replicate","diploid","triploid","treatment")

T=rbind(ET,CT)
sum=T$triploid+T$diploid
Tfbabies=data.frame(cbind(T$year,T$replicate,T$treatment,T$triploid,sum))
names(Tfbabies)=c("year","rep","treatment","y","z")

# adults
dfA<-subset(df,df$cohort=="adult")
# subset by ploidy and treatment
aC=subset(dfA,dfA$treatment==0 & dfA$ploidy==3)
aE=subset(dfA,dfA$treatment==1 & dfA$ploidy==3)
sC=subset(dfA,dfA$treatment==0 & dfA$ploidy==2)
sE=subset(dfA,dfA$treatment==1 & dfA$ploidy==2)

# control triploids
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
treatment=rep(1,length(LaET$triploid))
ET=cbind(LsET,LaET$triploid,treatment)
names(ET)=c("year","replicate","diploid","triploid","treatment")

treatment=rep(0,length(LaCT$triploid))
CT=cbind(LsCT,LaCT$triploid,treatment)
names(CT)=c("year","replicate","diploid","triploid","treatment")

T=rbind(ET,CT)
freq=T$triploid/(T$triploid+T$diploid)
Tfadults=cbind(T,freq)

asex=cbind(Tfbabies,Tfadults$freq)
names(asex)=c("year","rep","treatment","y","z","x")

detach(df)

asex$year<-as.factor(asex$year)
asex$treatment<-as.factor(asex$treatment)
asex$rep<-as.factor(asex$rep)

#sort by mpg (ascending) and cyl (descending)
asex <- asex[order(asex$treatment,asex$year),] 

attach(asex)

#####
# model selection
# open packages
library(bbmle)
library(emdbook)

# models
# model 1: No cost of sex
m1 = x

# model 2: 2-fold cost of sex
m2 = (x*2)/(1+x)

# model 3: c-fold cost of sex
m3 = function(F0, c){ # function of F0 (x in prior models) and c, the cost of sex
  (F0*c)/(1+F0*(c-1))}

# model 4: c-fold cost of sex with treatment variation
m4 = function(F0,c0,c.exp){ # c0 is our base cost in control, modifies for exposed
  cdiff=c(0,c.exp) # this value is added to our cost - nothing added for control
  cost=c0+cdiff[asex$treatment] # the cost that we'll use, with the addition of cdiff indexed by treatment
  m4.tr=(F0*cost)/(1+F0*(cost-1)) # the same function as above
}

# negative log-likelihood functions
# model 1
NLL.m1<--sum(dbinom(y,size=z,prob=m1,log=TRUE)) 
# model 1.bb
NLL.m1.bb<- function (theta){
  p.trip=m1
  -sum(dbetabinom(y,size=z,prob=p.trip,theta=theta,log=TRUE)) 
}

# model 2
NLL.m2<--sum(dbinom(y,size=z,prob=m2,log=TRUE))
# model 2.bb
NLL.m2.bb<- function (theta){
  p.trip=m2
  -sum(dbetabinom(y,size=z,prob=p.trip,theta=theta,log=TRUE)) 
}

# model 3
NLL.m3 = function (c){
  p.trip = m3(x,c)
  -sum(dbinom(y,size=z,prob=p.trip,log=TRUE))
}
# model 3.bb
NLL.m3.bb = function (c,theta){
  p.trip = m3(x,c)
  -sum(dbetabinom(y,size=z,prob=p.trip,theta=theta,log=TRUE))
}

# model 4
NLL.m4 = function (c0,c.exp){
  p.trip = m4(x,c0,c.exp)
  -sum(dbinom(y,size=z,prob=p.trip,log=TRUE))
}
# model 4.bb
NLL.m4.bb = function (c0,c.exp,theta){
  p.trip = m4(x,c0,c.exp)
  -sum(dbetabinom(y,size=z,prob=p.trip,theta=theta,log=TRUE))
}

# optimization
# model 1.bb 
FSP.m1.bb = mle2(NLL.m1.bb,start=list(theta=10)) 

# model 2.bb 
FSP.m2.bb = mle2(NLL.m2.bb,start=list(theta=10)) 

# model 3 
FSP.m3 = mle2(NLL.m3,start=list(c=2))
# model 3.bb 
FSP.m3.bb = mle2(NLL.m3.bb,start=list(c=2,theta=10))

# model 4 
FSP.m4 = mle2(NLL.m4,start=list(c0=2,c.exp=0))
# model 4.bb 
FSP.m4.bb = mle2(NLL.m4.bb,start=list(c0=2,c.exp=0,theta=10))

logL = c(NLL.m1,-logLik(FSP.m1.bb),NLL.m2,-logLik(FSP.m2.bb),-logLik(FSP.m3),-logLik(FSP.m3.bb),
         -logLik(FSP.m4),-logLik(FSP.m4.bb))
k = c(0,1,0,1,1,2,2,3)

# AICc values
# model 1
AICc1 = -2*(-NLL.m1)+2*0 + (2*0*(0+1))/(48-0-1)
# model 1.bb
AICc1.bb = -2*(logLik(FSP.m1.bb)[1])+2*1 + (2*1*(1+1))/(48-1-1)

# model 2
AICc2 = -2*(-NLL.m2)+2*0 + (2*0*(0+1))/(48-0-1)
# model 2.bb
AICc2.bb = -2*(logLik(FSP.m2.bb)[1])+2*1 + (2*1*(1+1))/(48-1-1)

# model 3
AICc3 = -2*(logLik(FSP.m3)[1])+2*1 + (2*1*(1+1))/(48-1-1)
# model 3.bb
AICc3.bb = -2*(logLik(FSP.m3.bb)[1])+2*2 + (2*2*(2+1))/(48-2-1)

# model 4
AICc4 = -2*(logLik(FSP.m4)[1])+2*2 + (2*2*(2+1))/(48-2-1)
# model 4.bb
AICc4.bb = -2*(logLik(FSP.m4.bb)[1])+2*3 + (2*3*(3+1))/(48-3-1)

AICc = c(AICc1,AICc1.bb,AICc2,AICc2.bb,AICc3,AICc3.bb,AICc4,AICc4.bb)

# delta AICs
#relative to model3
base = rep(AICc2.bb,8)

delta = AICc-base
weights = rep(0,8)
# Akaike weights - estimates of uncertainty
for (i in 1:8){
  weights[i]<-(exp(-0.5*delta[i])/(exp(-0.5*delta[1])+exp(-0.5*delta[2])+exp(-0.5*delta[3])+exp(-0.5*delta[4])+exp(-0.5*delta[5])+exp(-0.5*delta[6])+exp(-0.5*delta[7])+exp(-0.5*delta[8])))
}

# Summary table
summary = as.table(cbind(round(logL,2),k,round(AICc,2),round(delta,2),round(weights,2)))
colnames(summary) = c("logL", "K","AICc","delta AIC","Akaike weights")
rownames(summary)<-c("no cost","no cost, betabinom","2-fold cost","2-fold cost, betabinom","optimized cost",
                     "optimized cost, betabinom",
                     "cost by treatment","cost by treatment, betabinom")

# calculate AIC statistics for beta-binomial set only
# Beta-binomial values only!
logL = c(logLik(FSP.m1.bb),logLik(FSP.m2.bb),logLik(FSP.m3.bb),logLik(FSP.m4.bb))
k = c(1,1,2,3)
AICc = c(AICc1.bb,AICc2.bb,AICc3.bb,AICc4.bb)
base = rep(AICc4.bb,4)
delta = AICc-base
weights = rep(0,4)
for (i in 1:4){
  weights[i]<-(exp(-0.5*delta[i])/(exp(-0.5*delta[1])+exp(-0.5*delta[2])+exp(-0.5*delta[3])+exp(-0.5*delta[4])))
}
summary = as.table(cbind(round(logL,2),k,round(AICc,2),round(delta,2),round(weights,2)))
colnames(summary) = c("logL", "K","AICc","delta AIC","Akaike weights")
rownames(summary)<-c("no cost","2-fold cost","optimized cost",
                     "cost by treatment")

##############
# Model Selection Probabilities
m = 10000 
delta4 = rep(0,m)

delta = matrix(nrow=m,ncol=4)

for (i in 1:m) { 
  dboot<-asex[sample(1:48,48, replace=TRUE), ]
  NLL.m1.bb = function (theta){
    x=dboot$x
    p.trip=m1
    -sum(dbetabinom(dboot$y,size=dboot$z,prob=p.trip,theta=theta,log=TRUE)) #1 parameter to estimate
  }
  NLL.m2.bb = function (theta){
    x=dboot$x
    p.trip=m2
    -sum(dbetabinom(dboot$y,size=dboot$z,prob=p.trip,theta=theta,log=TRUE)) #1 parameter to estimate
  }
  NLL.m3.bb = function (c,theta){
    p.trip = m3(x,c)
    -sum(dbetabinom(dboot$y,size=dboot$z,prob=p.trip,theta=theta,log=TRUE))
  }
  NLL.m4.bb = function (c0,c.exp,theta){
    p.trip = m4(x,c0,c.exp)
    -sum(dbetabinom(dboot$y,size=dboot$z,prob=p.trip,theta=theta,log=TRUE))
  }
  
  FSP.m1.bb = mle2(NLL.m1.bb,start=list(theta=10)) 
  FSP.m2.bb = mle2(NLL.m2.bb,start=list(theta=10)) 
  FSP.m3.bb = mle2(NLL.m3.bb,start=list(c=2,theta=10)) # this result is sensitive to the starting value of theta
  FSP.m4.bb = mle2(NLL.m4.bb,start=list(c0=2,c.exp=0,theta=10))
  
  AICc1 = -2*(logLik(FSP.m1.bb)[1])+2*1 + (2*1*(1+1))/(48-1-1)
  AICc2 = -2*(logLik(FSP.m2.bb)[1])+2*1 + (2*1*(1+1))/(48-1-1)
  AICc3 = -2*(logLik(FSP.m3.bb)[1])+2*2 + (2*2*(2+1))/(48-2-1)
  AICc4 = -2*(logLik(FSP.m4.bb)[1])+2*3 + (2*3*(3+1))/(48-3-1)
  
  minimum = min(c(AICc1,AICc2,AICc3,AICc4))
  delta4[i] = AICc4 - minimum
  delta[i,] = c(AICc1-minimum,AICc2-minimum, AICc3-minimum,AICc4-minimum)
}

delta4 = sort(delta4,decreasing=FALSE)
delta4[m*0.95] 

###########3
# Parameter Estimations

# model3
m = 10000 
cboot=rep(0,m)

for (i in 1:m) { 
  dboot<-asex[sample(1:48,48, replace=TRUE), ]
  NLL.m3.bb = function (c,theta){
    p.trip = m3(x,c)
    -sum(dbetabinom(dboot$y,size=dboot$z,prob=p.trip,theta=theta,log=TRUE))
  }
  FSP.m3.bboot = mle2(NLL.m3.bb,start=list(c=2,theta=10)) # this result is sensitive to the starting value of theta
  cboot[i]=coef(FSP.m3.bboot)[1]
}

cboot = sort(cboot,decreasing=FALSE) 
low<-cboot[m*0.025]
high<-cboot[m*0.975] 

# model 4
m = 10000
cboot = rep(0,m)
c.expboot = rep(0,m) 

for (i in 1:m) { 
  control<-asex[sample(1:24,24, replace=TRUE), ]
  exposed<-asex[sample(25:48,24, replace=TRUE), ]
  
  dboot<-rbind(control, exposed)
  NLL.m4.bb = function (c0,c.exp,theta){
    p.trip = m4(x,c0,c.exp)
    -sum(dbetabinom(dboot$y,size=dboot$z,prob=p.trip,theta=theta,log=TRUE))
  }
  FSP.m4.bboot = mle2(NLL.m4.bb,start=list(c0=2,c.exp=0,theta=10))
  cboot[i]=coef(FSP.m4.bboot)[1]
  c.expboot[i]=coef(FSP.m4.bboot)[2]
}

cboot = sort(cboot,decreasing=FALSE) 
low<-cboot[m*0.025]
high<-cboot[m*0.975]  

c.expboot = sort(c.expboot,decreasing=FALSE)
c.explow<-c.expboot[m*0.025] 
c.exphigh<-c.expboot[m*0.975] 
