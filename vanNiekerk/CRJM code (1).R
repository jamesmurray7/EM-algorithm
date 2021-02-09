
#Competing risks simulated example
N=1000 #number of patients
N_i=round(runif(N,min=10,max=15)) #Number of longitudinal observations per person
#Generate time-dependent exponential event times with iid shared random effects term
#Survival component
u_i=1+rnorm(N,mean=0,sd=1) #iid random effect
coeff=0.5 #share parameter
lambda=exp(coeff*u_i)
s_i=rexp(N,rate=lambda) #exponential survival times
#library(purrr)
#c_i=rbernoulli(N,p=0.9)#censoring
c_i=rep(1,N) #no censoring
age=round(runif(N,15,75))
eta1=0.3*u_i+0.01*age
eta2=-0.1*u_i+0.02*age
eta3=0.2*u_i+0.0003*age
lambda1=exp(eta1)
lambda2=exp(eta2)
lambda3=exp(eta3)

IDs=1:N
time1=rexp(n=N,rate=lambda1)
Cause=rbinom(n=N,3,0.6)
time=rep(NA,N)
for (i in 1:N)
{if (Cause[i]==1) {time[i]=rexp(n=1,rate=lambda1[i])}
  if (Cause[i]==2) {time[i]=rexp(n=1,rate=lambda2[i])}
  if (Cause[i]==3) {time[i]=rexp(n=1,rate=lambda3[i])}
  if (Cause[i]==0) {time[i]=rexp(n=1,rate=1)}
}

dataS<-data.frame(ID=IDs,C=Cause,Time=time,Age=age)
data1<-dataS
dataE1<-data1
dataE1$event<-dataE1$C
dataE1$event[dataE1$event!=1]<-0
dataE2<-data1
dataE2$event<-dataE2$C
dataE2$event[dataE2$event!=2]<-0
dataE2$event<-dataE2$event/2
dataE3<-data1
dataE3$event<-dataE3$C
dataE3$event[dataE3$event!=3]<-0
dataE3$event<-dataE3$event/3

#Longitudinal component
t=rep(NA,sum(N_i))
ID=rep(NA,sum(N_i))
a=rep(NA,sum(N_i))
age_i=rep(NA,sum(N_i))
t[1:N_i[1]]=runif(N_i[1],min=0,max=s_i[1]) # observation times 
ID[1:N_i[1]]=1
a[1:N_i[1]]=u_i[1]
age_i[1:N_i[1]]=age[1]

for (i in 2:N){
  t[(sum(N_i[1:(i-1)])+1):(sum(N_i[1:(i-1)])+N_i[i])]=runif(N_i[i],min=0,max=s_i[i]) #observation times
  ID[(sum(N_i[1:(i-1)])+1):(sum(N_i[1:(i-1)])+N_i[i])]=i #person ID
  a[(sum(N_i[1:(i-1)])+1):(sum(N_i[1:(i-1)])+N_i[i])]=u_i[i] #iid random effect
  age_i[(sum(N_i[1:(i-1)])+1):(sum(N_i[1:(i-1)])+N_i[i])]=age[i] #iid random effect
}
y=rpois(n=sum(N_i),lambda=(exp(t^(1.2)+a))) #longitudinal response with N_i[i] observations per person and a t^2 trend
plot(t,y,pch=1,cex=0.3) #view longitudional data

#INLA on sim data
dataL<-data.frame(Time=t,Age=age_i,Y=y,ID=ID)
nL<-nrow(dataL)
nS<-nrow(dataS)

fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(2,nS),rep(3,nS),rep(4,nS))),
                      ageL=c(dataL$Age,rep(0,nS),rep(0,nS),rep(0,nS)),
                      age1=c(rep(0,nL),dataS$Age,rep(0,nS),rep(0,nS)),
                      age2=c(rep(0,nL),rep(0,nS),dataS$Age,rep(0,nS)),
                      age3=c(rep(0,nL),rep(0,nS),rep(0,nS),dataS$Age),
                      Ltime=c(dataL$Time,rep(0,nS),rep(0,nS),rep(0,nS)),
                      C1time=c(rep(0,nL),dataS$Time,rep(0,nS),rep(0,nS)),
                      C2time=c(rep(0,nL),rep(0,nS),dataS$Time,rep(0,nS)),
                      C3time=c(rep(0,nL),rep(0,nS),rep(0,nS),dataS$Time))
random.eff<-list(timeL=c(dataL$Time,rep(NA,nS),rep(NA,nS),rep(NA,nS)),
                 Lr1=c(dataL$ID,rep(NA,nS),rep(NA,nS),rep(NA,nS)),
                 Lr2=c(dataL$ID,rep(NA,nS),rep(NA,nS),rep(NA,nS)),
                 C1r1=c(rep(NA,nL),dataS$ID,rep(NA,nS),rep(NA,nS)),
                 C1r2=c(rep(NA,nL),dataS$ID,rep(NA,nS),rep(NA,nS)),
                 C2r1=c(rep(NA,nL),rep(NA,nS),dataS$ID,rep(NA,nS)),
                 C2r2=c(rep(NA,nL),rep(NA,nS),dataS$ID,rep(NA,nS)),
                 C3r1=c(rep(NA,nL),rep(NA,nS),rep(NA,nS),dataS$ID),
                 C3r2=c(rep(NA,nL),rep(NA,nS),rep(NA,nS),dataS$ID))


jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$Y,rep(NA, nS),rep(NA,nS),rep(NA,nS))
y.survC1 <- inla.surv(time = c(rep(NA, nL),dataS$Time,rep(NA,nS),rep(NA,nS)), event = c(rep(NA, nL),dataE1$event,rep(NA,nS),rep(NA,nS)))
y.survC2 <- inla.surv(time = c(rep(NA, nL), rep(NA,nS),dataS$Time, rep(NA,nS)), event = c(rep(NA, nL),rep(NA,nS),dataE2$event,rep(NA,nS)))
y.survC3 <- inla.surv(time = c(rep(NA, nL), rep(NA,nS), rep(NA,nS),dataS$Time), event = c(rep(NA, nL),rep(NA,nS),rep(NA,nS),dataE3$event))

y.joint<-list(y.long,y.survC1,y.survC2,y.survC3)

jointdata$Y=y.joint

#Model fit
formula.model=Y~-1+mu+age1+age2+age3+
  f(inla.group(timeL,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(Lr1, model="iid")+ 
  f(C1r1, copy="Lr1",fixed=FALSE)+
  f(C2r1, copy="Lr1",fixed=FALSE)+
  f(C3r1, copy="Lr1",fixed=FALSE)
 
Jointmodel= inla(formula.model, family = c("poisson","exponentialsurv","exponentialsurv","exponentialsurv"),
                 data = jointdata, verbose=TRUE)

summary(Jointmodel)

plot(t,y,pch=1,cex=0.3,ylab="Count",xlab="Time") #view longitudional data
lines(Jointmodel$summary.random$`inla.group(timeL, n = 50)`[,1],exp(Jointmodel$summary.random$`inla.group(timeL, n = 50)`[,2]),col="blue",lwd=2)

########################################################################
#SANAD - epileptic
library(joineR)
library(INLA)

mTime=max(max(epileptic$time),max(epileptic$with.time))
epileptic$time<-epileptic$time/mTime
epileptic$with.time <-epileptic$with.time/mTime

epileptic$interaction <- with(epileptic, time * (treat == "LTG"))
epileptic$interaction2 <- with(epileptic, time * (treat == "CBZ"))
epileptic$interaction[epileptic$interaction==0]<-NA
epileptic$interaction2[epileptic$interaction2==0]<-NA

data1<-epileptic
dataL<-data.frame(ID=data1$id,LDose=log(data1$dose+0.1),Time=data1$time,TimeSp=data1$time,Age=data1$age,Gender=data1$gender,LD=data1$learn.dis,Treatment=data1$treat,InteractionLTG=data1$interaction,InteractionCBZ=data1$interaction2,list(V=data1$id,W=data1$id),Dose=data1$dose)
dataS<-data.frame(ID=data1$id,Time=as.numeric(data1$with.time),Status=data1$with.status2,StatusISC=data1$with.status.isc,StatusUAE=data1$with.status.uae,Age=data1$age,Gender=data1$gender,LD=data1$learn.dis,Treatment=data1$treat)
dataS<-subset(dataS,!duplicated(dataS$ID))
summary(dataL)
summary(dataS)

#Visualize
plot(dataL$Time[dataL$ID==42],dataL$Dose[dataL$ID==42],type="n",xlim=c(0,1),ylim=c(0,10))
for (i in unique(dataL$ID)){
  lines(dataL$Time[dataL$ID==i],dataL$Dose[dataL$ID==i],col=i)
}

dataLCBZ=dataL[dataL$Treatment=="CBZ",]
dataLLTG=dataL[dataL$Treatment=="LTG",]

#Mean trajectories
modCBZ<-inla(formula=dataLCBZ$Dose~f(inla.group(dataLCBZ$Time,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01)))),data=dataLCBZ,family="gaussian",control.compute=list(dic=TRUE,config=TRUE))
modLTG<-inla(formula=dataLLTG$Dose~f(inla.group(dataLLTG$Time,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01)))),data=dataLLTG,family="gaussian",control.compute=list(dic=TRUE,config=TRUE))
#Plot means
plot(modCBZ$summary.random$`inla.group(dataLCBZ|S|Time, n = 50)`[,1]*2400,log(modCBZ$summary.fixed[1,1]+modCBZ$summary.random$`inla.group(dataLCBZ|S|Time, n = 50)`[,2]),type="l",lwd=2,ylim=c(0.5,1.2),xlab="Time",ylab="log(Dose)",xlim=c(0,2400))
lines(modLTG$summary.random$`inla.group(dataLLTG|S|Time, n = 50)`[,1]*2400,log(modLTG$summary.fixed[1,1]+modLTG$summary.random$`inla.group(dataLLTG|S|Time, n = 50)`[,2]),lty=2,lwd=2)
legend(x=2000,y=0.7,legend=c("CBZ","LTG"),lty=c(1,2),lwd=c(2,2))

#Joint model
#Pre-work for entire predictor
nL<-nrow(dataL)
nS<-nrow(dataS)
fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(1,nL),rep(2,nS),rep(3,nS))),
                      ageL=c(dataL$Age,dataL$Age,rep(0,nS),rep(0,nS)),
                      ageUAE=c(rep(0,nL),rep(0,nL),dataS$Age,rep(0,nS)),
                      ageISC=c(rep(0,nL),rep(0,nL),rep(0,nS),dataS$Age),
                      treatmentL=as.factor(c(dataL$Treatment,dataL$Treatment,rep(NA,nS),rep(NA,nS))),
                      treatmentUAE=as.factor(c(rep(NA,nL),rep(NA,nL),dataS$Treatment,rep(NA,nS))),
                      treatmentISC=as.factor(c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$Treatment)),
                      genderL=as.factor(c(dataL$Gender ,dataL$Gender ,rep(NA,nS),rep(NA,nS))),
                      genderUAE=as.factor(c(rep(NA,nL),rep(NA,nL),dataS$Gender ,rep(NA,nS))),
                      genderISC=as.factor(c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$Gender)))
random.eff<-list(timeL_LTG=c(dataL$InteractionLTG,rep(NA,nL),rep(NA,nS),rep(NA,nS)),
                 timeL_CBZ=c(dataL$InteractionCBZ,rep(NA,nL),rep(NA,nS),rep(NA,nS)),
                 linpredL=c(rep(NA,nL),dataL$ID,rep(NA,nS),rep(NA,nS)),
                 linpredL2=c(rep(NA,nL),rep(-1,nL),rep(NA,nS),rep(NA,nS)),
                 betaUAE=c(rep(NA,nL),rep(NA,nL),dataS$ID,rep(NA,nS)),
                 betaISC=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ID),
                 frailtyUAE=c(rep(NA,nL),rep(NA,nL),dataS$ID,rep(NA,nS)),
                 frailtyISC=c(rep(NA,nL),rep(NA,nL),rep(NA,nS),dataS$ID))
                      

jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$Dose,rep(NA,nL),rep(NA, nS),rep(NA,nS))
y.eta<-c(rep(NA,nL),rep(0,nL),rep(NA,nS),rep(NA,nS))
y.survUAE <- inla.surv(time = c(rep(NA, nL), rep(NA,nL),dataS$Time,rep(NA,nS)), event = c(rep(NA, nL),rep(NA,nL),dataS$StatusUAE,rep(NA,nS)))
y.survISC <- inla.surv(time = c(rep(NA, nL), rep(NA,nL),rep(NA,nS),dataS$Time), event = c(rep(NA, nL),rep(NA,nL),rep(NA,nS),dataS$StatusISC))
y.joint<-list(y.long,y.eta,y.survUAE,y.survISC)

jointdata$Y=y.joint

#Model fit
formula.model=Y~treatmentL+treatmentUAE+treatmentISC+f(inla.group(timeL_LTG,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
  f(inla.group(timeL_CBZ,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01))))+
f(linpredL, linpredL2, model="iid", hyper = list(prec = list(initial = -6, fixed=TRUE))) + 
  f(betaUAE, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))+
  f(betaISC, copy="linpredL", hyper = list(beta = list(fixed = FALSE)))
    

Jointmodel= inla(formula.model, family = c("gaussian","gaussian","weibullsurv","weibullsurv"),
                      data = jointdata, verbose=TRUE, control.compute=list(dic=TRUE),
                      control.family = list(
                        list(),
                        list(hyper = list(prec = list(initial = 10, fixed=TRUE))),
                        list(),
                        list()
                        )
                 )

summary(Jointmodel)

##Random effects only
fixed.eff<-data.frame(mu=as.factor(c(rep(1,nL),rep(2,nS),rep(3,nS))),
                      ageL=c(dataL$Age,rep(0,nS),rep(0,nS)),
                      ageUAE=c(rep(0,nL),dataS$Age,rep(0,nS)),
                      ageISC=c(rep(0,nL),rep(0,nS),dataS$Age),
                      treatmentL=as.factor(c(dataL$Treatment,rep(NA,nS),rep(NA,nS))),
                      treatmentUAE=as.factor(c(rep(NA,nL),dataS$Treatment,rep(NA,nS))),
                      treatmentISC=as.factor(c(rep(NA,nL),rep(NA,nS),dataS$Treatment)),
                      genderL=as.factor(c(dataL$Gender ,rep(NA,nS),rep(NA,nS))),
                      genderUAE=as.factor(c(rep(NA,nL),dataS$Gender ,rep(NA,nS))),
                      genderISC=as.factor(c(rep(NA,nL),rep(NA,nS),dataS$Gender)),
                      UAETime=c(rep(0,nL),dataS$Time,rep(0,nS)),
                      ISCTime=c(rep(0,nL),rep(0,nS),dataS$Time),
                      Ltime=c(dataL$Time,rep(0,nS),rep(0,nS)))
random.eff<-list(timeL=c(dataL$Time,rep(NA,nS),rep(NA,nS)),
                 Lr1=c(dataL$ID,rep(NA,nS),rep(NA,nS)),
                 Lr2=c(dataL$ID,rep(NA,nS),rep(NA,nS)),
                 UAEr1=c(rep(NA,nL),dataS$ID,rep(NA,nS)),
                 UAEr2=c(rep(NA,nL),dataS$ID,rep(NA,nS)),
                 ISCr1=c(rep(NA,nL),rep(NA,nS),dataS$ID),
                 ISCr2=c(rep(NA,nL),rep(NA,nS),dataS$ID),
                 frailtyUAE=c(rep(NA,nL),dataS$ID,rep(NA,nS)),
                 frailtyISC=c(rep(NA,nL),rep(NA,nS),dataS$ID))


jointdata<-c(fixed.eff,random.eff)
y.long <- c(dataL$Dose,rep(NA, nS),rep(NA,nS))
y.survUAE <- inla.surv(time = c(rep(NA, nL),dataS$Time,rep(NA,nS)), event = c(rep(NA, nL),dataS$StatusUAE,rep(NA,nS)))
y.survISC <- inla.surv(time = c(rep(NA, nL), rep(NA,nS),dataS$Time), event = c(rep(NA, nL),rep(NA,nS),dataS$StatusISC))
y.joint<-list(y.long,y.survUAE,y.survISC)

jointdata$Y=y.joint

#Model fit
formula.model=Y~treatmentL+treatmentUAE+treatmentISC+
  f(Lr1, model="iid2d",n=2*(nL+nS+nS)) + 
  f(Lr2, Ltime,copy="Lr1",fixed=TRUE)+
  f(UAEr1,copy="Lr1",fixed=FALSE)+f(UAEr2,UAETime,copy="UAEr1",fixed=FALSE)+
  f(ISCr1,copy="Lr1",fixed=FALSE)+f(ISCr2,ISCTime,copy="ISCr1",fixed=FALSE)

Jointmodel= inla(formula.model, family = c("gaussian","weibullsurv","weibullsurv"),
                 data = jointdata, verbose=TRUE, control.compute=list(dic=TRUE))

summary(Jointmodel)