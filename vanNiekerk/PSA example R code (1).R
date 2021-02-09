## JOINT MODELS
library(JointModel)
data1 <- prostate
data2 <- dropout
ng<-nrow(data1)
ns<-nrow(data2)
data1=data1[order(data1$VisitTime),] #Just for "nice" graphs

## prepare the response variable
library(INLA)
y.long <- c(data1$logPSA.postRT , rep(NA, ns))
y.surv <- inla.surv(time = c(rep(NA, ng), data2$DropTime ), event = c(rep(NA, ng),data2$Status))
Yjoint <- list(y.long, y.surv)
N<-length(unique(data1$ID))


linear.covariate <- data.frame(mu = as.factor(c(rep(1, ng), rep(2, ns))),
                               b13.PSAbase = c(data1$logPSA.base, rep(0, ns)),
                               b22.PSABase2 = c(rep(0, ng), data2$logPSA.base2 ), 
                               b12.time = c(data1$VisitTime,rep(0,ns)),
                               b21.time=c(rep(0,ng),data2$DropTime))


random.covariate <- list(U11 = c(data1$ID,rep(NA, ns)),
                         V1 = c(data1$VisitTime,rep(NA,ns)),
                         U21 = c(data1$ID,rep(NA, ns)),
                         U12=c(rep(NA,ng),data2$ID2),
                         U22=c(rep(NA,ng),data2$ID2))


#Model 1 - intslope
subject=c(data1$ID,data2$ID2)
idx=1:(ng+ns)
z=c(data1$VisitTime,data2$DropTime)
strata=c(rep(1,ng),rep(2,ns))
joint.data <- c(linear.covariate,random.covariate,list(subject=subject,idx=idx,z=z,strata=strata))
joint.data$Y <- Yjoint
formulaM1 = Y ~ mu+f(inla.group(V1,n=N),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01)))) + b13.PSAbase + b22.PSABase2 +
  f(idx, model = "intslope",args.intslope = list(subject = subject,strata = strata,covariates = z))


Jointmodelres1 = inla(formulaM1, family = c("gaussian","weibullsurv"),
                 data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE,config=TRUE))

summary(Jointmodelres1)

#Model 2 - separate coefficients
joint.data <- c(linear.covariate,random.covariate)
joint.data$Y <- Yjoint

formulaM2 = Y ~ mu + f(inla.group(V1,n=N),model="rw2")+ b13.PSAbase + b22.PSABase2 +
  f(U11 , model="iid2d", param = c(4,1,1,-0.5), initial = c(-2.7,0.9,-0.22), n=2*(ng+ns)) +
  f(U21, b12.time, copy="U11",fixed=TRUE)+f(U12, copy="U11",fixed=FALSE)+f(U22,b21.time,copy="U12",fixed=FALSE)

Jointmodelres2 = inla(formulaM2, family = c("gaussian","weibullsurv"),
            data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))

summary(Jointmodelres2)

#Spline only
Jointmodelres3 = inla(Y ~ mu + f(inla.group(V1,n=N),model="rw2"), family = c("gaussian","weibullsurv"),
                      data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE))



#plot PSA
data1=data1[order(data1$VisitTime),]
plot(data1$VisitTime,data1$logPSA.postRT,xlab="VisitTime",ylab="log(PSA)") 
points(data1$VisitTime,Jointmodelres1$summary.fitted.values[1:ng,1],col='blue',lwd=1)
lines(data1$VisitTime,Jointmodelres3$summary.fitted.values[1:ng,1],col='red',lwd=5)

### Jointmodel with jplm
# joint fit of a partially linear model and a proportional hazards model 
# with a subject-specific random intercept and random slope
attach(prostate)
attach(dropout)
fit0 <- jplm(logPSA.postRT ~ logPSA.base + (1 + VisitTime|ID), 
             nlm.par=prostate$VisitTime, data.y=prostate, 
             Surv(DropTime, Status) ~ logPSA.base2, 
             formula.frailty= ~ 1 + DropTime, 
             id.vec=dropout$ID2, transf.par=0, data.surv=dropout, n.knots=5)
summary(fit0)
fit0$coef.lm.y
fit0$coef.nlm.y
fit0$coef.fixed.surv
fit0$coef.frailty.surv
fit0$loglik
fit0$var.resid
fit0$raneff.vcomp
pts <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)*max(prostate$VisitTime)  
out <- pred.jplm.nonlinear(fit0, at=pts, CI=TRUE, nlm.par=prostate$VisitTime)

#Comparison
plot(data1$VisitTime,data1$logPSA.postRT,xlab="VisitTime",ylab="log(PSA)",col="lightgrey") 
lines(data1$VisitTime,Jointmodelres3$summary.fitted.values[1:ng,1],col='red',lwd=5)
lines(pts,out$Value,col="purple",lwd=5)

