## JOINT MODELS AS LGMS
#Toy example
N=100 #number of patients
N_i=round(runif(N,min=5,max=10)) #Number of longitudinal observations per person
#Generate time-dependent exponential event times with iid shared random effects term
#Survival component
age=round(runif(N,min=18,max=75))
u_i=1+rnorm(N,mean=0,sd=1) #iid random effect
coeff=0.5 #share parameter
lambda=exp(coeff*u_i+0.01*age)
s_i=rexp(N,rate=lambda) #exponential survival times
#library(purrr)
#c_i=rbernoulli(N,p=0.9)#censoring
c_i=rep(1,N) #no censoring


#Longitudinal component
t=rep(NA,sum(N_i))
ID=rep(NA,sum(N_i))
a=rep(NA,sum(N_i))
t[1:N_i[1]]=runif(N_i[1],min=0,max=s_i[1]) # observation times 
ID[1:N_i[1]]=1
a[1:N_i[1]]=u_i[1]

for (i in 2:N){
t[(sum(N_i[1:(i-1)])+1):(sum(N_i[1:(i-1)])+N_i[i])]=runif(N_i[i],min=0,max=s_i[i]) #observation times
ID[(sum(N_i[1:(i-1)])+1):(sum(N_i[1:(i-1)])+N_i[i])]=i #person ID
a[(sum(N_i[1:(i-1)])+1):(sum(N_i[1:(i-1)])+N_i[i])]=u_i[i] #iid random effect
}
y=rnorm(sum(N_i),mean=(t^(-0.25)+a),sd=0.1) #longitudinal response with N_i[i] observations per person and a t^2 trend
plot(t,y) #view longitudional data

#Set-up for joint model
ng=length(y)
ns=length(s_i)
linear.covariate <- data.frame(muJ = as.factor(c(rep(1, ng), rep(2, ns))),
                               ltime = c(t,rep(0,ns)),
                               stime=c(rep(0,ng),s_i),
                               ages=c(rep(0,ng),age))

random.covariate<-list(IDs=c(rep(NA,ng),1:N),
                       IDl=c(ID,rep(NA,ns)),
                       V1 = c(t,rep(NA,ns)))


#Pad the covariates with NA for random effects
joint.data <- c(linear.covariate,random.covariate)
Yjoint=list(c(y,rep(NA,ns)),inla.surv(time=c(rep(NA,ng),s_i),event=c(rep(NA,ng),c_i)))
joint.data$Y <- Yjoint

#model with INLA
formulaJ = Y ~ muJ+ages+f(inla.group(V1,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01)))) +
  f(IDl, model = "iid",hyper=list(prec=list(prior="pc.prec",param=c(1,0.01))))+f(IDs,copy="IDl",fixed=FALSE)
inlares = inla(formulaJ, family = c("gaussian","exponentialsurv"),
               data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE,config=TRUE))
#non-linear trend only
inlares1 = inla(Y ~ muJ+ages+f(inla.group(V1,n=50),model="rw2",scale.model = TRUE,hyper = list(prec = list(prior="pc.prec", param=c(1, 0.01)))), family = c("gaussian","exponentialsurv"),
               data = joint.data, verbose=TRUE, control.compute=list(dic=TRUE,config=TRUE))
#jplm
library(JointModel)
fit0 <- jplm(y ~  (1 + t|ID), 
             nlm.par=t, data.y=data.frame(y=y,t=t,ID=ID), 
             Surv(s_i, c_i) ~ age, 
             formula.frailty= ~ 1 + s_i, 
             id.vec=(1:N), transf.par=0, data.surv=data.frame(s_i=s_i,c_i=c_i,age=age), n.knots=10)
summary(fit0)
pts <- seq(0.05,0.95,0.02)*max(t)  
out <- pred.jplm.nonlinear(fit0, at=pts, CI=TRUE, nlm.par=t)


#Check linear predictors
plot(a,inlares$summary.linear.predictor$mean[1:ng]) #longitudinal predictor
plot((coeff*u_i),inlares$summary.linear.predictor$mean[(ng+1):(ng+ns)]) #survival predictor

summary(inlares) #final estimates of joint model

#Data illustration
plot(t,y,col="lightgrey")
points(t,inlares1$summary.fitted.values[1:ng,1],col="red",cex=0.5,pch=2)
points(t,t^(-0.25)+0.01*mean(age),col="purple",cex=0.5,pch=3)
points(pts,out[,2],col="blue",cex=0.5,pch=5)
legend(x=1.5,y=11,legend=c("Data","Spline - INLA","True trajectory","Spline - JPLM"),col=c("lightgrey","red","purple","blue"),
pch=c(1,2,3,5))
