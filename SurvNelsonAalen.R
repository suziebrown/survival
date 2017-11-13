library(asaur)
library(KMsurv)
library(survival)
library(mice)

mess <- read.csv("MESSdat.csv", sep="\t", header=T)
data(drughiv)   #34 obs
data(bcdeter)   #95 obs
data(alloauto)
data(bfeed)
data(prostateSurvival)  # 14294 (but we sould check the status variable)

###MESS DATA############################################################################

mess$trt<-as.factor(mess$trt)
mess$sex<-as.factor(mess$sex)

surv1 <- with(mess, Surv(intseiz1, indseiz1 == 1))
mfit <- survfit(surv1 ~ 1, conf.type = "log-log", data = mess)
#options(survfit.print.mean =T)
#plot(mfit, mark.time=T, xlab =factor(c(1,1)) "Days since randomisation",ylab = "Survival")
#summary(mfit,censored=T)
n.years <- seq(365.25, max(mess$intseiz1), by=365.25)
mfit.trt <- survfit(surv1 ~ trt,  data = mess)  #conf.type = "log-log",
plot(mfit.trt, conf.int = T, mark.time=F, col =  2:3, xlab = "Years since randomisation", ylab = "Survival probability" ,xscale = 365.25, cex=0.8)#,xaxt = 'n')
#axis(side=1, at=n.years,labels = 1:length(n.years))
legend("topright", legend = c("Immediate","Deferred"), col = 2:3, lty=1)

mess$intseiz1<-as.numeric(mess$intseiz1)
#mess <- mess[order(mess$intseiz1),] 
mess_trt0<-mess[mess$trt==0,]
mess_trt1<-mess[mess$trt==1,]
NAAestim0<-nelsonaalen(mess_trt0, intseiz1, indseiz1)
NAAestim1<-nelsonaalen(mess_trt1, intseiz1, indseiz1)
plot(mess_trt0$intseiz1,NAAestim0,ylab='Cumulative hazard', ylim=c(0,1),xlab='Time (in days)',col=2)
points(mess_trt1$intseiz1,NAAestim1,col=3)
legend("bottomright", legend = c("Immediate","Deferred"), col = 2:3, lty=1)




###NELSON-AALEN####################################################

surv2 <- with(alloauto, Surv(time, delta))
surv3 <- with(drughiv, Surv(time, delta))

surv.NAA<-function(survobj,...){
mfit <- survfit(survobj ~ 1, conf.type = "log-log")
fitsumm <- summary(mfit)
H.hat <- -log(fitsumm$surv)
H.hat <- c(H.hat, H.hat[length(H.hat)])
H.tilde<-cumsum(fitsumm$n.event/fitsumm$n.risk)
H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
nobs <- mfit$n
plot(c(fitsumm$time,tail(fitsumm$time,1)+1), H.hat, xlab='time', 
     ylab='cumulative hazard',
     main= paste("Dataset,", nobs ,"obs"), type='s',...)
points(c(fitsumm$time,tail(fitsumm$time,1)+1), H.tilde, lty=2, type='s')
legend("bottomright", legend=c("H.hat","Nelson-Aalen"), lty=1:2)
}

par(mfrow=c(1,3))
surv.NAA(surv3)
surv.NAA(surv2)
surv.NAA(surv1)


surv.NAA.confint<-function(survobj,...){
  mfit <- survfit(survobj ~ 1, conf.type = "log-log")
  fitsumm <- summary(mfit)
  H.hat <- -log(fitsumm$surv)
  H.hat <- c(H.hat, H.hat[length(H.hat)])
  H.tilde<-cumsum(fitsumm$n.event/fitsumm$n.risk)
  H.tilde <- c(H.tilde, H.tilde[length(H.tilde)])
  nobs <- mfit$n
  mfit.NAA <- survfit(survobj ~ 1, type='fleming')  #conf.type = "log-log",
  plot(mfit.NAA, conf.int = T, mark.time=F, col = 3, fun="cumhaz", xlab='time', 
       ylab='cumulative hazard', main= paste("Dataset,", nobs ,"obs"))
  points(c(fitsumm$time,tail(fitsumm$time,1)+1), H.tilde, col=3, type='s')
  points.default(c(fitsumm$time,tail(fitsumm$time,1)+1), H.hat, type='s',col=2,...)
  legend("bottomright", legend=c("Kaplan-Meier","Nelson-Aalen"), col=2:3,lty=1)
}


data(drughiv)  
survobj <- with(drughiv, Surv(time, delta))
surv.NAA(survobj)
mfit.drug <- survfit(Surv(time, delta) ~ 1,  data = drughiv, type='fleming')  #conf.type = "log-log",
time.ties<-summary(mfit.drug)$time[summary(mfit.drug)$n.event>1]
vec.timeties<-rep(NA,length(time.ties))
for(i in 1:length(time.ties)) vec.timeties[i]<-sample(which(drughiv$time == time.ties[i]),1)
drughiv.noties<-drughiv[-vec.timeties,]
mfit.drug.noties <- survfit(Surv(time, delta) ~ 1,  data = drughiv.noties, type='fleming')  #conf.type = "log-log",
plot(mfit.drug.noties, conf.int = T, mark.time=F, col =  1, fun="cumhaz",xlab = "Days since randomisation", ylab = "Cumulative Hazard",ylim=c(0,2.1))
surv.NAA(with(drughiv.noties, Surv(time, delta)))

par(mfrow=c(1,3))
surv.NAA.confint(with(drughiv.noties, Surv(time, delta)))
surv.NAA.confint(surv2)
surv.NAA.confint(surv1)

### Test ##################################################################

library(YPmodel)
data(gastric)
names(gastric)<-c("time","delta","trt")
survgastric <- with(gastric, Surv(time, delta))
surv.NAA(survgastric)
mfit.gastric.all <- survfit(survgastric ~ 1,  data = gastric, type='fleming')  #conf.type = "log-log",
mfit.gastric.trt <- survfit(survgastric ~ trt,  data = gastric, type='fleming')  #conf.type = "log-log",
par(mfrow=c(1,2))
plot(mfit.gastric, conf.int = F, mark.time=F, col = 1:2, xlab='time', 
     ylab='Survival')
legend("topright", legend=c("CT only","CT + RT"), col=1:2,lty=1)
plot(mfit.gastric, conf.int = F, mark.time=F, col = 1:2, fun="cumhaz", xlab='time', 
     ylab='cumulative hazard')
legend("bottomright", legend=c("CT only","CT + RT"), col=1:2,lty=1)


###renyi.test function#################################################################

renyi.test<-function(summary.pooled,summary.strata,type.test=c("two-sided","lower")){
library(plyr)
library(dplyr)
summ.all<-as.data.frame(with(summary(summary.pooled),cbind(time,n.risk,n.event)))
summ.trt1<-filter(as.data.frame(with(summary(summary.strata),cbind(time,n.risk,n.event,strata))),strata==1)[,1:3]
summ.trt2<-filter(as.data.frame(with(summary(summary.strata),cbind(time,n.risk,n.event,strata))),strata==2)[,1:3]
dftab<-join_all(list(summ.all, summ.trt1, summ.trt2), by='time', type='left')
names(dftab)[-1]<-c("n.riskALL","n.eventALL","n.risk1","n.event1","n.risk2","n.event2")
if(is.na(dftab[1,4])) dftab[1,4:5]<-c(summ.trt1$n.risk[1],0)
if(is.na(dftab[1,6])) dftab[1,6:7]<-c(summ.trt2$n.risk[1],0)
for(i in 2:dim(dftab)[1]){
  if(is.na(dftab[i,4])) dftab[i,4:5]<-c(dftab[i-1,4]-dftab[i-1,5],0)
  if(is.na(dftab[i,6])) dftab[i,6:7]<-c(dftab[i-1,6]-dftab[i-1,7],0)
}
  
z.ti<-cumsum(with(dftab,n.event1-n.risk1*(n.eventALL/n.riskALL)))
if(type.test=="two-sided") z.ti<-abs(z.ti)

plot(dftab$time,z.ti,type="l",xlab="Time",ylab="Z-statistic")
supZt<-max(z.ti)
supt<-dftab$time[which.max(z.ti)]
abline(v=supt,lty=3,col=2)

tau.endpoint<-tail(which(with(dftab,n.risk1*n.risk2>0)),1)
var.estim.endpoint<-cumsum(with(dftab,((n.risk1*n.risk2)/n.riskALL^2)*n.eventALL*(n.riskALL-n.eventALL)/(n.riskALL-1)))[tau.endpoint]
Q.stat<-supZt/sqrt(var.estim.endpoint)

prob.Qstat<-function(Q.stat){
  k<-0:100
  temp<-((-1)^k)/(2*k+1)
  temp<-temp*exp(-pi^2*(2*k+1)^2/(8*Q.stat^2))    ###dftab$time[tau.endpoint]
  1-4*sum(temp)/pi
}
p.val<-ifelse(type.test=="two-sided",prob.Qstat(Q.stat),2*pnorm(-Q.stat))
cat("\nRenyi test statistic Q =",Q.stat,"\nAlternative hypothesis:",type.test,"\np-value =",p.val)
}

renyi.test(summary.pooled=mfit.gastric.all,summary.strata=mfit.gastric.trt,type.test="two-sided")
renyi.test(summary.pooled=mfit.gastric.all,summary.strata=mfit.gastric.trt,type.test="lower")
