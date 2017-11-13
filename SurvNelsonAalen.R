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



# diff.vec0<-c(NAAestim0[1],diff(NAAestim0))
# var.NAA0<-cumsum(diff.vec0*(diffvec0-1))
# diff.vec1<-c(NAAestim1[1],diff(NAAestim1))
# var.NAA1<-cumsum(diff.vec1^2)
# diff.vec.all<-c(NAAestimator[1],diff(NAAestimator))
# var.NAA.all<-cumsum(diff.vec.all^2)
# 
# matplot(c(NAAestim0,NAAestim0-2*sqrt(var.NAA),NAAestim0+2*sqrt(var.NAA)))
# 
# 
# plot(mess_trt0$intseiz1,NAAestim0,ylab='Cumulative hazard', ylim=c(0,1),xlab='Time (in days)',col=2)
# points(NAAestim0*exp(-2*sqrt(var.NAA)/NAAestim0))
# points(NAAestim0*exp(2*sqrt(var.NAA)/NAAestim0))


###scaling#############################################

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

par(mfrow=c(1,3))
surv.NAA.confint(with(drughiv.noties, Surv(time, delta)))
surv.NAA.confint(surv2)
surv.NAA.confint(surv1)



mfit2 <- survfit(survobj ~ 1, conf.type = "log-log")
my.fit2 <- summary(mfit2)
h.sort.of2 <- my.fit2$n.event / my.fit2$n.risk

scaledt<-my.fit2$time/max(my.fit2$time)
n<-sum(my.fit2$n.event)
n*scaledt
