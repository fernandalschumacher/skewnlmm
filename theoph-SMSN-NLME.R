####################################################################################################
## SM(S)N Nonlinear Mixed Effects Model for Theoph data ##
####################################################################################################
## last updated in 2020-07-26
## downloaded from https://github.com/fernandalschumacher/skewnlmm

#loading packages
library(tidyverse)
library(reshape2)
library(ggplot2)
library(skewlmm) #optional. This package fits the linear version of the model, its documentation could be helpful to understand this code
library(mvtnorm)
library(nlme) #used for initial values

#loading SMSN-NLME functions
source("smsn-nlmm.R")
source("smn-nlmm.R")
source("loglik-NL.R")

## creating logistic function and derivative
derivlog1 <- function(x1,x2,beta1,beta2,beta3) x1*exp(-beta1+beta2+beta3)*
  (exp(-exp(beta3)*x2)-exp(-exp(beta2)*x2))/((exp(beta2)-exp(beta3))) #function without random effects, used for initial values

# nonlinear function to be fitted, 
# must have arguments x1,...,xp,beta1,...,betar,b1,...,bq
logist1<- function(x1,x2,beta1,beta2,beta3,b1,b2) { 
  x1*exp(-(beta1+b1)+(beta2+b2)+beta3)*
    (exp(-exp(beta3)*x2)-exp(-exp(beta2+b2)*x2))/((exp(beta2+b2)-exp(beta3)))
}
# derivative of the nonlinear function to be used in the linearization, 
# must have the same arguments as the previous function
der1logist <- deriv( ~x1*exp(-(beta1+b1)+(beta2+b2)+beta3)*
                       (exp(-exp(beta3)*x2)-exp(-exp(beta2+b2)*x2))/((exp(beta2+b2)-exp(beta3))),
                     c("beta1","beta2","beta3","b1","b2"), function(x1,x2,beta1,beta2,beta3,b1,b2){})

###################################################################################################
## loading the data set
data("Theoph")
help(Theoph)
theophdf <- as.data.frame(Theoph)
#
ggplot(theophdf,aes(x=Time,y=conc,group=Subject,alpha=I(.4)))+geom_line()+geom_point()+
  ylab("Theophylline concentration (mg/L)")+xlab("Time since drug administration (hr)")

#####################################################################################
# fitting
#####################################################################################
# fitting nlme for getting initial values
fit_nlme=nlme(conc~derivlog1(Dose,Time,beta1,beta2,beta3),data=Theoph,
              fixed=beta1+beta2+beta3~1,random=beta1+beta2~1|Subject,
              start=c(beta1=-3,beta2=.8,beta3=-2.45))
summary(fit_nlme)
# initial values
beta1=as.numeric(fit_nlme$coefficients$fixed)
sigmae=fit_nlme$sigma^2
D1=var(ranef(fit_nlme))
lambda=rep(3,2)*as.numeric(sign(moments::skewness(ranef(fit_nlme))))

##############################################################################
# fitting the model
##############################################################################
xmat<-cbind(theophdf$Dose,theophdf$Time)

############################# normal fit
set.seed(112)
fitN<- EM.sim.NL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,
                 nlfder=der1logist,nlf=logist1,beta1=beta1+runif(3,-.5,.5),
                 sigmae=sigmae,
                 D1=D1,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                 precisao=1e-4,max.iter=500,showiter=T,showerroriter=T)
rbind(fitN$theta,
      fitN$std.error)
#
colnames(fitN$random.effects) <- c("b0","b1")
ggplot(as.data.frame(fitN$random.effects),aes(sample = b0)) + stat_qq() + stat_qq_line()+
  ylab(expression(b[1]))
ggplot(as.data.frame(fitN$random.effects),aes(sample = b1)) + stat_qq() + stat_qq_line()+
  ylab(expression(b[2]))
#

fitlogN<- data.frame(value=theophdf$conc,time=theophdf$Time,
                     ind=theophdf$Subject,
                     fitted=fitN$fitted)
ggplot(fitlogN,aes(x=time,y=value))+geom_line()+geom_point() +
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~ind,scales = "free_y")

############################# skew-normal fit
fitSN<- EM.Skew.NL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,
                   nlfder=der1logist,nlf=logist1,beta1=fitN$estimates$beta,
                   sigmae=sigmae,lambda=lambda,
                   D1=D1,distr="sn",nu=NULL,lb=NULL,lu=NULL,
                   precisao=1e-4,max.iter=500,showiter=T,showerroriter=T)
rbind(fitSN$theta,
      fitSN$std.error)

############################# skew-t fit
fitST<- EM.Skew.NL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,
                   nlfder=der1logist,nlf=logist1,beta1=fitSN$estimates$beta,
                   sigmae=fitSN$estimates$sigma2,lambda=fitSN$estimates$lambda,
                   D1=fitSN$estimates$D,distr="st",nu=10,lb=1.1,lu=100,
                   #bi=fitSN$random.effects,
                   precisao=1e-2,max.iter=500,showiter=T,showerroriter=T)
rbind(fitST$theta,
      fitST$std.error)

fitlogST<- data.frame(value=theophdf$conc,time=theophdf$Time,
                      ind=theophdf$Subject,
                      fitted=fitST$fitted)
ggplot(fitlogST,aes(x=time,y=value))+geom_line()+geom_point() +
  geom_line(aes(y=fitted,color=I("blue"))) + facet_wrap(~ind)

# creating Healy-type plots
healy.plot(fitN,ind=theophdf$Subject,y=theophdf$conc,x=xmat,distr = "sn",nlfder = der1logist) + ggtitle("N-NLME model")
healy.plot(fitST,ind=theophdf$Subject,y=theophdf$conc,x=xmat,distr = "st",nlfder = der1logist)+ ggtitle("ST-NLME model")

############################# ssl fit
fitSSL<- EM.Skew.NL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,
                    nlfder=der1logist,nlf=logist1,beta1=fitSN$estimates$beta,
                    sigmae=fitSN$estimates$sigma2,D1=fitSN$estimates$D,lambda=fitSN$estimates$lambda,
                    distr="ss",nu=5,lb=.6,lu=50,
                    precisao=1e-2,max.iter=500,showiter=T,showerroriter=T)
rbind(fitSSL$theta,
      fitSSL$std.error)

############################# scn fit
fitSCN<- EM.Skew.NL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,
                    nlfder=der1logist,nlf=logist1,beta1=fitSN$estimates$beta,
                    sigmae=fitSN$estimates$sigma2,D1=fitSN$estimates$D,lambda=fitSN$estimates$lambda,
                    distr="scn",nu=c(.05,.8),lb=rep(.01,2),lu=rep(.99,2),
                    precisao=1e-2,max.iter=500,showiter=T,showerroriter=T)
rbind(fitSCN$theta,
      fitSCN$std.error)


##############################################
#selection criteria
(l1 <- logveroNL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,nlfder=der1logist,
                 fitObj = fitN,distr="sn")) #normal
(l2 <- logveroNL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,nlfder=der1logist,
                 fitObj = fitSN,distr="sn")) #sn
(l3 <- logveroNL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,nlfder=der1logist,
                 fitObj = fitST,distr="st")) #st
(l4 <- logveroNL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,nlfder=der1logist,
                 fitObj = fitSSL,distr="ss")) #ssl
(l5 <- logveroNL(y=theophdf$conc,x=xmat,ind=theophdf$Subject,nlfder=der1logist,
                 fitObj = fitSCN,distr="scn")) #scn

N=nrow(theophdf)
data.frame(distr=c("N","SN","ST","SSL","SCN"),
                   loglik=c(l1,l2,l3,l4,l5) #loglik
                   ,AIC=c( #AIC
                     -2*l1+2*length(fitN$theta),
                     -2*l2+2*length(fitSN$theta),
                     -2*l3+2*length(fitST$theta),
                     -2*l4+2*length(fitSSL$theta),
                     -2*l5+2*length(fitSCN$theta)),
                   # BIC
                   BIC=c(-2*l1+log(N)*length(fitN$theta),
                         -2*l2+log(N)*length(fitSN$theta),
                         -2*l3+log(N)*length(fitST$theta),
                         -2*l4+log(N)*length(fitSSL$theta),
                         -2*l5+log(N)*length(fitSCN$theta)))
