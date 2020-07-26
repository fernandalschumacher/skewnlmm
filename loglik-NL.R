####################################################################################################
## Log-likelihood function for SM(S)N Nonlinear Mixed Effects Model ##
####################################################################################################
## last updated in 2020-07-26
## downloaded from https://github.com/fernandalschumacher/skewnlmm

#function for computing the approximated log-likelihood function

logveroNL = function(y,#response vector (N x 1)
                     x,#covariates matrix (N x p1)
                     ind,#factor whose levels indicates the subjects or groups (N x 1) 
                     nlfder,#first derivative of the nonlinear function 
                     fitObj,#object returned from EM.Skew.NL or EM.sim.NL function
                     distr#distribution to be used: "sn" (normal), "st" (student's-t), "ss" (slash), or "scn" (contaminated normal)
                     ){
  fittedVal=fitObj$fitted
  biest=fitObj$random.effects
  beta1=fitObj$estimates$beta
  sigmae=fitObj$estimates$sigma2
  D1=fitObj$estimates$D
  lambda=fitObj$estimates$lambda
  nu=fitObj$estimates$nu
  #
  if (is.null(lambda)) lambda=rep(0,nrow(D1))
  
  m<-n_distinct(ind)
  N<-length(ind)
  p1<-dim(x)[2]
  p<-length(beta1)
  q1<- nrow(D1)
  
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda));
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #
  if (distr=="sn") {c. = -sqrt(2/pi)}
  if (distr=="st") {c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  
  ###
  Wtil<-matrix(nrow=N,ncol=p)
  Htil<-matrix(nrow=N,ncol=q1)
  ytil<-numeric(N)
  centerbi <- numeric(N)
  for (j in 1:nlevels(ind)) {
    jseq = ind==levels(ind)[j]
    y1=y[jseq]
    x1=matrix(x[jseq,  ],ncol=p1)
    ub1 <- biest[j,]#ub0[(((j-1)*q1)+1) : (j*q1), j]
    nj = length(y1)
    HWmat <- matrix(0,nrow=nj,ncol=(p+q1))
    fxtil <- matrix(0,nrow=nj,ncol=1)
    for(i in 1:nj)
    {
      farg <- as.list(c(x1[i,],beta1,ub1))
      names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep=""))
      formals(nlfder) <- farg
      fx <- nlfder()
      fxtil[i,] <- fx[1]
      HWmat[i,] <- attr(fx,"gradient")
    }
    
    Wtil[jseq,] <- matrix(HWmat[,1:p],nrow=nj,ncol=p)            # W til
    Htil[jseq,] <- matrix(HWmat[,(p+1):(p+q1)],nrow=nj,ncol=q1)    # H til
    ytil[jseq] <- y1 - fxtil + Wtil[jseq,]%*%beta1 + Htil[jseq,]%*%ub1
    centerbi[jseq] <- Htil[jseq,]%*%(matrix(ub1,ncol=1)-c.*Deltab)
  }
  centerpar <- fittedVal -centerbi
  
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalNL,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtNL,nu=nu,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsNL,nu=nu,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnNL,nu=nu,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  lv
}


logveroNL2 = function(y,#response vector (N x 1)
                     x,#covariates matrix (N x p1)
                     ind,#factor whose levels indicates the subjects or groups (N x 1) 
                     nlfder,#first derivative of the nonlinear function 
                     fitObj,#object returned from EM.Skew.NL or EM.sim.NL function
                     distr,#distribution to be used: "sn" (normal), "st" (student's-t), "ss" (slash), or "scn" (contaminated normal)
                     theta=fitObj$estimates#fitObj$estimates
){ 
  fittedVal=fitObj$fitted
  biest=fitObj$random.effects
  beta1=theta$beta
  sigmae=theta$sigma2
  D1=theta$D
  lambda=theta$lambda
  nu=theta$nu
  #
  if (is.null(lambda)) lambda=rep(0,nrow(D1))
  
  m<-n_distinct(ind)
  N<-length(ind)
  p1<-dim(x)[2]
  p<-length(beta1)
  q1<- nrow(D1)
  
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda));
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  #
  if (distr=="sn") {c. = -sqrt(2/pi)}
  if (distr=="st") {c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  
  ###
  Wtil<-matrix(nrow=N,ncol=p)
  Htil<-matrix(nrow=N,ncol=q1)
  ytil<-numeric(N)
  centerbi <- numeric(N)
  for (j in 1:nlevels(ind)) {
    jseq = ind==levels(ind)[j]
    y1=y[jseq]
    x1=matrix(x[jseq,  ],ncol=p1)
    ub1 <- biest[j,]#ub0[(((j-1)*q1)+1) : (j*q1), j]
    nj = length(y1)
    HWmat <- matrix(0,nrow=nj,ncol=(p+q1))
    fxtil <- matrix(0,nrow=nj,ncol=1)
    for(i in 1:nj)
    {
      farg <- as.list(c(x1[i,],beta1,ub1))
      names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep=""))
      formals(nlfder) <- farg
      fx <- nlfder()
      fxtil[i,] <- fx[1]
      HWmat[i,] <- attr(fx,"gradient")
    }
    
    Wtil[jseq,] <- matrix(HWmat[,1:p],nrow=nj,ncol=p)            # W til
    Htil[jseq,] <- matrix(HWmat[,(p+1):(p+q1)],nrow=nj,ncol=q1)    # H til
    ytil[jseq] <- y1 - fxtil + Wtil[jseq,]%*%beta1 + Htil[jseq,]%*%ub1
    centerbi[jseq] <- Htil[jseq,]%*%(matrix(ub1,ncol=1)-c.*Deltab)
  }
  centerpar <- fittedVal -centerbi
  
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormalNL,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljtNL,nu=nu,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljsNL,nu=nu,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcnNL,nu=nu,y=y,centerpar=centerpar,z=Htil,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  lv
}
###############################################################
##### auxiliary functions
###############################################################

################################################################
#Log-likelihood NL
################################################################
ljnormalNL <-function(j,y,centerpar,z,Gammab,Deltab,sigmae){
  y1=y[j]
  q1=ncol(Gammab)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-centerpar[j]#fiti - z1%*%(biesti- c.*Deltab)
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj) #z1 D1 z1^t + sig2*I_nj ??
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljtNL <-function(j,nu,y,centerpar,z,Gammab,Deltab,sigmae){
  y1=y[j]
  q1=ncol(Gammab)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-centerpar[j]#fiti - z1%*%(biesti- c.*Deltab)
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  #log(2*dmvt(y1,delta = med, sigma = Psi, df = nu,log=F)*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))#veroST1(Psi,Ajj,dj,nu,pp=njj))
  log(2*dtj*pt(sqrt(nu+njj)*Ajj/sqrt(dj+nu),nu+njj))
}
#
ljsNL <-function(j,nu,y,centerpar,z,Gammab,Deltab,sigmae){
  y1=y[j]
  q1=ncol(Gammab)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-centerpar[j]#fiti - z1%*%(biesti- c.*Deltab)
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-as.numeric(sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med))
  #f <- function(u) u^(nu - 1)*dmvnorm(y1,med,Psi/u)*pnorm(u^(1/2)*Ajj)
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
  resp <- integrate(Vectorize(f2),0,1)$value
  log(2*nu*resp)
}
#
ljcnNL <-function(j,nu,y,centerpar,z,Gammab,Deltab,sigmae){
  y1=y[j]
  q1=ncol(Gammab)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-centerpar[j]#fiti - z1%*%(biesti- c.*Deltab)
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}
