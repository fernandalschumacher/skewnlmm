####################################################################################################
## SMSN Nonlinear Mixed Effects Model ##
####################################################################################################
## last updated in 2020-07-26
## downloaded from https://github.com/fernandalschumacher/skewnlmm

#function for estimating a SMSN-NLME model

EM.Skew.NL<- function(y,#response vector (N x 1)
                      x,#covariates matrix (N x p1)
                      ind,#factor whose levels indicates the subjects or groups (N x 1) 
                      nlfder,#first derivative of the nonlinear function 
                      nlf,#nonlinear function with args: x1,...,xp1, beta1,...,betap, b1,...,bq
                      distr,#distribution to be used: "sn" (skew-normal), "st" (skew-student's-t), "ss" (skew-slash), or "scn" (skew-contaminated normal)
                      beta1,#p x 1 vector with initial values for beta
                      sigmae,#initial value for sigma2
                      D1,#q x q matrix with initial values for D 
                      lambda,#q x 1 vector with initial values for the skewness parameter
                      nu=NULL,#initial values for nu, when distr != "sn"
                      lb=NULL,#lower bound to be used inside optim for updating nu, when distr != "sn"
                      lu=NULL,#upper bound to be used inside optim for updating nu, when distr != "sn"
                      bi=NULL,#optional. Matrix of initial random effects (nlevels(ind) x q) 
                      precisao=1e-4,#tolarance for the convergence criterion
                      max.iter=600,#maximum number of iterations for the EM algorithm
                      showiter=T,#logical. Should the iteration message be printed?
                      showerroriter=F,#logical. Should some convergence information be printed?
                      informa=T,#logical. Should standard errors be computed?
                      estimatenu=T#logical. Should nu be estimated? If FALSE, the value passed in nu is kept fixed
                      ){
  ti = Sys.time()
  #
  m<-n_distinct(ind)
  N<-length(ind)
  p1<-ncol(x)
  p <- length(beta1)
  q1 <- dim(D1)[2]
  #
  if (is.null(bi)) bi <- matrix(0,ncol=q1,nrow=m)
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  
  if (estimatenu) {
    teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,nu)
  } else teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab)
  
  criterio<-10
  count<-0
  llji=1
  #
  loglikVec<- numeric(max.iter)
  while((criterio > precisao)&(count<max.iter)){
    
    count <- count + 1
    # linearization step
    Wtil<-matrix(nrow=N,ncol=p)
    Htil<-matrix(nrow=N,ncol=q1)
    ytil<-numeric(N)
    for (j in 1:nlevels(ind)) {
      jseq = ind==levels(ind)[j]
      y1=y[jseq]
      x1=matrix(x[jseq,  ],ncol=p1)
      ub1 <- matrix(bi[j,],ncol=1)
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
    }
    
    # expectation step
    res_emj = revert_list(tapply(1:N,ind,emj,y=ytil, x=Wtil, z=Htil, beta1=beta1, Gammab=Gammab,
                                 Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    sum5 = Reduce("+",res_emj$sum5)
    ut2j = unlist(res_emj$ut2j,use.names = F)
    uj = unlist(res_emj$uj,use.names = F)
    bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
    
    # maximization step
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    Gammab<-sum4/m
    Deltab<-sum5/sum(ut2j)
    #
    D1<-Gammab+Deltab%*%t(Deltab)
    lambda<-matrix.sqrt(solve(D1))%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%solve(D1)%*%Deltab))
    #
    param <- teta
    #
    if (estimatenu) {
      logvero1<-function(nu){logvero(ytil, Wtil, Htil, ind, beta1, sigmae, D1, lambda, distr, nu)}
  
      if (distr=="sn"){ nu<-NULL} else
      {
       nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
      }
      teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab,nu)
    } else teta <- c(beta1,sigmae,Gammab[upper.tri(Gammab, diag = T)],Deltab)
    loglikVec[count] <- logvero(ytil, Wtil, Htil, ind, beta1, sigmae, D1, lambda, distr, nu)
    if (count>9){ #only updates after 10 iterations
      at<- (loglikVec[count]-loglikVec[count-1])/(loglikVec[count-1]-loglikVec[count-2])
      criterio<-abs((loglikVec[count]-loglikVec[count-1])/(1-at))
      #print(loglik[count])
    }
    if (is.nan(criterio)) criterio=10
    if (all(is.nan(teta))) stop("NaN values")
    if (showiter&!showerroriter) cat("Iteration ",count," of ",max.iter,"\r") 
    if (showerroriter) cat("Iteration ",count," of ",max.iter," - criterium =",criterio," - loglik =",loglikVec[count],"\r") #  criterium ",criterio," or ",criterio2,"\r")
  }
  if (count==max.iter) message("\n maximum number of iterations reachead")
  cat("\n")
  ####update bi and ui
  Wtil<-matrix(nrow=N,ncol=p)
  Htil<-matrix(nrow=N,ncol=q1)
  ytil<-numeric(N)
  for (j in 1:nlevels(ind)) {
    jseq = ind==levels(ind)[j]
    y1=y[jseq]
    x1=matrix(x[jseq,  ],ncol=p1)
    ub1 <- matrix(bi[j,],ncol=1)
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
  }
  
  #
  res_emj = revert_list(tapply(1:N,ind,emj,y=ytil, x=Wtil, z=Htil, beta1=beta1, Gammab=Gammab,
                               Deltab=Deltab, sigmae=sigmae, distr=distr,nu=nu))
  uj = unlist(res_emj$uj,use.names = F)
  bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
  ####
  # creating object to return
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,dd,lambda,nu)
  if (distr=="sn") names(theta)<-c( paste0("beta",1:p),"sigma2",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1))
  else names(theta)<- c( paste0("beta",1:p),"sigma2",paste0("Dsqrt",1:length(dd)),paste0("lambda",1:q1),paste0("nu",1:length(nu)))
  
  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           D=D1,lambda=as.numeric(lambda)),
                  uhat=unlist(res_emj$uj))
  if (distr != "sn") obj.out$estimates$nu = nu
  obj.out$random.effects<- bi
  
  obj.out$loglik <-loglikVec[count]
  # fitted values
  fittedval<- numeric(length(ind))
  ind_levels <- levels(ind)
  for (i in seq_along(ind_levels)) {
    seqi <- which(ind==ind_levels[i])
    xfiti <- matrix(x[seqi,],ncol=ncol(x))
    ub1 <- matrix(bi[i,],ncol=1)
    for (j in seq_along(seqi)) 
      {
      farg <- as.list(c(xfiti[j,],beta1,ub1))
      names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep=""))
      formals(nlf) <- farg
      fittedval[seqi[j]] <- nlf()
    }
  }
  obj.out$fitted <- fittedval
  #
  if (informa) {
    desvios<-try(
      Infmatrix(ytil,Wtil,Htil,ind,beta1,sigmae,D1,lambda,distr = distr,nu = nu)
      ,silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      desvios[substr(names(desvios),1,1)=="l"]=NA
      obj.out$std.error=desvios
    }
  }
  
  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}

######

# function to make prediction from a SMSN-NLME model
# for SMN-NLME predictions, set lambda=rep(0,q)

predictf.skew.NL<- function(yfit,#response vector from fitted data
                            xfit,#covariates matrix from fitted data
                            indfit,#subjects or groups factor from fitted data
                            xpred,#covariates matrix from data to be predicted
                            indpred,#subjects or groups factor from data to be predicted
                            distr,#distribution to be used: "sn", "st", "ss", or "scn"
                            nlfder,#first derivative of the nonlinear function 
                            nlf,#nonlinear function with args: x1,...,xp1, beta1,...,betap, b1,...,bq
                            theta,#fitted parameters c(beta,sigma,dd,lambda,nu)
                            p,#lenght(beta)
                            estbi #estimated random effects
                            ){
  p1 <- ncol(xfit)
  q1 <- ncol(estbi)
  q2 <- q1*(q1+1)/2
  beta1 <- matrix(theta[1:p],ncol=1)
  sigmae <- as.numeric(theta[p+1])
  dd <- theta[(p+2):(p+1+q2)]
  lambda <- matrix(theta[(p+q2+2):(p+1+q2+q1)],ncol=1)
  if (distr=="sn") {nu<- NULL; c. = -sqrt(2/pi)}
  if (distr=="st") {nu<- theta[p+q2+q1+2]; c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)}
  if (distr=="ss") {nu<- theta[p+q2+q1+2]; c.=-sqrt(2/pi)*nu/(nu-.5)}
  if (distr=="scn") {nu<- theta[(p+q2+q1+2):(p+q2+q1+3)]; c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))}
  if ((p+1+q2+q1+length(nu))!=length(theta)) stop("theta misspecified")
  D1sqrt <- Dmatrix(dd)
  D1 <- D1sqrt%*%D1sqrt
  #
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-D1sqrt%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  ypred <- numeric(length = length(indpred))
  indpred <- droplevels(indpred)
  #
  for (indj in levels(indpred)) {
    y1 <- yfit[indfit==indj]
    xfit1 <- matrix(xfit[indfit==indj,],ncol=p1)
    xpred1 <- matrix(xpred[indpred==indj,],ncol=p1)
    xPlus <- rbind(xfit1,xpred1)
    njFit = nrow(xfit1)
    njPred = nrow(xpred1)
    nj = nrow(xPlus)
    seqFit = 1:njFit
    seqPred = njFit+1:njPred
    # linearizacao
    ub1 = estbi[row.names(estbi)==indj,]
    HWmat <- matrix(0,nrow=nj,ncol=(p+q1))
    fxtil <- fxnltil <- matrix(0,nrow=nj,ncol=1)
    for(i in 1:nj){
      farg <- as.list(c(xPlus[i,],beta1,ub1))
      names(farg) = c(paste("x",1:p1,sep=""),paste("beta",1:p,sep=""),paste("b",1:q1,sep=""))
      formals(nlfder) <- farg
      fx <- nlfder()
      fxtil[i,] <- fx[1]
      HWmat[i,] <- attr(fx,"gradient")
      formals(nlf) <- farg
      fxnltil[i,] <- nlf()
    }
    Wtil <- matrix(HWmat[,1:p],nrow=nj,ncol=p)            # W til
    Htil <- matrix(HWmat[,(p+1):(p+q1)],nrow=nj,ncol=q1)    # H til
    #
    xPlus1 <- Wtil#model.matrix(formFixed,data=dataPlus)
    zPlus1<-Htil#model.matrix(formRandom,data=dataPlus)
    z1 <- matrix(zPlus1[seqFit,],ncol=ncol(zPlus1))
    x1 <- matrix(xPlus1[seqFit,],ncol=ncol(xPlus1))
    z1Pred <- matrix(zPlus1[seqPred,],ncol=ncol(zPlus1))
    x1Pred <- matrix(xPlus1[seqPred,],ncol=ncol(xPlus1))
    #
    medFit <- fxnltil[seqFit]-z1%*%ub1 +c.*z1%*%Deltab
    medPred <- fxnltil[seqPred] -z1Pred%*%ub1 + c.*z1Pred%*%Deltab
    #
    SigmaPlus = sigmae*diag(nj)
    PsiPlus<-(zPlus1)%*%(D1)%*%t(zPlus1)+SigmaPlus
    dj<-as.numeric(t(y1-medFit)%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit))
    LambdaPlus <- solve(solve(D1)+ t(zPlus1)%*%solve(SigmaPlus)%*%zPlus1)
    sPsiPlus <- solve(PsiPlus)
    lambdaBarPlus <- matrix.sqrt(sPsiPlus)%*%zPlus1%*%D1%*%zeta/as.numeric(sqrt(1+t(zeta)%*%LambdaPlus%*%zeta))
    vj <- matrix.sqrt(sPsiPlus)%*%lambdaBarPlus
    Psi22.1 <- PsiPlus[seqPred,seqPred]- PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]
    Psi22.1 <- (Psi22.1+t(Psi22.1))/2
    vjtil1<-(vj[seqFit] + solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]%*%vj[seqPred])
    vjtil <- (vj[seqFit] + solve(PsiPlus[seqFit,seqFit])%*%PsiPlus[seqFit,seqPred]%*%vj[seqPred])/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))
    tau2.1 <- as.numeric(t(vjtil1)%*%(y1-medFit))
    Ajj<-as.numeric(t(vjtil)%*%(y1-medFit)) 
    mu2.1 <- medPred + PsiPlus[seqPred,seqFit]%*%solve(PsiPlus[seqFit,seqFit])%*%(y1-medFit)
    if (distr=="sn") tau1 <- dnorm(Ajj)/pnorm(Ajj)
    if (distr=="st") {
      dtj = gamma((nu+njFit)/2)/gamma(nu/2)/pi^(njFit/2)/sqrt(det(PsiPlus[seqFit,seqFit]))*nu^(-njFit/2)*
        (dj/nu+1)^(-(njFit+nu)/2)
      denST = 2*dtj*pt(sqrt(nu+njFit)*Ajj/sqrt(dj+nu),nu+njFit)
      tau1 <- as.numeric(dmvt(y1,delta=medFit,sigma=PsiPlus[seqFit,seqFit],df=nu,log=F)*gamma((nu+njFit-1)/2)*(nu+dj)^((nu+njFit)/2)/
                           (denST*pi^.5*gamma((nu+njFit)/2)*(nu+dj+Ajj^2)^((nu+njFit-1)/2)))
    }
    if (distr=="ss") {
      f2 <- function(u) u^(nu - 1)*((2*pi)^(-njFit/2))*(u^(njFit/2))*((det(PsiPlus[seqFit,seqFit]))^(-1/2))*
        exp(-0.5*u*dj)*pnorm(u^(1/2)*Ajj)
      denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
      tau1<-2^(nu)*nu*gamma(nu-.5+njFit/2)*pgamma(1,nu-.5+njFit/2,(dj+Ajj^2)/2)/
        (denSS*(dj+Ajj^2)^(nu-.5+njFit/2)*pi^(njFit/2+.5)*det(PsiPlus[seqFit,seqFit])^.5)
    }
    if (distr=="scn") {
      fy<-as.numeric(2*(nu[1]*dmvnorm(y1,medFit,(PsiPlus[seqFit,seqFit]/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                          (1-nu[1])*dmvnorm(y1,medFit,PsiPlus[seqFit,seqFit])*pnorm(Ajj,0,1)))
      tau1<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,medFit,PsiPlus[seqFit,seqFit]/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,medFit,PsiPlus[seqFit,seqFit])*dnorm(Ajj,0,1))/fy)
    }
    ypredj <- mu2.1 + tau1/as.numeric(sqrt(1+t(vj[seqPred])%*%Psi22.1%*%vj[seqPred]))*Psi22.1%*%vj[seqPred]
    ypred[indpred==indj] <- ypredj
  }
  data.frame(groupVar=indpred,xpred,ypred)
}


## Healy's plot
healy.plot <- function(object,#fitted object from EM.Skew.NL or EM.sim.NL functions
                       y,x,ind,nlfder,distr,dotsize=0.4,...) {
  ind_levels <- levels(ind)
  njvec <- tapply(ind,ind,length)
  nj1 <- as.numeric(njvec[1])
  N <- length(ind)
  p<- length(object$estimates$beta)
  q1<- nrow(object$estimates$D)
  p1 <- ncol(x)
  if (!all(njvec==nj1)) stop("for this plot all subjects must have the same number of observations, you can predict the missing values and enter the full dataset using the dataPlus argument")
  #
  nu <- object$estimates$nu
  #
  if (distr=="sn") {c.=-sqrt(2/pi)}
  if (distr=="st") {c.=-sqrt(object$estimates$nu/pi)*
    gamma((object$estimates$nu-1)/2)/gamma(object$estimates$nu/2)}
  if (distr=="ss") {c.=-sqrt(2/pi)*object$estimates$nu/(object$estimates$nu-.5)}
  if (distr=="scn") {c.=-sqrt(2/pi)*(1+object$estimates$nu[1]*
                                       (object$estimates$nu[2]^(-.5)-1))}
  #
  beta1<- object$estimates$beta
  Dest <- object$estimates$D
  if (is.null(object$estimates$lambda)) object$estimates$lambda=rep(0,q1) 
  delta = object$estimates$lambda/as.numeric(
    sqrt(1+t(object$estimates$lambda)%*%(object$estimates$lambda)))
  Deltab = matrix.sqrt(Dest)%*%delta
  #
  #computing approx. from Theorem 1
  Htil<-matrix(nrow=N,ncol=q1)
  centerbi <- numeric(N)
  for (j in 1:nlevels(ind)) {
    jseq = ind==levels(ind)[j]
    y1=y[jseq]
    x1=matrix(x[jseq,  ],ncol=p1)
    ub1 <- object$random.effects[j,]
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
    Htil[jseq,] <- matrix(HWmat[,(p+1):(p+q1)],nrow=nj,ncol=q1)   
    centerbi[jseq] <- Htil[jseq,]%*%(matrix(ub1,ncol=1)-c.*Deltab)
  }
  #
  centerpar <- object$fitted -centerbi
  # computing Mahalanobis distance
  mahaldist <- numeric(length(ind_levels))
  for (i in seq_along(ind_levels)) {
    seqi <- ind==ind_levels[i]
    zfiti <- Htil[seqi,]
    Sigmaest <- object$estimates$sigma2*diag(sum(seqi))
    Psiy <- Sigmaest+ (zfiti)%*%Dest%*%t(zfiti)
    #
    ztil <- (y-centerpar)[seqi]
    mahaldist[i] <- t(ztil)%*%solve(Psiy)%*%ztil
  }
  # generating a Healy-type plot
  mahaldist <- sort(mahaldist)
  #
  m <- length(mahaldist)
  empprob <- (1:m)/(m)
  #
  if (distr=="sn") theoprob <- pchisq(mahaldist,nj1)
  if (distr=="st") theoprob <- pf(mahaldist/nj1,nj1,nu)
  if (distr=="ss") {
    theoprob <- pchisq(mahaldist,nj1) - 2^nu*gamma(nj1/2+nu)/(mahaldist^nu)/gamma(nj1/2)*
      pchisq(mahaldist,nj1+2*nu)
  }
  if (distr=="scn") {
    theoprob <- nu[1]*pchisq(nu[2]*mahaldist,nj1) + (1-nu[1])*pchisq(mahaldist,nj1)
  }
  ####
  #print(cor(empprob,theoprob))
  qplot(empprob,theoprob,size=I(dotsize),...=...) + geom_abline(intercept=0,slope=1) +
    ylab(NULL) + xlab(NULL) + theme_minimal() + 
    theme(plot.title = element_text( face="italic", size=10))+
    expand_limits(y = c(0, 1),x=c(0,1))
}

###############################################################
##### auxiliary functions
###############################################################

is.wholenumber <- function(x, tol1 = .Machine$double.eps^0.5)  abs(x - round(x)) < tol1
################################################################
#Root of a symmetric matrix
################################################################
matrix.sqrt <- function(A)
{
  if (length(A)==1) return(sqrt(A))
  else{
    sva <- svd(A)
    if (min(sva$d)>=0) {
      Asqrt <- sva$u%*%diag(sqrt(sva$d))%*%t(sva$v) # svd 
      if (all(abs(Asqrt%*%Asqrt-A)<1e-4)) return(Asqrt)
      else stop("Matrix square root is not defined/not real")
    }
    else stop("Matrix square root is not defined/not real")
  }
}
################################################################
#trace of a matrix of dim >=1
################################################################
traceM <- function(Mat){
  if(length(Mat)==1) tr<- as.numeric(Mat)
  else tr<-sum(diag(Mat))
  tr
}
################################################################
#inverter ordem de hierarquia de uma lista com nomes
################################################################
revert_list <- function(ls) { # @Josh O'Brien
  # get sub-elements in same order
  x <- lapply(ls, `[`, names(ls[[1]]))
  # stack and reslice
  apply(do.call(rbind, x), 2, as.list)
}

Dmatrix <- function(dd) {#function to construct D matrix based on distinct elements vector dd
  q2 <- length(dd)
  q1 <- -.5+sqrt(1+8*q2)/2
  if (q1%%1 != 0) stop("wrong dimension of dd")
  D1 <- matrix(nrow = q1,ncol = q1)
  D1[upper.tri(D1,diag = T)] <- as.numeric(dd)
  D1[lower.tri(D1)] <- t(D1)[lower.tri(D1)]
  return(D1)
}


################################################################
#Log-likelihood - independent
################################################################
ljnormal <-function(j,y,x,z,beta1,Gammab,Deltab,sigmae){
  c. = -sqrt(2/pi)
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj) #z1 D1 z1^t + sig2*I_nj ??
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))
}
#
ljt <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
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
ljs <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(2/pi)*nu/(nu-.5)
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
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
ljcn <-function(j,nu,y,x,z,beta1,Gammab,Deltab,sigmae){
  c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  njj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  Ajj<-sqrt(Mtj2)*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(njj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  log(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(sqrt(nu[2])*Ajj,0,1)+
           (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
}
logvero = function(y,x,z,ind,beta1,sigmae,D1,lambda,distr,nu){ #ind = indicadora de individuo
  m<-n_distinct(ind)
  N<-length(ind)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda));
  Deltab<-matrix.sqrt(D1)%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormal,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljt,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljs,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcn,nu=nu,y=y,x=x,z=z,beta1=beta1,Gammab=Gammab,Deltab=Deltab,sigmae=sigmae))
  lv
}

##############################################################################
# EM - independent
##############################################################################
emj = function(jseq, y, x, z, beta1, Gammab, Deltab, sigmae,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  #
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  nj = length(y1)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+sigmae*diag(nj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(sigmae*diag(nj)+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(sigmae*diag(nj)+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ajj<-as.numeric(mutj/sqrt(Mtj2))
  D1<- Gammab+Deltab%*%t(Deltab)
  #
  mediab<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)+c.*Deltab
  Lambda<-solve(solve(D1)+t(z1)%*%z1/sigmae)
  lambda<-matrix.sqrt(solve(D1))%*%Deltab/as.numeric(sqrt(1-t(Deltab)%*%solve(D1)%*%Deltab))
  zeta<-matrix.sqrt(solve(D1))%*%lambda
  #
  if  (distr=="sn"){
    uj<-1
    esper<-as.numeric(dnorm(Ajj,0,1)/pnorm(Ajj,0,1))
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*as.numeric(dnorm(Ajj,0,1))/as.numeric(pnorm(Ajj,0,1))
  }
  
  if (distr=="st"){
    dtj = gamma((nu+nj)/2)/gamma(nu/2)/pi^(nj/2)/sqrt(det(Psi))*nu^(-nj/2)*(dj/nu+1)^(-(nj+nu)/2)
    denST = 2*dtj*pt(sqrt(nu+nj)*Ajj/sqrt(dj+nu),nu+nj)
    uj<-as.numeric(2^2*(nu+dj)^(-(nu+nj+2)/2)*gamma((nj+nu+2)/2)*nu^(nu/2)*
                     pt(sqrt((nj+nu+2)/(dj+nu))*Ajj,nj+nu+2)/(gamma(nu/2)*det(Psi)^.5*denST*pi^(nj/2)))
    esper <-as.numeric(2*nu^(nu/2)*gamma((nj+nu+1)/2)/(denST*pi^((nj+1)/2)*gamma(nu/2)*det(Psi)^.5*(nu+dj+Ajj^2)^((nu+nj+1)/2)))
    esper2<- as.numeric(dmvt(y1,delta=med,sigma=Psi,df=nu,log=F)*gamma((nu+nj-1)/2)*(nu+dj)^((nu+nj)/2)/
                          (denST*pi^.5*gamma((nu+nj)/2)*(nu+dj+Ajj^2)^((nu+nj-1)/2)))
    bi<-mediab+Lambda%*%zeta/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))*esper2
  }
  
  if (distr=="ss"){
    f2esp <- function(s) pnorm(s^.5*Ajj)*dgamma(s,nu+1+nj/2,dj/2)/pgamma(1,nu+1+nj/2,dj/2)
    EspVal <- integrate(f2esp,0,1)$value#mean(pnorm(S^(1/2)*Ajj))#
    f2 <- function(u) u^(nu - 1)*((2*pi)^(-nj/2))*(u^(nj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))*pnorm(u^(1/2)*Ajj)
    denSS <- 2*nu*integrate(Vectorize(f2),0,1)$value
    uj<-2^(2+nu)*nu*gamma(nu+1+nj/2)*pgamma(1,nu+1+nj/2,dj/2)*EspVal/
      (denSS*dj^(nu+1+nj/2)*pi^(nj/2)*det(Psi)^.5)
    esper <- 2^(1+nu)*nu*gamma(nu+.5+nj/2)*pgamma(1,nu+.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu+.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
    esper2<-2^(nu)*nu*gamma(nu-.5+nj/2)*pgamma(1,nu-.5+nj/2,(dj+Ajj^2)/2)/
      (denSS*(dj+Ajj^2)^(nu-.5+nj/2)*pi^(nj/2+.5)*det(Psi)^.5)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }
  
  if (distr=="scn"){
    fy<-as.numeric(2*(nu[1]*dmvnorm(y1,med,(Psi/nu[2]))*pnorm(nu[2]^(1/2)*Ajj,0,1)+
                        (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1)))
    uj<-2*(nu[1]*nu[2]*dmvnorm(y1,med,Psi/nu[2])*pnorm(nu[2]^(1/2)*Ajj,0,1)+
             (1-nu[1])*dmvnorm(y1,med,Psi)*pnorm(Ajj,0,1))/fy
    esper<-as.numeric(2*(nu[1]*nu[2]^(1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                           (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
    esper2<-as.numeric(2*(nu[1]*nu[2]^(-1/2)*dmvnorm(y1,med,Psi/nu[2])*dnorm(nu[2]^(1/2)*Ajj,0,1)+
                            (1-nu[1])*dmvnorm(y1,med,Psi)*dnorm(Ajj,0,1))/fy)
    bi<-mediab+Lambda%*%zeta*esper2/as.numeric(sqrt(1+t(zeta)%*%Lambda%*%zeta))
  }
  
  Tbj<-solve(solve(Gammab)+t(z1)%*%z1/sigmae)
  r<-Tbj%*%t(z1)%*%(y1-x1%*%beta1)/sigmae
  s1<-(diag(q1)-Tbj%*%t(z1)%*%z1/sigmae)%*%Deltab
  utj<-as.numeric(uj*(mutj+c.)+sqrt(Mtj2)*esper)
  ut2j<-as.numeric(uj*(mutj+c.)^2+Mtj2+sqrt(Mtj2)*(mutj+2*c.)*esper)
  ub<-uj*r+s1*utj
  utbj<- r*utj+s1*ut2j
  ub2j<-Tbj+uj*r%*%t(r)+s1%*%t(r)*utj+r%*%t(s1)*utj+s1%*%t(s1)*ut2j
  #
  sum1<-uj*t(x1)%*%x1 #denom beta
  sum2<-(t(x1)%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%z1%*%ub-
    t(ub)%*%t(z1)%*%(y1-x1%*%beta1)+traceM(ub2j%*%t(z1)%*%z1) #soma do sig2
  sum4<-ub2j-utbj%*%t(Deltab)-Deltab%*%t(utbj)+
    ut2j*Deltab%*%t(Deltab) #soma do Gamma
  sum5<-utbj #num do delta
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,sum5=sum5,ut2j=ut2j,uj=uj)
  obj.out$bi=bi
  return(obj.out)
}


#######

#Information matrix for SMSN-LMM and SMSN-LMM-AR(p) with E(bi)=0
IPhi <- function(w,di,Ai,distr,nu) {
  if (distr=="sn") intval <- exp(-.5*di)*pnorm(Ai)
  else if (distr=="st") intval<- 2^w*nu^(nu/2)*gamma(w+nu/2)*
      pt(Ai*sqrt((2*w+nu)/(di+nu)),2*w+nu)/((di+nu)^(w+nu/2)*gamma(nu/2))
  else if (distr=="ss") {
    U <- runif(5000)
    V <- pgamma(1,w+nu, di/2)*U
    S <- qgamma(V,w+nu, di/2)
    intval<- nu*2^(nu+w)*gamma(nu+w)*pgamma(1,nu+w,di/2)*mean(pnorm(S^(1/2)*Ai))/(di^(nu+w))
  }
  else if (distr=="scn") intval <- nu[1]*nu[2]^w*exp(-.5*nu[2]*di)*pnorm(nu[2]^.5*Ai)+
      (1-nu[1])*exp(-.5*di)*pnorm(Ai)
  return(intval)
}
Iphi <- function(w,di,Ai,distr,nu) {
  if (distr=="sn") intval <- exp(-.5*di)*dnorm(Ai)
  else if (distr=="st") intval<- 2^w*nu^(nu/2)*gamma(w+nu/2)/
      (sqrt(2*pi)*gamma(nu/2)*(di+Ai^2+nu)^(w+nu/2))
  else if (distr=="ss") intval<- nu*2^(nu+w)*gamma(nu+w)*pgamma(1,nu+w,(di+Ai^2)/2)/
      (sqrt(2*pi)*(di+Ai^2)^(nu+w))
  else if (distr=="scn") intval <- nu[1]*nu[2]^w*exp(-.5*nu[2]*di)*dnorm(nu[2]^.5*Ai)+
      (1-nu[1])*exp(-.5*di)*dnorm(Ai)
  return(intval)
}

F.r <- function(r,q1){
  Fmat. <- matrix(0,ncol=q1,nrow=q1)
  Fmat.[upper.tri(Fmat.,diag=T)][r] = 1
  Fmat.[lower.tri(Fmat.)] = t(Fmat.)[lower.tri(Fmat.)]
  Fmat.
}

autocovsAR <- function(phi,n) {
  p <- length(phi)
  if (n==1) Rn <- 1
  else Rn<- ARMAacf(ar=phi, ma=0, lag.max = n-1)
  rhos <- ARMAacf(ar=phi, ma=0, lag.max = p)[-1]
  return(Rn/(1-sum(rhos*phi)))
}

scorei <- function(jseq,y,x,z,beta1,sigmae,D1,lambda,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  q2 = q1*(q1+1)/2
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  ni = length(y1)
  Fmat = matrix.sqrt(D1)
  delta<-lambda/as.numeric(sqrt(1+t(lambda)%*%lambda))
  Deltab<-Fmat%*%delta
  Gammab<-D1-Deltab%*%t(Deltab)
  med<-x1%*%beta1+ c.*z1%*%Deltab
  Sigma <- sigmae*diag(ni)
  Psi<-(z1)%*%(Gammab+Deltab%*%t(Deltab))%*%t(z1)+Sigma
  sPsi <- solve(Psi)
  di<-as.numeric(t(y1-med)%*%sPsi%*%(y1-med))
  Mtj2<-(1+t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%z1%*%Deltab)^(-1)
  mutj<-Mtj2*t(Deltab)%*%t(z1)%*%solve(Sigma+z1%*%Gammab%*%t(z1))%*%(y1-med)
  Ai<-as.numeric(mutj/sqrt(Mtj2))
  sFmat = solve(Fmat)
  Lambda = solve(solve(D1)+ t(z1)%*%z1/sigmae)
  F.lista <- lapply(1:q2,F.r,q1=q1)
  #theta = c(beta1,sigmae,dd,lambda,nu) - para independente
  indpar = c(rep("beta",p),"sigma",rep("dd",q2),rep("lambda",q1))
  lpar = length(indpar)
  ##### derivadas de log(det(Psi))
  dlogdpsi = numeric(lpar)
  dlogdpsi[indpar=="sigma"] =traceM(sPsi)
  for (i in 1:q2) dlogdpsi[indpar=="dd"][i] = traceM(sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+
                                                                    Fmat%*%F.lista[[i]])%*%t(z1))
  
  ##### derivadas de Ai
  dAi = numeric(lpar)
  ai = as.numeric((1+t(lambda)%*%sFmat%*%Lambda%*%sFmat%*%lambda)^.5)
  bi = as.numeric((1+t(lambda)%*%lambda)^.5)
  Bi = as.numeric(t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%Fmat%*%lambda)
  dAi[indpar=="beta"] = -1/ai*t(x1)%*%sPsi%*%z1%*%Fmat%*%lambda
  dAi[indpar=="lambda"] = 1/ai*Fmat%*%t(z1)%*%sPsi%*%(y1-x1%*%beta1- 2*c.*z1%*%Deltab)-
    1/ai^2*Ai*sFmat%*%Lambda%*%sFmat%*%lambda + c.*Bi/ai/(bi^3)*lambda
  dAi[indpar=="sigma"] = -1/ai*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%sPsi%*%(y1-med)-
    Ai/(2*ai^2*sigmae^2)*t(lambda)%*%sFmat%*%Lambda%*%t(z1)%*%z1%*%Lambda%*%sFmat%*%lambda
  for (i in 1:q2) dAi[indpar=="dd"][i] = 1/ai*(t(lambda)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med)-
                                                 t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med) -
                                                 c./bi*t(lambda)%*%Fmat%*%t(z1)%*%sPsi%*%z1%*%F.lista[[i]]%*%lambda)+
    1/ai^2*Ai/2*t(lambda)%*%sFmat%*%(F.lista[[i]]%*%sFmat%*%Lambda+Lambda%*%sFmat%*%F.lista[[i]]-Lambda%*%sFmat%*%(F.lista[[i]]%*%sFmat+sFmat%*%F.lista[[i]])%*%sFmat%*%Lambda)%*%sFmat%*%lambda
  
  ##### derivadas de di
  ddi = numeric(lpar)
  ddi[indpar=="beta"] =-2*t(x1)%*%sPsi%*%(y1-med)
  ddi[indpar=="lambda"] = -2*c./bi*(Fmat-delta%*%t(Deltab))%*%t(z1)%*%sPsi%*%(y1-med)
  ddi[indpar=="sigma"] = -t(y1-med)%*%sPsi%*%sPsi%*%(y1-med)
  for (i in 1:q2) ddi[indpar=="dd"][i] = -2*c.*t(delta)%*%F.lista[[i]]%*%t(z1)%*%sPsi%*%(y1-med) -
    t(y1-med)%*%sPsi%*%z1%*%(F.lista[[i]]%*%Fmat+Fmat%*%F.lista[[i]])%*%t(z1)%*%sPsi%*%(y1-med)
  
  ##### derivadas de ki
  ki = IPhi(ni/2,di=di,Ai=Ai,distr = distr,nu=nu)
  dki = -.5*IPhi(ni/2+1,di=di,Ai=Ai,distr = distr,nu=nu)*ddi+
    Iphi(ni/2+.5,di=di,Ai=Ai,distr = distr,nu=nu)*dAi
  
  sihat = -.5*dlogdpsi+1/ki*dki
  sihat
  
}


Infmatrix <- function(y,x,z,ind,beta1,sigmae,D1,lambda,distr,nu){
  N <-length(y)
  score_list=tapply(1:N,ind,scorei,y=y, x=x, z=z, beta1=beta1, sigmae=sigmae,D1=D1,lambda=lambda,distr=distr,nu=nu)
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt,ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-20*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

