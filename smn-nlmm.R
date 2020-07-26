####################################################################################################
## SMN Nonlinear Mixed Effects Model ##
####################################################################################################
## last updated in 2020-07-26
## downloaded from https://github.com/fernandalschumacher/skewnlmm

#function for estimating a SMN-NLME model

EM.sim.NL<- function(y,#response vector (N x 1)
                     x,#covariates matrix (N x p1)
                     ind,#factor whose levels indicates the subjects or groups (N x 1) 
                     nlfder,#first derivative of the nonlinear function 
                     nlf,#nonlinear function with args: x1,...,xp1, beta1,...,betap, b1,...,bq
                     distr,#distribution to be used: "sn" (normal), "st" (student's-t), "ss" (slash), or "scn" (contaminated normal)
                     beta1,#p x 1 vector with initial values for beta
                     sigmae,#initial value for sigma2
                     D1,#q x q matrix with initial values for D 
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
  if (estimatenu){
    teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
  } else teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)])
  #
  if (is.null(bi)) bi <- matrix(0,ncol=q1,nrow=m)
  #
  criterio<-10
  count<-0
  llji = 1
  
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
    res_emj = revert_list(tapply(1:N,ind,emjs,y=ytil, x=Wtil, z=Htil, beta1=beta1, D1=D1,
                                 sigmae=sigmae, distr=distr,nu=nu))
    sum1 = Reduce("+",res_emj$sum1)
    sum2 = Reduce("+",res_emj$sum2)
    sum3 = sum(unlist(res_emj$sum3))
    sum4 = Reduce("+",res_emj$sum4)
    uj = unlist(res_emj$uj,use.names = F)
    bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
    
    # maximization step
    beta1<-solve(sum1)%*%sum2
    sigmae<-as.numeric(sum3)/N
    D1<-sum4/m
    #
    param <- teta
    #
    if (estimatenu){
      logvero1<-function(nu){logveros(ytil, Wtil, Htil, ind, beta1, sigmae, D1, distr, nu)}
      #
      if (distr=="sn"){ nu<-NULL} else
      {
        nu <- optim(nu,(logvero1),gr = NULL,method = "L-BFGS-B", lower =lb, upper = lu,control = list(fnscale=-1))$par
      }
      teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)],nu)
    } else teta <- c(beta1,sigmae,D1[upper.tri(D1, diag = T)])
    loglikVec[count] <- logveros(ytil, Wtil, Htil, ind, beta1, sigmae, D1, distr, nu)
    if (count>2){
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
  ### update bi and ui
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
  res_emj = revert_list(tapply(1:N,ind,emjs,y=ytil, x=Wtil, z=Htil, beta1=beta1, D1=D1,
                               sigmae=sigmae, distr=distr,nu=nu))
  uj = unlist(res_emj$uj,use.names = F)
  bi = t(bind_cols(res_emj$bi))#t(matrix(unlist(res_emj$bi),nrow=q1))
  ###
  # creating object to return
  dd<-matrix.sqrt(D1)[upper.tri(D1, diag = T)]
  theta = c(beta1,sigmae,dd,nu)
  if (distr=="sn") names(theta)<-c(paste0("beta",1:p),"sigma2",paste0("Dsqrt",1:length(dd)))
  else names(theta)<- c(paste0("beta",1:p),"sigma2",paste0("Dsqrt",1:length(dd)),paste0("nu",1:length(nu)))
  
  obj.out <- list(theta=theta, iter = count,estimates=list(beta=as.numeric(beta1),sigma2=sigmae,
                                                           D=D1),
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
      Infmatrixs(ytil,Wtil,Htil,ind,beta1,sigmae,D1,distr = distr,nu = nu)
      ,silent = T)
    if (class(desvios)=="try-error") {
      warning("Numerical error in calculating standard errors")
      obj.out$std.error=NULL
    } else{
      desvios <- c(desvios,rep(NA,length(nu)))
      names(desvios) <- names(theta)
      obj.out$std.error=desvios
    }
  }
  
  tf = Sys.time()
  obj.out$elapsedTime = as.numeric(difftime(tf,ti,units="secs"))
  obj.out$error=criterio
  obj.out
}


###############################################################
##### auxiliary functions
###############################################################

################################################################
#Log-likelihood - independent
################################################################
ljnormals <-function(j,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log(dmvnorm(y1,med,Psi))
}

ljts <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  dtj = gamma((nu+njj)/2)/gamma(nu/2)/pi^(njj/2)/sqrt(det(Psi))*nu^(-njj/2)*(dj/nu+1)^(-(njj+nu)/2)
  log(dtj)
}

ljss <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  f2 <- function(u) u^(nu - 1)*((2*pi)^(-njj/2))*(u^(njj/2))*((det(Psi))^(-1/2))*exp(-0.5*u*t(y1-med)%*%solve(Psi)%*%(y1-med))
  resp <- integrate(Vectorize(f2),0,1)$value
  log(nu*resp)
}

ljcns <-function(j,nu,y,x,z,beta1,D1,sigmae){
  y1=y[j]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[j,  ],ncol=p)
  z1=matrix(z[j ,  ],ncol=q1)
  med<-x1%*%beta1
  njj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(njj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  log((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
         (1-nu[1])*dmvnorm(y1,med,Psi)))
}
logveros = function(y,x,z,ind,beta1,sigmae,D1,distr,nu){ #ind = indicadora de individuo
  m<-n_distinct(ind)
  N<-length(ind)
  p<-dim(x)[2]
  q1<-dim(z)[2]
  
  if (distr=="sn") lv = sum(tapply(1:N,ind,ljnormals,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="st") lv = sum(tapply(1:N,ind,ljts,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="ss") lv = sum(tapply(1:N,ind,ljss,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  else if (distr=="scn") lv = sum(tapply(1:N,ind,ljcns,nu=nu,y=y,x=x,z=z,beta1=beta1,D1=D1,sigmae=sigmae))
  lv
}

##############################################################################
# EM - independent
##############################################################################
emjs = function(jseq, y, x, z, beta1,D1, sigmae,distr,nu) {
  y1=y[jseq]
  p= ncol(x);q1=ncol(z)
  x1=matrix(x[jseq,  ],ncol=p)
  z1=matrix(z[jseq,  ],ncol=q1)
  med<-x1%*%beta1
  nj = length(y1)
  Psi<-(z1)%*%(D1)%*%t(z1)+sigmae*diag(nj)
  dj<-as.numeric(t(y1-med)%*%solve(Psi)%*%(y1-med))
  #
  if  (distr=="sn"){
    uj<-1
  }
  if (distr=="st"){
    uj<-(nj+nu)/(dj+nu)
  }
  
  if (distr=="ss"){
    uj<-pgamma(1,nj/2+nu+1,dj/2)/pgamma(1,nj/2+nu,dj/2)*(nj+2*nu)/dj
  }
  
  if (distr=="scn"){
    fy<-as.numeric((nu[1]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))
    uj<-as.numeric((nu[1]*nu[2]*dmvnorm(y1,med,(Psi/nu[2]))+
                      (1-nu[1])*dmvnorm(y1,med,Psi)))/fy
  }
  
  bi<-D1%*%t(z1)%*%solve(Psi)%*%(y1-med)
  
  Tbj<-solve(solve(D1)+t(z1)%*%z1/sigmae)
  r<-Tbj%*%t(z1)%*%(y1-x1%*%beta1)/sigmae
  ub<-uj*r
  ub2j<-Tbj+uj*r%*%t(r)
  #
  sum1<-uj*t(x1)%*%x1 #denom beta
  sum2<-(t(x1)%*%(uj*y1-z1%*%ub)) #num beta
  sum3<-uj*t(y1-x1%*%beta1)%*%(y1-x1%*%beta1)-t(y1-x1%*%beta1)%*%z1%*%ub-
    t(ub)%*%t(z1)%*%(y1-x1%*%beta1)+traceM(ub2j%*%t(z1)%*%z1) #soma do sig2
  sum4<-ub2j #soma do Gamma
  obj.out = list(sum1=sum1,sum2=sum2,sum3=sum3,sum4=sum4,uj=uj)
  obj.out$bi=bi
  return(obj.out)
}

####
#inf mat (using first derivative and codes from the asymmetric case)
#Information matrix for SMN-LMM and SMN-LMM-AR(p) with E(bi)=0
scoreis <- function(jseq,y,x,z,beta1,sigmae,D1,distr,nu) {
  if (distr=="sn") c.=-sqrt(2/pi)
  if (distr=="st") c.=-sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
  if (distr=="ss") c.=-sqrt(2/pi)*nu/(nu-.5)
  if (distr=="scn") c.=-sqrt(2/pi)*(1+nu[1]*(nu[2]^(-.5)-1))
  lambda=rep(0,nrow(D1))
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


Infmatrixs <- function(y,x,z,ind,beta1,sigmae,D1,distr,nu){
  N <-length(y)
  score_list=tapply(1:N,ind,scoreis,y=y, x=x, z=z, beta1=beta1, sigmae=sigmae,D1=D1,distr=distr,nu=nu)
  mi_list = lapply(score_list,function(tt) {xm = matrix(tt,ncol=1);xm%*%t(xm)})
  infmat <- Reduce("+",mi_list)
  npar <- nrow(infmat);nd<-nrow(D1)
  infmat<- infmat[1:(npar-nd),1:(npar-nd)]
  if (abs(det(infmat))<1e-5) infmat= infmat+1e-20*diag(nrow(infmat))
  sqrt(diag(solve(infmat)))
}

