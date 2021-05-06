bspline=function(x,xl=min(x),xr=max(x),ndx=ndx,bdeg=3){
  dx=(xr-xl)/ndx
  knots=seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
  B=splineDesign(knots,x,bdeg+1,0*x,outer.ok=T)
  
}

N <-function(t,t0){
  n=length(t)
  r=length(t0)
  N=matrix(0,n,r)
  for (i in 1:n){
    for (k in 1:r){
      if (t[i]==t0[k]) N[i,k]=1
    }
  }
  return(N)
}



R <-function(t0){
  r=length(t0)
  h=rep(0,r-1)
  for (i in 1:(r-1))  h[i]=t0[i+1]-t0[i]
  r0=matrix(0,r-1,r-1)
  for (i in 2:(r-2)){
    for (k in 1:(r-2)){
      if (abs(i-k)<2){
        r0[i,i]=1/3*(h[i-1]+h[i])
        r0[i,i+1]=h[i]/6
        r0[i+1,i]=h[i]/6
      }
    }
  }
  R=r0[2:(r-1),2:(r-1)]
  return(R)
}




Q <- function(t0){
  r=length(t0)
  h=rep(0,r-1)
  for (i in 1:(r-1))  h[i]=t0[i+1]-t0[i]
  q1=matrix(0,r,r-1)
  for (i in 1:r){
    for (k in 2:(r-1)){
      if (abs(i-k)<2){
        q1[k-1,k]=1/h[k-1]
        q1[k,k]=-(1/h[k-1]+1/h[k])
        q1[k+1,k]=1/h[k]
      }
    }
  }
  Q=q1[1:r,2:(r-1)]
  return(Q)
}

getK=function(method,t,ndx,n){
if(method=="green"){
    nk=1
    for(i in 2:n){      
      aux=t[i]-t[i-1]      
      nki=i*(aux>0)
      nk=c(nk,nki)
    }
    t0=t[nk>0]
    Nc=N(t,t0)
    Rc=R(t0)
    Qc=Q(t0)
    K=Qc%*%solve(Rc)%*%t(Qc)
  }
  else
  if(method=="eilers"){
    require(splines)
    Nc=bspline(t, ndx=ndx)
    Dff<-diag(ncol(Nc)) 
    for (k in 1:2) {
      Dff<-diff(Dff)
    }
    K=t(Dff)%*%Dff
  }
  else{
     stop("Method error")
  }
  obj=list(K=K,Nc=Nc)
  return(obj)
  }


GLMSPCV <-function(Y,X,t,alpha,phi=NULL,family = c("normal", "binomial", "gamma", "inverse_gaussian", "poisson", "nb"),link=c("identity", "log", "inverse", "logit", "probit", "cauchit", "cloglog", "1mu2", "sqrt"), method = c("green","eilers"), ndx=10, optim=TRUE, by=(max(alpha)-min(alpha))/100){ 
  X=as.matrix(X)  
  stopifnot(!missing(Y),!missing(t),!missing(X),nrow(X)==nrow(as.matrix((Y))))
  if(suppressWarnings(all(family==c("normal", "binomial", "gamma", "inverse_gaussian", "poisson", "nb"))))
  family="normal"
 
  if(suppressWarnings(all(link==c("identity", "log", "inverse", "logit", "probit", "cauchit", "cloglog", "1mu2", "sqrt"))))
  link="identity"
  
  if(length(alpha)==1)
  optim=FALSE
  
  require(MASS)
  n=nrow(as.matrix((Y)))
  
  if(is.null(phi)){
    H=X%*%solve((t(X)%*%X))%*%t(X)
    phi=as.numeric(1/n*t(Y)%*%(diag(n)-H)%*%Y)
  }
  
  thetaest_greed_gcv=c()
  tempo_decorrido=c()
  alphagreed=c()
  
  
  auxKNc=getK(method,t,ndx,n)
  K=auxKNc$K
  Nc=auxKNc$Nc
  if(family=="normal"){
     if(link=="identity"){
        source(paste(getSrcDirectory(GLMSPCV), "/normal_identity.R", sep=""))
        alphaest=alpha
        if(optim){       
         alphagreed = seq(min(alpha),max(alpha),by=by)  
         thetaest_greed_gcv=c()
         tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
   thetaest_greed_gcv[contador_greed] = CVspalphaNID(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
  )[3])
         GCV=min(thetaest_greed_gcv)
         alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
         }
  theta=mlg_estNID(Y,X,Nc,K,alphaest,phi)
 GCV=theta$GCV
     
     }
     else
     if(link=="log"){
        source(paste(getSrcDirectory(GLMSPCV), "/normal_log.R", sep=""))  
        alphaest=alpha
        if(optim){ 
       alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaNL(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}
theta=mlg_estNL(Y,X,Nc,K,alphaest,phi)
GCV=theta$GCV 
         }
     else
     if(link=="inverse"){
        source(paste(getSrcDirectory(GLMSPCV), "/normal_inverse.R", sep=""))   
        alphaest=alpha
        if(optim){ 
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaNIV(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}
theta=mlg_estNIV(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
  
     }    
     else{
        stop("Normal link erro")
     }
          
  }
  ##binomial
  else
  if (family=="binomial"){     
     if(link=="logit"){
        source(paste(getSrcDirectory(GLMSPCV), "/binomial_logit.R", sep="")) 
        alphaest=alpha
        if(optim){   
        alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaBL(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}
 
  theta=mlg_estBL(Y,X,Nc,K,alphaest,phi)
 GCV=theta$GCV
     }
     else
     if(link=="probit"){
        source(paste(getSrcDirectory(GLMSPCV), "/binomial_probit.R", sep=""))   
        alphaest=alpha
        if(optim){
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaBP(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}
 theta=mlg_estBP(Y,X,Nc,K,alphaest,phi)
 GCV=theta$GCV
     }
     else
     if(link=="cauchit"){
        source(paste(getSrcDirectory(GLMSPCV), "/binomial_cauchit.R", sep=""))  
        alphaest=alpha
        if(optim){
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaBC(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}
   theta=mlg_estBC(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }            
     else
     if(link=="cloglog"){
        source(paste(getSrcDirectory(GLMSPCV), "/binomial_cloglog.R", sep="")) 
        alphaest=alpha
        if(optim){ 
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaBCLL(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}

   theta=mlg_estBCLL(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }     
     else{
        stop("Binomial link erro")
     } 
  }
  ##gamma
  else
  if (family=="gamma"){
  if(link=="identity"){
        source(paste(getSrcDirectory(GLMSPCV), "/gamma_identity.R", sep=""))   
        alphaest=alpha
        if(optim){ 
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaGID(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}
  theta=mlg_estGID(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }
     else
     if(link=="log"){
        source(paste(getSrcDirectory(GLMSPCV), "/gamma_log.R", sep=""))    
        alphaest=alpha
        if(optim){  
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaGL(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
}
  theta=mlg_estGL(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }
     else                
     if(link=="inverse"){
        source(paste(getSrcDirectory(GLMSPCV), "/gamma_inverse.R"))     
        alphaest=alpha
        if(optim){
        alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaGIV(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]   
         }
 
   theta=mlg_estGIV(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }    
     else{
        stop("Gamma link erro")
     }
  }  
  #InverseGaussian
  else
  if(family=="inverse.gaussian"|family=="inverse_gaussian"){
     if(link=="identity"){
        source(paste(getSrcDirectory(GLMSPCV), "/inverse_gaussian_identity.R", sep=""))
        alphaest=alpha
        if(optim){
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaIGID(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                  
         }
 
   theta=mlg_estIGID(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }
     else
     if(link=="log"){
        source(paste(getSrcDirectory(GLMSPCV), "/inverse_gaussian_log.R", sep=""))
        alphaest=alpha
        if(optim){
        alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaIGL(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                  
         }
 
   theta=mlg_estIGL(Y,X,Nc,K,alphaest,phi)
 GCV=theta$GCV
     }
     else
     if(link=="inverse"){
        source(paste(getSrcDirectory(GLMSPCV), "/inverse_gaussian_inverse.R", sep=""))
        alphaest=alpha
        if(optim){
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaIGIV(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                  
         }
 
   theta=mlg_estIGIV(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }
     else
     if(link=="1/mu^2"|link=="1mu2"){
        source(paste(getSrcDirectory(GLMSPCV), "/inverse_gaussian_1mu2.R", sep=""))
        alphaest=alpha
        if(optim){
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaIG1MU2(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                  
         }
 
   theta=mlg_estIG1MU2(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }    
     else{
        stop("Gaussian Inverse link erro")
     }
  }   
  #Poisson   
  else
  if(family=="poisson"){
     if(link=="log"){
        source(paste(getSrcDirectory(GLMSPCV), "/poisson_log.R", sep=""))   
        alphaest=alpha
        if(optim){   
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaPL(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
         }
 
   theta=mlg_estPL(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }
     else
     if(link=="identity"){
        source(paste(getSrcDirectory(GLMSPCV), "/poisson_identity.R", sep="")) 
        alphaest=alpha
        if(optim){ 
        alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaPID(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                  
         }
 
   theta=mlg_estPID(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }
     else
     if(link=="sqrt"){
        source(paste(getSrcDirectory(GLMSPCV), "/poisson_sqrt.R", sep=""))      
        alphaest=alpha
        if(optim){
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaPSQRT(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                 
         }
 
   theta=mlg_estPSQRT(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
     }    
     else{
        stop("Poisson link erro")
     }
  }
  #Binomial Negativa
  else
  if(family=="nb"){
     if(link=="log"){
        source(paste(getSrcDirectory(GLMSPCV), "/nb_log.R", sep=""))   
        alphaest=alpha
        if(optim){                                       
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaNBL(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                 
         }
 
  theta=mlg_estNBL(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
  
     }
     else
     if(link=="identity"){
        source(paste(getSrcDirectory(GLMSPCV), "/nb_identity.R", sep=""))     
        alphaest=alpha
        if(optim){ 
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaNBID(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                  
         }
 
   theta=mlg_estNBID(Y,X,Nc,K,alphaest,phi)
  GCV=theta$GCV
    
     }
     else
     if(link=="sqrt"){
        source(paste(getSrcDirectory(GLMSPCV), "/nb_sqrt.R", sep=""))  
        alphaest=alpha
        if(optim){    
         alphagreed = seq(min(alpha),max(alpha),by=by)  
thetaest_greed_gcv=c()
tempo_decorrido=as.numeric(system.time( 
  for(contador_greed in 1:length(alphagreed)){   
    thetaest_greed_gcv[contador_greed] = CVspalphaNBSQRT(alphagreed[contador_greed],Y,X,t, Nc, K, phi, method, ndx)      
  }
)[3])
GCV=min(thetaest_greed_gcv)
alphaest=alphagreed[which.min(thetaest_greed_gcv)]                
         }
 
   theta=mlg_estNBSQRT(Y,X,Nc,K,alphaest,phi)
 GCV=theta$GCV
     }    
     else{
        stop("Binomial link erro")
     }
  }
  else{
     stop("Family erro")
     } 
        
    p=ncol(X)
        
    MIf=theta$MI                          
   
  ep=sqrt(diag(solve(MIf, tol=1e-23)))  #erro padrão
  betaest=theta$thetafinal[1:p]
  fest=theta$thetafinal[(p+1):length(theta$thetafinal)]
  epf=ep[(p+1):length(ep)]
  epbeta=ep[1:p]
 
  mu=theta$mu
 
  
  
  obj=list(thetafinal = theta$thetafinal, beta=betaest, f=fest, W = theta$W, mu = theta$mu, phi = theta$phi, N = Nc, alpha=alphaest, epb=epbeta, epf=epf, MIf=MIf, GCV=GCV,  dev=theta$dev, GL=theta$GL)
  if(optim){
    obj=list(thetafinal = theta$thetafinal, beta=betaest, f=fest, W = theta$W, mu = theta$mu, phi = theta$phi, N = Nc, alpha=alphaest, epb=epbeta, epf=epf, MIf=MIf, GCV=GCV,  dev=theta$dev, GL=theta$GL, greed_gcv=thetaest_greed_gcv, greed_alpha=alphagreed, time=tempo_decorrido)
  }
  return(obj)
  }