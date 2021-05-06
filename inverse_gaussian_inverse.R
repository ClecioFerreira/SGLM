mlg_estIGIV=function(Y,X,Nc,K,alpha,phi){
  X=as.matrix(X)
  q1=ncol(Nc)
  p=ncol(X)
  lteta=p+q1
  n=length(Y)
  ob=glm(Y~X,family = inverse.gaussian(link = "inverse"))
   beta0=as.matrix(as.vector(ob$coefficients[2:(p+1)]))
  f0=solve(t(Nc)%*%Nc+alpha*K)%*%t(Nc)%*%as.matrix(as.vector(ob$residuals)) 
  eta=as.vector(X%*%beta0+Nc%*%f0)
  mu=1/eta
  dev0=sum((as.vector(Y)-as.vector(mu))^2)
  thetachute=rbind(beta0,f0)
  criterio=1
  cont=0
  MI=matrix(0,p+q1,p+q1)
  U=matrix(0,p+q1,1)
  H=X%*%solve((t(X)%*%X))%*%t(X)
  H1=diag(n)-H
   mu=1/eta
    W=diag(1/eta)
   V=diag(mu^3)
    Wm=diag(sqrt(1/eta))
    Vm=diag(sqrt(mu^3))
    Vminv=diag(1/sqrt(mu^3))
    
  while((criterio>1e-4)&&(cont<1000)){
    cont=cont+1
    
     I_beta2=t(X)%*%W%*%X                                           
    I_beta_f=t(Nc)%*%W%*%X                                          
    I_f2=t(Nc)%*%W%*%Nc+alpha*K/phi                      
    MI[1:p,1:p]=I_beta2                       
    MI[1:p,(p+1):lteta]=t(I_beta_f)                               
    MI[(p+1):lteta,1:p]=I_beta_f   
    MI[(p+1):lteta,(p+1):lteta]=I_f2 
    U[1:p,1]=t(X)%*%Wm%*%Vminv%*%(Y-mu) 
    U[(p+1):(p+q1),1]=t(Nc)%*%Wm%*%Vminv%*%(Y-mu) - alpha*K%*%f0/phi
    thetafinal=thetachute+solve(MI)%*%U   
    beta1=as.matrix(thetafinal[1:p])        
    f0=as.matrix(thetafinal[(p+1):(p+q1)])  
    eta=as.vector(X%*%beta1+Nc%*%f0)
    mu=1/eta
    W=diag(1/eta)
   V=diag(mu^3)
    Wm=diag(sqrt(1/eta))
    Vm=diag(sqrt(mu^3))
    Vminv=diag(1/sqrt(mu^3))
    
   dev=sum(((as.vector(Y)-as.vector(mu))^2)/(as.vector(Y)*(as.vector(mu)^2))) 
                         
    phi=n/dev                                  
    criterio=sum(abs(dev-dev0))                                              
    thetachute=thetafinal                                                                    
    dev0=dev                
  }
   H=X%*%solve((t(X)%*%X))%*%t(X)
  trH=sum(diag(H))
   Z=solve(Vm%*%Wm)%*%(Y-mu)+X%*%beta1+Nc%*%f0 
  
  auxCV=as.numeric(t(Z-X%*%beta1)%*%W%*%(Z-X%*%beta1))
  GCV=as.numeric(auxCV/(1-trH/n)^2)
   Hn=Nc%*%solve(t(Nc)%*%Nc+alpha/phi*K)%*%t(Nc)
  trHn=sum(diag(Hn))
   GL=trH+trHn  
   
   obj = list(thetafinal = thetafinal, W = W, mu = mu, phi = phi, V=V, GCV=GCV, Wm=Wm, Vm=Vm, MI=(phi*MI), dev=dev, GL=GL)
}

CVspalphaIGIV<-function(alpha,Y,X,t,Nc,K,phi,method,ndx){
  p=ncol(X)
   n=length(Y)
  param=mlg_estIGIV(Y,X,Nc,K,alpha,phi)
  theta=param$thetafinal
  W=param$W
  V=param$V
  Wm=param$Wm
  Vm=param$Vm
  mu=param$mu
  beta1=theta[1:p]
  f=theta[(p+1):length(theta)]

   Z=solve(Vm%*%Wm)%*%(Y-mu)+X%*%beta1+Nc%*%f 
   
   MI=param$MI/param$phi
   H=cbind(X,Nc)%*%solve(MI)%*%rbind(t(X),t(Nc))%*%W
   trH=sum(diag(H))
   
  Z=solve(Vm%*%Wm)%*%(Y-mu)+X%*%beta1+Nc%*%f 
  
  auxCV=as.numeric(t(Z-X%*%beta1)%*%W%*%(Z-X%*%beta1))
  GCV=as.numeric(auxCV/(1-trH/n)^2)

  }