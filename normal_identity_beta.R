mlg_estNID_beta=function(Y,X,Nc,K,alpha,phi,beta,posbeta){
  X=as.matrix(X)
  q1=ncol(Nc)
  p=ncol(X)
  lteta=p+q1
  n=length(Y)
  beta=as.vector(beta)
  posbeta=as.vector(posbeta)
  
  auxbeta0=matrix(0,p,1)
  beta0=rep(NA,ncol(X))
   
  beta0[posbeta]=beta
  j=1
  for(i in 1:ncol(X)){
  if(is.na(beta0[i])){
  beta0[i]=auxbeta0[j]
  j=j+1
  }
  }
  beta0=as.matrix(beta0)
  f0=solve(t(Nc)%*%Nc+alpha*K)%*%t(Nc)%*%(Y-X%*%beta0) 
  eta=as.vector(X%*%beta0+Nc%*%f0)
  mu=eta
  dev0=sum((as.vector(Y)-as.vector(mu))^2)
  thetachute=rbind(beta0,f0)
  criterio=1
  cont=0
  MI=matrix(0,p+q1,p+q1)
  U=matrix(0,p+q1,1)
  H=X%*%solve((t(X)%*%X))%*%t(X)
  H1=diag(n)-H
  ##W, Wm, V, Vm dentro do while nas outras
  W=diag(n)
  Wm=diag(n)
  V=diag(n)
  Vm=V
  Vminv=V 
    mu=eta
  while((criterio>1e-4)&&(cont<1000)){
    cont=cont+1
  
  
  
    
    
    #Z=solve(Vm%*%Wm)%*%(Y-mu)+X%*%beta0+Nc%*%f0                           #2
    #beta1=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%(Z-Nc%*%f0)                    #2
    #f1=solve(t(Nc)%*%W%*%Nc+alpha*K/phi)%*%t(Nc)%*%W%*%(Z-X%*%beta0)      #2
    
    
    
     I_beta2=t(X)%*%W%*%X                                                     #1
    I_beta_f=t(Nc)%*%W%*%X                                                   #1
    I_f2=t(Nc)%*%W%*%Nc+alpha*K/phi                                          #1
    MI[1:p,1:p]=I_beta2                                                      #1
    MI[1:p,(p+1):lteta]=t(I_beta_f)                                          #1
    MI[(p+1):lteta,1:p]=I_beta_f                                             #1
    MI[(p+1):lteta,(p+1):lteta]=I_f2                                         #1
    U[1:p,1]=t(X)%*%Wm%*%Vminv%*%(Y-mu)                                         #1
    U[(p+1):(p+q1),1]=t(Nc)%*%Wm%*%Vminv%*%(Y-mu) - alpha*K%*%f0/phi  
    
    
    
    auxthetafinal=thetachute[-posbeta]+solve(MI[-posbeta,-posbeta])%*%U[-posbeta,]
    
     thetafinal=rep(NA,length(thetachute))
  thetafinal[posbeta]=beta
  j=1
  for(i in 1:length(thetachute)){
  if(is.na(thetafinal[i])){
  thetafinal[i]=auxthetafinal[j]
  j=j+1
  }
  }
    
    
    beta1=as.matrix(thetafinal[1:p])                                         #1
    f0=as.matrix(thetafinal[(p+1):(p+q1)])                                   #1 
    eta=as.vector(X%*%beta1+Nc%*%f0)
     mu=eta
     
    dev=sum((as.vector(Y)-as.vector(mu))^2)                                  
    #print(beta1)                                  
    criterio=sum(abs(dev-dev0))                                              
                                                                             
    
     
   
    phi=n/dev                                                                
    
    thetachute=thetafinal                                                    #1
    
    #beta0=beta1                                                             #2
    #f0=f1                                                                   #2
    
                   

    dev0=dev                                                                
  }
   #obj modificado 10/9
   H=X%*%solve((t(X)%*%X))%*%t(X)
  trH=sum(diag(H))
   Z=solve(Vm%*%Wm)%*%(Y-mu)+X%*%beta1+Nc%*%f0 
  
  auxCV=as.numeric(t(Z-X%*%beta1)%*%W%*%(Z-X%*%beta1))
  GCV=as.numeric(auxCV/(1-trH/n)^2)
   Hn=Nc%*%solve(t(Nc)%*%Nc+alpha/phi*K)%*%t(Nc)
  trHn=sum(diag(Hn))
   GL=trH+trHn  
   
   obj = list(thetafinal = thetafinal, W = W, mu = mu, phi = phi, V=V, GCV=GCV, Wm=Wm, Vm=Vm, MI=(phi*MI), dev=dev, GL=GL)
    return(obj)
}

CVspalphaNID_beta<-function(alpha,Y,X,t,Nc,K,phi,method,ndx,kb){
  p=ncol(X)
   n=length(Y)
  param=mlg_estNID_beta(Y,X,Nc,K,alpha,phi,beta,posbeta)#[beta,f]
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
  # Hn=Nc%*%solve(t(Nc)%*%Nc+alpha/phi*K)%*%t(Nc)
  #trHn=sum(diag(Hn))
   #GL=trH+trHn
   
  #auxCV=as.numeric(t(Z-X%*%beta1)%*%W%*%(Z-X%*%beta1))
  #CV=auxCV/n*(t(1-H)%*%(1-H))                              
  
  }