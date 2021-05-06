source('main_function.r')
require(splines)

###################### Simulation
# True values
n=100
sigma2=0.1
phi_real=100
b1=4
b2=-1
family="gamma"
link="log"
method="eilers"
curve="coseno" #coseno or doppler
ndx=15
limits_alpha=c(0,10)

if(curve=="coseno"){
t0=seq(-3*pi/2,5*pi/2, length=n/2)
t=sort(c(t0,t0))
ft=matrix(cos(t))
}

if(curve=="doppler"){
t0=t0=seq(0.66,1.66,length=n/2)
t=sort(c(t0,t0))
ft=cos(4*pi*t)*exp(-t^2/2)
}


nk=1
for(i in 2:n){ 
  aux=t[i]-t[i-1]
  nki=i*(aux>0)
  nk=c(nk,nki)
}
t0=t[nk>0]

x1=runif(n,0.1,1)
x2=runif(n,0,2)

mup=b1*x1+x2*b2
etaV= mup + ft 
X=cbind(x1,x2)

if(link=="identity"){
mu=etaV 
}
if(link=="log"){
  mu=exp(etaV)
}
if(link=="inverse"){
  mu=1/etaV 
}
if(link=="cauchit"){
  mu=pcauchy(etaV) 
}
if(link=="cloglog"){
  mu=1-exp(-exp(etaV))
}
if(link=="logit"){
  mu=exp(etaV)/(1+exp(etaV)) 
}
if(link=="probit"){
  mu=pnorm(etaV)
}
if(link=="sqrt"){
  mu=etaV^2 
}
if(link=="1mu2"){
  mu=1/sqrt(etaV)
}


if(family=="normal"){
    y=rnorm(length(mu),mu,sigma2)
  }
if(family=="binomial"){
    y=rbinom(length(mu),1,mu)
  }
if(family=="poisson"){
    y=rpois(length(mu),mu)
  }
if(family=="gamma"){
     y=rgamma(length(mu),phi_real, phi_real/mu)
  }
if(family=="inverse_gaussian"){
    library(statmod)
    y=rinvgauss(length(mu),mu,dispersion = 1/sigma2)
}
if(family=="nb"){
    y=rnbinom(length(mu),1,mu=mu)
  }
  
Y=matrix(y)
  
thetaest= GLMSPCV(Y,X,t,limits_alpha,family = family,link = link,method = method, ndx = ndx)
  
 
p=ncol(X)
betaest=thetaest$beta
epbeta=thetaest$epb
muest=as.matrix(thetaest$mu)
phi=thetaest$phi
alpha=thetaest$alpha
Var_theta=thetaest$MIf
Varf= Var_theta[(p+1):ncol(Var_theta),(p+1):ncol(Var_theta)]  
fest=as.matrix(thetaest$f)
Nc=as.matrix(thetaest$N)

# Parametric estimates
P<-cbind(betaest,epbeta); colnames(P)<-c("estimate","s.e.")
rownames(P)<-c(paste("beta",1:p,sep=""))
P
print(phi)

# Plot of the true and estimated non-parametric curves
  
if(curve=="cosseno"){
xlim=c(-4,7) 
ylim=c(-2,2) 
}else
if(curve=="doppler"){
xlim=c(0.66,1.66) 
ylim=c(-1,1) 
}

munp=as.vector(Nc%*%fest)
Varynp=Nc%*%solve(Varf)%*%t(Nc)
epf=sqrt(diag(Varynp))
bup=munp+2*epf
bdo=munp-2*epf
  
plot(t,munp,ylab="f(Time) Predicted",xlab="Time",font.lab=2,ps=0.1, cex.axis=1.5, cex.lab=1.5,type="n", ylim=c(min(bdo),max(bup)), xlim=c(min(t),max(t)))
lines(t,bdo,cex=2,lty=3,lwd=3)
lines(t,bup,cex=2,lty=3,lwd=3)
polygon(c(t,rev(t)), c(bup,rev(bdo)),col="gray")
lines(t,munp,cex=2,lty=1,lwd=3)


 
