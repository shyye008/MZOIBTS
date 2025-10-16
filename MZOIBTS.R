library(copula)
#####Beta(mu,phi)#####
{
  dbeta0<-function(x,mu,phi){
    a<-mu*phi
    b<-(1-mu)*phi
    dbeta(x,a,b)
  }
  
  pbeta0<-function(q,mu,phi){
    a<-mu*phi
    b<-(1-mu)*phi
    pbeta(q,a,b)
  }
  
  qbeta0<-function(p,mu,phi){
    a<-mu*phi
    b<-(1-mu)*phi
    qbeta(p,a,b)
  }
  
  rbeta0<-function(n,mu,phi){
    a<-mu*phi
    b<-(1-mu)*phi
    rbeta(n,a,b)
  }
}

#####ZOIB(p1,p2,mu,phi)#####
{
  dZOIB<-function(x,p1,p2,mu,phi){
    ifelse(x==0,1-p1,ifelse(x==1,p1*p2,p1*(1-p2)*dbeta0(x,mu,phi)))
  }
  
  pZOIB<-function(q,p1,p2,mu,phi){
    ifelse(q==0,1-p1,ifelse(q==1,1,1-p1+p1*(1-p2)*pbeta0(q,mu,phi)))
  }
  
  qZOIB<-function(p,p1,p2,mu,phi){
    ifelse(p<=1-p1,0,ifelse(p>1-p1*p2,1,qbeta0((p-1+p1)/(p1*(1-p2)),mu,phi)))
  }
  
  rZOIB<-function(n,p1,p2,mu,phi){
    d1<-rbinom(n,1,p1);d2<-rbinom(n,1,p2);x0<-rbeta0(n,mu,phi)
    ifelse(d1==0,0,ifelse(d2==1,1,x0))
  }
}

#####Generate ZOIBTS#####
rAR1<-function(n,copula,alpha){
  if(copula=="Gaussian"){
    cc<-normalCopula(alpha)
  } else if(copula=="Clayton"){
    cc<-claytonCopula(alpha)
  } else if(copula=="Gumbel"){
    cc<-gumbelCopula(alpha)
  } else if(copula=="Frank"){
    cc<-frankCopula(alpha)
  } else if(copula=="AMH"){
    cc<-amhCopula(alpha)
  }
  u<-runif(n+1)
  v<-u
  cCopula(cbind(u[-(n+1)],v[-1]),cc,inverse=T)[,2]
}

rZOIBTS<-function(n,copula,alpha,p1,p2,mu,phi){
  u<-rAR1(n,copula,alpha)
  qZOIB(u,p1,p2,mu,phi)
}

#####Generate MZOIBTS#####
rMZOIBTS<-function(n,copula,alpha,X1,beta1,X2,beta2,X3,beta3,X4,beta4){
  if(beta1==Inf){
    p1<-1
  }else{
    a1<-as.vector(X1%*%beta1)
    p1<-exp(a1)/(1+exp(a1))
  }
  
  if(beta2==-Inf){
    p2<-0
  }else{
    a2<-as.vector(X2%*%beta2)
    p2<-exp(a2)/(1+exp(a2))
  }
  
  a3<-as.vector(X3%*%beta3)
  
  v<-exp(a3)/(1+exp(a3))
  mu<-(v/p1-p2)/(1-p2)
  
  phi<-exp(as.vector(X4%*%beta4))
  
  rZOIBTS(n,copula,alpha,p1,p2,mu,phi)
}

library(Matrix)
#####Fit ZOIB#####
est.ILL<-function(X1,X2,X3,X4,y){
  ILL_all<-function(theta,X1,X2,X3,X4,y){
    dim4<-ncol(X4)
    if(sum(y==0)<1 & sum(y==1)<1){
      dim3<-ncol(X3)
      p1<-1;p2<-0
      mu<-exp(as.vector(X3%*%theta[1:dim3]))/(1+exp(as.vector(X3%*%theta[1:dim3])))
      phi<-exp(as.vector(X4%*%theta[(dim3+1):(dim3+dim4)]))
      l<-sum(log(dZOIB(y,p1,p2,mu,phi)))
    }else if(sum(y==0)<1 & sum(y==1)>=1){
      dim2<-ncol(X2);dim3<-ncol(X3)
      p1<-1
      p2<-exp(as.vector(X2%*%theta[1:dim2]))/(1+exp(as.vector(X2%*%theta[1:dim2])))
      v<-exp(as.vector(X3%*%theta[(dim2+1):(dim2+dim3)]))/(1+exp(as.vector(X3%*%theta[(dim2+1):(dim2+dim3)])))
      mu<-(v-p2)/(1-p2)
      phi<-exp(as.vector(X4%*%theta[(dim2+dim3+1):(dim2+dim3+dim4)]))
      l<-sum(log(dZOIB(y,p1,p2,mu,phi)))
    }else if(sum(y==0)>=1 & sum(y==1)<1){
      dim1<-ncol(X1);dim3<-ncol(X3)
      p1<-exp(as.vector(X1%*%theta[1:dim1]))/(1+exp(as.vector(X1%*%theta[1:dim1])))
      p2<-0
      v<-exp(as.vector(X3%*%theta[(dim1+1):(dim1+dim3)]))/(1+exp(as.vector(X3%*%theta[(dim1+1):(dim1+dim3)])))
      mu<-v/p1
      phi<-exp(as.vector(X4%*%theta[(dim1+dim3+1):(dim1+dim3+dim4)]))
      l<-sum(log(dZOIB(y,p1,p2,mu,phi)))
    }else{
      dim1<-ncol(X1);dim2<-ncol(X2);dim3<-ncol(X3)
      p1<-exp(as.vector(X1%*%theta[1:dim1]))/(1+exp(as.vector(X1%*%theta[1:dim1])))
      p2<-exp(as.vector(X2%*%theta[(dim1+1):(dim1+dim2)]))/(1+exp(as.vector(X2%*%theta[(dim1+1):(dim1+dim2)])))
      v<-exp(as.vector(X3%*%theta[(dim1+dim2+1):(dim1+dim2+dim3)]))/
        (1+exp(as.vector(X3%*%theta[(dim1+dim2+1):(dim1+dim2+dim3)])))
      mu<-(v/p1-p2)/(1-p2)
      phi<-exp(as.vector(X4%*%theta[(dim1+dim2+dim3+1):(dim1+dim2+dim3+dim4)]))
      l<-sum(log(dZOIB(y,p1,p2,mu,phi)))
    }
    -l
  }
  
  d1<-ifelse(y>0,1,0)
  y2<-y[d1==1];d2<-ifelse(y2==1,1,0);X2_2<-X2[d1==1,]
  if(sum(y==0)<1 & sum(y==1)<1){
    if(dim(X3)[2]==1){
      b3<-coef(lm(y~1))
    } else {
      b3<-coef(lm(y~X3[,-1]))
    }
    
    mu<-exp(as.vector(X3%*%b3))/(1+exp(as.vector(X3%*%b3)))
    y0<-y[y!=0 & y!=1];mu0<-mu[y!=0 & y!=1]
    phi<-mean(mu0*(1-mu0))/var(y0-mu0)-1
    dim4<-ncol(X4)
    initial<-c(b3,log(phi*2),rep(0,dim4-1))
  }else if(sum(y==0)<1 & sum(y==1)>=1){
    if(dim(X2)[2]==1){
      fit2<-glm(d2~1,family=binomial)
      b2<-coef(fit2)
      p2<-mean(predict(fit2,type="response"))
      
    } else {
      fit2<-glm(d2~X2_2[,-1],family=binomial)
      b2<-coef(fit2)
      p2<-predict(fit2,newdata=X2,type="response")
    }
    if(dim(X3)[2]==1){
      b3<-coef(lm(y~1))
    } else {
      b3<-coef(lm(y~X3[,-1]))
    }
    
    v<-exp(as.vector(X3%*%b3))/(1+exp(as.vector(X3%*%b3)))
    mu<-(v-p2)/(1-p2)
    y0<-y[y!=0 & y!=1];mu0<-mu[y!=0 & y!=1]
    phi<-mean(mu0*(1-mu0))/var(y0-mu0)-1
    dim4<-ncol(X4)
    initial<-c(b2,b3,log(phi*2),rep(0,dim4-1))
  }else if(sum(y==0)>=1 & sum(y==1)<1){
    if(dim(X1)[2]==1){
      fit1<-glm(d1~1,family=binomial)
      b1<-coef(fit1)
      p1<-mean(predict(fit1,type="response"))
    } else {
      fit1<-glm(d1~X1[,-1],family=binomial)
      b1<-coef(fit1)
      p1<-predict(fit1,type="response")
    }
    if(dim(X3)[2]==1){
      b3<-coef(lm(y~1))
    } else {
      b3<-coef(lm(y~X3[,-1]))
    }
    
    v<-exp(as.vector(X3%*%b3))/(1+exp(as.vector(X3%*%b3)))
    mu<-v/p1
    y0<-y[y!=0 & y!=1];mu0<-mu[y!=0 & y!=1]
    phi<-mean(mu0*(1-mu0))/var(y0-mu0)-1
    dim4<-ncol(X4)
    initial<-c(b1,b3,log(phi*2),rep(0,dim4-1))
  }else{
    if(dim(X1)[2]==1){
      fit1<-glm(d1~1,family=binomial)
      b1<-coef(fit1)
      p1<-mean(predict(fit1,type="response"))
    } else {
      fit1<-glm(d1~X1[,-1],family=binomial)
      b1<-coef(fit1)
      p1<-predict(fit1,type="response")
    }
    if(dim(X2)[2]==1){
      fit2<-glm(d2~1,family=binomial)
      b2<-coef(fit2)
      p2<-mean(predict(fit2,type="response"))
    } else {
      fit2<-glm(d2~X2_2[,-1],family=binomial)
      b2<-coef(fit2)
      p2<-predict(fit2,newdata=X2,type="response")
    }
    if(dim(X3)[2]==1){
      b3<-coef(lm(y~1))
    } else {
      b3<-coef(lm(y~X3[,-1]))
    }
    
    v<-exp(as.vector(X3%*%b3))/(1+exp(as.vector(X3%*%b3)))
    mu<-(v/p1-p2)/(1-p2)
    y0<-y[y!=0 & y!=1];mu0<-mu[y!=0 & y!=1]
    phi<-mean(mu0*(1-mu0))/var(y0-mu0)-1
    dim4<-ncol(X4)
    initial<-c(b1,b2,b3,log(phi*2),rep(0,dim4-1))
  }
  
  fit<-optim(initial,ILL_all,X1=X1,X2=X2,X3=X3,X4=X4,y=y,
             method="BFGS",hessian=T)
  list(fit$par,fit$value,fit$hessian)
}

#####Fit copula parameter#####
qest.cop<-function(theta,X1,X2,X3,X4,y,copula){
  qll.cop<-function(alpha,theta,X1,X2,X3,X4,y,copula){
    cop.like<-function(u2,u1,p21,p22,p11,p12,alpha,copula){
      if(copula=="Gaussian"){
        cc<-normalCopula(alpha)
      } else if(copula=="Clayton"){
        cc<-claytonCopula(alpha)
      } else if(copula=="Gumbel"){
        cc<-gumbelCopula(alpha)
      } else if(copula=="Frank"){
        cc<-frankCopula(alpha)
      } else if(copula=="AMH"){
        cc<-amhCopula(alpha)
      }
      
      if(u1==0){
        if(u2==0){
          cd<-pCopula(c(1-p11,1-p21),cc)
        } else if(u2==1){
          cd<-(1-p11-pCopula(c(1-p21*p22,1-p11),cc))
        } else{
          cd<-(cCopula(c(u2,1-p11),cc)[,2])
        }
      } else if(u1==1){
        if(u2==0){
          cd<-(1-p21-pCopula(c(1-p11*p12,1-p21),cc))
        } else if(u2==1){
          cd<-(p11*p12+p21*p22-1+pCopula(c(1-p21*p22,1-p11*p12),cc))
        } else{
          cd<-(1-cCopula(c(u2,1-p11*p12),cc)[,2])
        }
      } else{
        if(u2==0){
          cd<-cCopula(c(u1,1-p21),cc)[,2]
        } else if(u2==1){
          cd<-1-cCopula(c(u1,1-p21*p22),cc)[,2]
        } else{
          cd<-dCopula(c(u1,u2),cc)
        }
      }
      cd
    }
    
    n<-length(y)
    dim1<-ncol(X1);dim2<-ncol(X2);dim3<-ncol(X3);dim4<-ncol(X4)
    if(sum(y==0)<1 & sum(y==1)<1){
      p1<-rep(1,n);p2<-rep(0,n)
      v<-exp(as.vector(X3%*%theta[1:dim3]))/
        (1+exp(as.vector(X3%*%theta[1:dim3])))
      phi<-exp(as.vector(X4%*%theta[(dim3+1):(dim3+dim4)]))
    }else if(sum(y==0)<1 & sum(y==1)>=1){
      p1<-rep(1,n)
      p2<-exp(as.vector(X2%*%theta[1:dim2]))/(1+exp(as.vector(X2%*%theta[1:dim2])))
      v<-exp(as.vector(X3%*%theta[(dim2+1):(dim2+dim3)]))/
        (1+exp(as.vector(X3%*%theta[(dim2+1):(dim2+dim3)])))
      phi<-exp(as.vector(X4%*%theta[(dim2+dim3+1):(dim2+dim3+dim4)]))
    }else if(sum(y==0)>=1 & sum(y==1)<1){
      p1<-exp(as.vector(X1%*%theta[1:dim1]))/(1+exp(as.vector(X1%*%theta[1:dim1])))
      p2<-rep(0,n)
      v<-exp(as.vector(X3%*%theta[(dim1+1):(dim1+dim3)]))/
        (1+exp(as.vector(X3%*%theta[(dim1+1):(dim1+dim3)])))
      phi<-exp(as.vector(X4%*%theta[(dim1+dim3+1):(dim1+dim3+dim4)]))
    }else{
      p1<-exp(as.vector(X1%*%theta[1:dim1]))/(1+exp(as.vector(X1%*%theta[1:dim1])))
      p2<-exp(as.vector(X2%*%theta[(dim1+1):(dim1+dim2)]))/(1+exp(as.vector(X2%*%theta[(dim1+1):(dim1+dim2)])))
      v<-exp(as.vector(X3%*%theta[(dim1+dim2+1):(dim1+dim2+dim3)]))/
        (1+exp(as.vector(X3%*%theta[(dim1+dim2+1):(dim1+dim2+dim3)])))
      phi<-exp(as.vector(X4%*%theta[(dim1+dim2+dim3+1):(dim1+dim2+dim3+dim4)]))
    }
    mu<-(v/p1-p2)/(1-p2)
    
    u<-pZOIB(y,p1,p2,mu,phi)
    ll<-0
    for (i in 1:(n-1)) {
      ll<-ll+log(cop.like(u[i+1],u[i],p1[i+1],p2[i+1],p1[i],p2[i],alpha,copula))
    }
    -ll
  }
  
  if(copula=="Gaussian" | copula=="AMH"){
    a0<-0.4
  } else {
    a0<-2
  }
  initial<-a0
  
  if(copula=="Gaussian"){
    fit<-optim(initial,qll.cop,theta=theta,X1=X1,X2=X2,X3=X3,X4=X4,y=y,copula=copula,
               lower=-0.999,upper=0.999,method="Brent")
  } else if(copula=="Clayton"){
    fit<-optim(initial,qll.cop,theta=theta,X1=X1,X2=X2,X3=X3,X4=X4,y=y,copula=copula,
               lower=-0.999,upper=100,method="Brent")
  } else if(copula=="Gumbel"){
    fit<-optim(initial,qll.cop,theta=theta,X1=X1,X2=X2,X3=X3,X4=X4,y=y,copula=copula,
               lower=1.0001,upper=100,method="Brent")
  } else if(copula=="Frank"){
    fit<-optim(initial,qll.cop,theta=theta,X1=X1,X2=X2,X3=X3,X4=X4,y=y,copula=copula,
               lower=-100,upper=100,method="Brent")
  } else if(copula=="AMH"){
    fit<-optim(initial,qll.cop,theta=theta,X1=X1,X2=X2,X3=X3,X4=X4,y=y,copula=copula,
               lower=-0.999,upper=0.999,method="Brent")
  }
  list(fit$par,fit$value)
}

#####Fit MZOIBTS via parametric bootstrap#####
twostep.MZOIBTS<-function(R=100,X1,X2,X3,X4,y,copula,sig.level){
  dim1<-ncol(X1);dim2<-ncol(X2);dim3<-ncol(X3);dim4<-ncol(X4);n<-length(y)
  if(sum(y==0)<1 & sum(y==1)<1){
    name<-c("initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
  }else if(sum(y==0)<1 & sum(y==1)>=1){
    name<-c(rep("b2",dim2),"initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
  }else if(sum(y==0)>=1 & sum(y==1)<1){
    name<-c(rep("b1",dim1),"initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
  }else{
    name<-c(rep("b1",dim1),rep("b2",dim2),"initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
  }
  
  that<-est.ILL(X1,X2,X3,X4,y)[[1]]
  alpha<-qest.cop(that,X1,X2,X3,X4,y,copula)[[1]]
  bootEst<-c()
  for (i in 1:R) {
    if(sum(y==0)<1 & sum(y==1)<1){
      y_booti<-rMZOIBTS(n,copula,alpha,X1,Inf,X2,-Inf,X3,
                        that[1:dim3],X4,that[(dim3+1):(dim3+dim4)])
    }else if(sum(y==0)<1 & sum(y==1)>=1){
      y_booti<-rMZOIBTS(n,copula,alpha,X1,Inf,X2,that[1:dim2],X3,
                        that[(dim2+1):(dim2+dim3)],
                        X4,that[(dim2+dim3+1):(dim2+dim3+dim4)])
    }else if(sum(y==0)>=1 & sum(y==1)<1){
      y_booti<-rMZOIBTS(n,copula,alpha,X1,that[1:dim1],X2,-Inf,X3,
                        that[(dim1+1):(dim1+dim3)],
                        X4,that[(dim1+dim3+1):(dim1+dim3+dim4)])
    }else{
      y_booti<-rMZOIBTS(n,copula,alpha,X1,that[1:dim1],X2,that[(dim1+1):(dim1+dim2)],X3,
                        that[(dim1+dim2+1):(dim1+dim2+dim3)],
                        X4,that[(dim1+dim2+dim3+1):(dim1+dim2+dim3+dim4)])
    }
    
    that_booti<-try(est.ILL(X1,X2,X3,X4,y_booti),silent=T)
    while(class(that_booti)=="try-error" | length(that_booti[[1]])!=length(that)){
      if(sum(y==0)<1 & sum(y==1)<1){
        y_booti<-rMZOIBTS(n,copula,alpha,X1,Inf,X2,-Inf,X3,
                          that[1:dim3],X4,that[(dim3+1):(dim3+dim4)])
      }else if(sum(y==0)<1 & sum(y==1)>=1){
        y_booti<-rMZOIBTS(n,copula,alpha,X1,Inf,X2,that[1:dim2],X3,
                          that[(dim2+1):(dim2+dim3)],
                          X4,that[(dim2+dim3+1):(dim2+dim3+dim4)])
      }else if(sum(y==0)>=1 & sum(y==1)<1){
        y_booti<-rMZOIBTS(n,copula,alpha,X1,that[1:dim1],X2,-Inf,X3,
                          that[(dim1+1):(dim1+dim3)],
                          X4,that[(dim1+dim3+1):(dim1+dim3+dim4)])
      }else{
        y_booti<-rMZOIBTS(n,copula,alpha,X1,that[1:dim1],X2,that[(dim1+1):(dim1+dim2)],X3,
                          that[(dim1+dim2+1):(dim1+dim2+dim3)],
                          X4,that[(dim1+dim2+dim3+1):(dim1+dim2+dim3+dim4)])
      }
      
      that_booti<-try(est.ILL(X1,X2,X3,X4,y_booti),silent=T)
    }
    bootEst<-rbind(bootEst,that_booti[[1]])
  }
  p<-length(that)
  res<-matrix(0,p,5,dimnames=list(name, c('est','CI_lower','CI_upper','SE','p.value')))
  for (i in 1:p){
    res[i,1]<-signif(that[i],digits=4)
    res[i,2]<-signif(that[i]-qnorm(1-sig.level/2)*sd(bootEst[,i]),digits=4)
    res[i,3]<-signif(that[i]+qnorm(1-sig.level/2)*sd(bootEst[,i]),digits=4)
    res[i,4]<-signif(sd(bootEst[,i]),digits=4)
    ts<-(that[i]^2)/(res[i,4]^2)
    res[i,5]<-signif(1-pchisq(ts,1),digits=4)
  }
  list(res,bootEst)
}

#####Calculate score functions#####
fit.score<-function(X1,X2,X3,X4,y,b1=NULL,b2=NULL,b3,b4){
  n<-length(y)
  s3<-as.vector(X3%*%b3)
  s4<-as.vector(X4%*%b4)
  if(sum(y==0)<1 & sum(y==1)<1){
    mu<-1/(1+exp(-s3))
    dmu3<--exp(-s3)/(1+exp(-s3))^2
    dh3<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu3*X3
    dh4<-(digamma(exp(s4))-digamma(mu*exp(s4))*mu-digamma((1-mu)*exp(s4))*(1-mu))*exp(s4)*X4
    score3<-dh3+dmu3*exp(s4)*log(y/(1-y))*X3
    score4<-dh4+mu*exp(s4)*log(y/(1-y))*X4+exp(s4)*log(1-y)*X4
    score<-cbind(score3,score4)
  }else if(sum(y==0)<1 & sum(y==1)>=1){
    ytilde<-ifelse(y==1,y-0.0001,y)
    if(length(b2)==1){
      s2<-as.vector(b2*X2)
    } else{
      s2<-as.vector(X2%*%b2)
    }
    mu<-(1+exp(s2))/(1+exp(-s3))-exp(s2)
    dmu2<-exp(s2)/(1+exp(-s3))-exp(s2)
    dmu3<--(1+exp(s2))*exp(-s3)/(1+exp(-s3))^2
    dh2<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu2*X2
    dh3<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu3*X3
    dh4<-(digamma(exp(s4))-digamma(mu*exp(s4))*mu-digamma((1-mu)*exp(s4))*(1-mu))*exp(s4)*X4
    score2<--exp(s2)/(1+exp(s2))*X2+ifelse(y==1,1,0)*X2+
      ifelse(y<1,1,0)*(dh2+dmu2*exp(s4)*log(ytilde/(1-ytilde))*X2)
    score3<-ifelse(y<1,1,0)*(dh3+dmu3*exp(s4)*log(ytilde/(1-ytilde))*X3)
    score4<-ifelse(y<1,1,0)*(dh4+mu*exp(s4)*log(ytilde/(1-ytilde))*X4+exp(s4)*log(1-ytilde)*X4)
    score<-cbind(score2,score3,score4)
  }else if(sum(y==0)>=1 & sum(y==1)<1){
    ytilde<-ifelse(y==0,y+0.0001,y)
    if(length(b1)==1){
      s1<-as.vector(b1*X1)
    } else{
      s1<-as.vector(X1%*%b1)
    }
    mu<-(1+exp(-s1))/(1+exp(-s3))
    dmu1<--exp(-s1)/(1+exp(-s3))
    dmu3<--(1+exp(-s1))*exp(-s3)/(1+exp(-s3))^2
    dh1<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu1*X1
    dh3<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu3*X3
    dh4<-(digamma(exp(s4))-digamma(mu*exp(s4))*mu-digamma((1-mu)*exp(s4))*(1-mu))*exp(s4)*X4
    score1<--exp(s1)/(1+exp(s1))*X1+ifelse(y>0,1,0)*X1+
      ifelse(y>0,1,0)*(dh1+dmu1*exp(s4)*log(ytilde/(1-ytilde))*X1)
    score3<-ifelse(y>0,1,0)*(dh3+dmu3*exp(s4)*log(ytilde/(1-ytilde))*X3)
    score4<-ifelse(y>0,1,0)*(dh4+mu*exp(s4)*log(ytilde/(1-ytilde))*X4+exp(s4)*log(1-ytilde)*X4)
    score<-cbind(score1,score3,score4)
  }else{
    ytilde<-ifelse(y==0,y+0.0001,ifelse(y==1,y-0.0001,y))
    if(length(b1)==1){
      s1<-as.vector(b1*X1)
    } else{
      s1<-as.vector(X1%*%b1)
    }
    if(length(b2)==1){
      s2<-as.vector(b2*X2)
    } else{
      s2<-as.vector(X2%*%b2)
    }
    mu<-(1+exp(s2))*(1+exp(-s1))/(1+exp(-s3))-exp(s2)
    dmu1<--(1+exp(s2))*exp(-s1)/(1+exp(-s3))
    dmu2<-exp(s2)*(1+exp(-s1))/(1+exp(-s3))-exp(s2)
    dmu3<--(1+exp(s2))*(1+exp(-s1))*exp(-s3)/(1+exp(-s3))^2
    dh1<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu1*X1
    dh2<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu2*X2
    dh3<--(digamma(mu*exp(s4))-digamma((1-mu)*exp(s4)))*exp(s4)*dmu3*X3
    dh4<-(digamma(exp(s4))-digamma(mu*exp(s4))*mu-digamma((1-mu)*exp(s4))*(1-mu))*exp(s4)*X4
    score1<--exp(s1)/(1+exp(s1))*X1+ifelse(y>0,1,0)*X1+
      ifelse(y>0 & y<1,1,0)*(dh1+dmu1*exp(s4)*log(ytilde/(1-ytilde))*X1)
    score2<--ifelse(y>0,1,0)*exp(s2)/(1+exp(s2))*X2+ifelse(y==1,1,0)*X2+
      ifelse(y>0 & y<1,1,0)*(dh2+dmu2*exp(s4)*log(ytilde/(1-ytilde))*X2)
    score3<-ifelse(y>0 & y<1,1,0)*(dh3+dmu3*exp(s4)*log(ytilde/(1-ytilde))*X3)
    score4<-ifelse(y>0 & y<1,1,0)*(dh4+mu*exp(s4)*log(ytilde/(1-ytilde))*X4+exp(s4)*log(1-ytilde)*X4)
    score<-cbind(score1,score2,score3,score4)
  }
  score
}

#####HAC estimator#####
HAC<-function(score,L=10){
  n<-dim(score)[1]
  Jhat<-0
  for (i in 1:n) {
    for (j in 1:n) {
      l<-abs(i-j)
      wl<-max(0,1-l/(L+1))
      Jhat<-Jhat+wl*score[i,]%*%t(score[j,])
    }
  }
  Jhat/n
}

#####Fit MZOIBTS via robust sandwich#####
HAC.MZOIBTS<-function(X1,X2,X3,X4,y,sig.level,L=10){
  fit<-est.ILL(X1,X2,X3,X4,y)
  bhat<-fit[[1]]
  dim1<-ncol(X1);dim2<-ncol(X2);dim3<-ncol(X3);dim4<-ncol(X4);n<-length(y)
  if(sum(y==0)<1 & sum(y==1)<1){
    name<-c("initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
    b3<-bhat[1:dim3];b4<-bhat[(dim3+1):(dim3+dim4)]
    score<-fit.score(X1,X2,X3,X4,y,b3=b3,b4=b4)
  }else if(sum(y==0)<1 & sum(y==1)>=1){
    name<-c(rep("b2",dim2),"initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
    b2<-bhat[1:dim2];b3<-bhat[(dim2+1):(dim2+dim3)]
    b4<-bhat[(dim2+dim3+1):(dim2+dim3+dim4)]
    score<-fit.score(X1,X2,X3,X4,y,b2=b2,b3=b3,b4=b4)
  }else if(sum(y==0)>=1 & sum(y==1)<1){
    name<-c(rep("b1",dim1),"initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
    b1<-bhat[1:dim1];b3<-bhat[(dim1+1):(dim1+dim3)]
    b4<-bhat[(dim1+dim3+1):(dim1+dim3+dim4)]
    score<-fit.score(X1,X2,X3,X4,y,b1=b1,b3=b3,b4=b4)
  }else{
    name<-c(rep("b1",dim1),rep("b2",dim2),"initial.level","initial.trend","level.change","trend.change",rep("phi",dim4))
    b1<-bhat[1:dim1];b2<-bhat[(dim1+1):(dim1+dim2)]
    b3<-bhat[(dim1+dim2+1):(dim1+dim2+dim3)]
    b4<-bhat[(dim1+dim2+dim3+1):(dim1+dim2+dim3+dim4)]
    score<-fit.score(X1,X2,X3,X4,y,b1=b1,b2=b2,b3=b3,b4=b4)
  }
  Jhat<-HAC(score,L)/n
  Hhat<-fit[[3]]/n
  vcov<-solve(Hhat%*%solve(Jhat)%*%Hhat)
  ese<-sqrt(diag(vcov))
  p<-length(bhat)
  res<-matrix(0,p,5,dimnames=list(name, c('est','CI_lower','CI_upper','SE','p.value')))
  for (i in 1:p){
    res[i,1]<-signif(bhat[i],digits=4)
    res[i,2]<-signif(bhat[i]-qnorm(1-sig.level/2)*ese[i],digits=4)
    res[i,3]<-signif(bhat[i]+qnorm(1-sig.level/2)*ese[i],digits=4)
    res[i,4]<-signif(ese[i],digits=4)
    ts<-(bhat[i]^2)/(res[i,4]^2)
    res[i,5]<-signif(1-pchisq(ts,1),digits=4)
  }
  list(res,vcov)
}

#####Model selection for ITS#####
selection<-function(X1,X2,t,y,candidates,L=10){
  dim1<-ncol(X1);dim2<-ncol(X2);dim3<-ncol(X3);dim4<-ncol(X4);n<-length(y)
  bic<-c()
  for (inter in candidates) {
    X3<-cbind(rep(1,n),t,ifelse(t>inter,1,0),ifelse(t>inter,1,0)*(t-inter))
    X4<-cbind(rep(1,n),ifelse(t>inter,1,0))
    fiti<-est.ILL(X1,X2,X3,X4,y)
    bhat<-fiti[[1]]
    lli<-est.ILL(X1,X2,X3,X4,y)[[2]]
    if(sum(y==0)<1 & sum(y==1)<1){
      b3<-bhat[1:dim3];b4<-bhat[(dim3+1):(dim3+dim4)]
      score<-fit.score(X1,X2,X3,X4,y,b3=b3,b4=b4)
    }else if(sum(y==0)<1 & sum(y==1)>=1){
      b2<-bhat[1:dim2];b3<-bhat[(dim2+1):(dim2+dim3)]
      b4<-bhat[(dim2+dim3+1):(dim2+dim3+dim4)]
      score<-fit.score(X1,X2,X3,X4,y,b2=b2,b3=b3,b4=b4)
    }else if(sum(y==0)>=1 & sum(y==1)<1){
      b1<-bhat[1:dim1];b3<-bhat[(dim1+1):(dim1+dim3)]
      b4<-bhat[(dim1+dim3+1):(dim1+dim3+dim4)]
      score<-fit.score(X1,X2,X3,X4,y,b1=b1,b3=b3,b4=b4)
    }else{
      b1<-bhat[1:dim1];b2<-bhat[(dim1+1):(dim1+dim2)]
      b3<-bhat[(dim1+dim2+1):(dim1+dim2+dim3)]
      b4<-bhat[(dim1+dim2+dim3+1):(dim1+dim2+dim3+dim4)]
      score<-fit.score(X1,X2,X3,X4,y,b1=b1,b2=b2,b3=b3,b4=b4)
    }
    Jhat<-HAC(score,L)
    Hhat<-fiti[[3]]
    bici<-2*lli+log(n)*sum(diag(solve(Hhat)%*%Jhat))
    bic<-c(bic,bici)
  }
  list(candidates[which.min(bic)],bic)
}

