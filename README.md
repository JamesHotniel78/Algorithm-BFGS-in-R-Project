# JamesHotniel.github.io
>  OPTIMIZATION OF BINARY LOGISTIC REGRESSION PARAMETERS USING THE BROYDEN-FLETCHER-GOLDFARB-SHANNO ALGORITHM IN THE QUASI-NEWTON METHOD
>  The following is the syntax of the BFGS algorithm in Binary Logistic Regression
>  Case Study: Study Period Status of FMIPA Undergraduate Program Graduates Mulawarman University in 2021
>  The research variables consist of a dependent variable (Y) and an independent variable (X). The dependent variable is the study period statistics of graduates of the FMIPA Mulawarman University Undergraduate Program in 2021 (Y). Independent variables include Grade Point Average (GPA) (X1), gender (X2), TOEFL score (X3), study program (X4), high school origin status (X5), region of origin (X6), and age (X7) . Variable notation, data types, operational definitions of variables.
>  The binary logistic regression parameter estimates were calculated with the help of R software using the "maxLik" package.
>  Syntax:

##Statistic Deskription
>library("pastecs")
>library("dplyr")
>data=read.csv(file.choose(), header=TRUE, sep= ";")
>View(data)
>names(data)
>data$X2<-as.factor(data$X2)
>data$X4<-as.factor(data$X4)
>data$X41<-ifelse(data$X4==1,1,0)
>data$X41<-as.factor(data$X41)
>data$X42<-ifelse(data$X4==2,1,0)
>data$X42<-as.factor(data$X42)
>data$X43<-ifelse(data$X4==3,1,0)
>data$X43<-as.factor(data$X43)
>data$X5<-as.factor(data$X5)
>data$X6<-as.factor(data$X6)
>str(data)
>attach(data)
>summary(data$X1)
>summary(data$X3)
>summary(data$X7)

##Multikolinieritas Detection
>library("car")
>reg=lm(V~X1+X3+X7,data=data)
>vif(reg)
#Tabulasi X2 dengan X41
>tab.x2.x4=table(data$X2,data$X4)
>tab.x2.x4
>ind.test.x2.x4=chisq.test(tab.x2.x4)
>ind.test.x2.x4
>#Tabulasi X2 dengan X5
>tab.x2.x5=table(data$X2,data$X5)
>tab.x2.x5
>ind.test.x2.x5=chisq.test(tab.x2.x5)
>ind.test.x2.x5
#Tabulasi X2 dengan X6
>tab.x2.x6=table(data$X2,data$X6)
>tab.x2.x6
>ind.test.x2.x6=chisq.test(tab.x2.x6)
>ind.test.x2.x6
#Tabulasi X4 dengan X5
>tab.x4.x5=table(data$X4,data$X5)
>tab.x4.x5
>ind.test.x4.x5=chisq.test(tab.x4.x5)
>ind.test.x4.x5
#Tabulasi X4 dengan X6
>tab.x4.x6=table(data$X6,data$X4)
>tab.x4.x6
>ind.test.x4.x6=chisq.test(tab.x4.x6)
>ind.test.x4.x6
#Tabulasi X5 dengan X6
tab.x5.x6=table(data$X5,data$X6)
>tab.x5.x6
>ind.test.x5.x6=chisq.test(tab.x5.x6)
>ind.test.x5.x6

##Parameter Estimation
>library("maxLik")
>library("optimx")
>library("pscl")

#Likelihood Function
>ll<- function(par){
  y<- as.vector(data$Y)
  x<- as.matrix(cbind(1, data$X1, data$X2, data$X3, data$X5, data$X6, data$X7))
  beta <- par[1:7]
  m = length(par)
  n = length(y)
  loglik = rep(0,n)
  for(i in 1:n){
    xbeta= as.numeric(x[i,]%*%beta)
    yd = y[i]*xbeta
    loglik[i]=yd-log(1+exp(xbeta))
  }
  return(loglik)
}  

#Gradien
>gl<- function(par){
  y<- as.matrix(data$Y)
  x <- as.matrix(cbind(1, data$X1, data$X2, data$X3, data$X5, data$X6, data$X7))
  beta <- par[1:7]
  n = length(y)
  m = length(par)
  gg <- matrix(0,n,m)
  p<- matrix(n,1)
  xbeta<- matrix(n,1)
  for(i in 1:n){
    for(j in 1:m){
      xbeta[i] <- as.numeric(x[i,]%*%beta)
      p[i]<- exp(xbeta[i])/(1+exp(xbeta[i]))
    gg[i,j] <- x[i,j]%*%(y[i]-p[i])
    }
  }
  return(gg)
}

#MaxLikBFGS
>sv<- c(Intercept=0, B1=0, B2=0, B3=0, B5=0, B6=0, B7=0)
>max <- maxControl(tol=1e-3,print.level=3,iterlim=200)
>mle<-maxLik(logLik=ll, grad = gl, start= sv, method="BFGS", control=max)
>mle
>gradient(mle)
>hessian(mle)

##Simultan Test
>library("MASS")
>llfull<-as.matrix(ll(c(coef(mle))))
>llreduced<-as.matrix(ll(c(0,0,0,0,0,0,0,0,0,0)))
>model_full <- as.numeric(colSums(llfull))
>model_full
>model_reduced <-as.numeric(colSums(llreduced))
>model_reduced
>LRT <- -2*(model_reduced-model_full)
>LRT
>p.val <-pchisq(LRT, df = 9, lower.tail=FALSE)
>p.val

##Parsial Test
>BFGS=coef(mle)
>BFGS
>se=stdEr(mle)
>W<-BFGS/se
>W

##Odd Ratio
>Ob1<-4.4899213
>Ob7<--1.9893193
#Odd Ratio B1
>OR1<-exp(Ob1)*0.01
>OR1
#Odd Ratio B7
>OR7<-exp(Ob7)*1
>OR7
