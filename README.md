# JamesHotniel.github.io
OPTIMIZATION OF BINARY LOGISTIC REGRESSION PARAMETERS USING THE BROYDEN-FLETCHER-GOLDFARB-SHANNO ALGORITHM IN THE QUASI-NEWTON METHOD
The following is the syntax of the BFGS algorithm in Binary Logistic Regression
(Case Study: Study Period Status of FMIPA Undergraduate Program Graduates Mulawarman University in 2021)
The research variables consist of a dependent variable (Y) and an independent variable (X). The dependent variable is the study period statistics of graduates of the FMIPA Mulawarman University Undergraduate Program in 2021 (Y). Independent variables include Grade Point Average (GPA) (X1), gender (X2), TOEFL score (X3), study program (X4), high school origin status (X5), region of origin (X6), and age (X7) . Variable notation, data types, operational definitions of variables.
The binary logistic regression parameter estimates were calculated with the help of R software using the "maxLik" package.
Syntax:

##Estimasi Parameter
> library("maxLik")
#Fungsi Likelihood
>  ll<- function(par){
+      y<- as.vector(data$Y)
+     x<- as.matrix(cbind(1, data$X1, data$X2, data$X3, data$X5))
+     n = length(y)
+     m = length(par)
+     beta <- par[1:5]   
+     loglik = rep(0,n)
+     for(i in 1:n){
+     xbeta= as.numeric(x[i,]%*%beta)
+     yd = y[i]*xbeta
+     loglik[i]=yd-log(1+exp(xbeta))
+     }
+     return(loglik)
+     }
#Gradien
>     gl<- function(par){
+     y<- as.matrix(data$Y)
+     x <- as.matrix(cbind(1, data$X1, data$X2, data$X3, data$X5))
+     n = length(y)
+     m = length(par)
+     beta <- par[1:5]
+     gg <- matrix(0,n,m)
+     p<- matrix(n,1)
+     xbeta<- matrix(n,1)
+     for(i in 1:n){
+     for(j in 1:m){
+       xbeta[i] <- as.numeric(x[i,]%*%beta)
+       p[i]<- exp(xbeta[i])/(1+exp(xbeta[i]))
+     
+       gg[i,j] <- x[i,j]%*%(y[i]-p[i])
+     }}
+     return(gg)
+     }

#MaxLikBFGS
>     sv<- c(Intercept=0, B1=0, B2=0, B3=0, B4=0)
>     max <- maxControl(tol=1e-8,print.level=3,iterlim=200)
>     mle<-maxLik(logLik=ll, grad = gl, start= sv, method="BFGS", control=max)

##Uji Simultan
>     library("MASS")
>     llfull<-as.matrix(ll(c(coef(mle))))
>     llreduced<-as.matrix(ll(c(0,0,0,0,0,0,0,0,0,0)))
>     model_full <- as.numeric(colSums(llfull))
>     model_full
>     model_reduced <-as.numeric(colSums(llreduced))
>     model_reduced
>     LRT <- -2*(model_reduced-model_full)
>     LRT
>     p.val <-pchisq(LRT, df = 9, lower.tail=FALSE)
>     p.val

##Uji Parsial
>     BFGS=coef(mle)
>     BFGS
>     se=stdEr(mle)
>     W<-BFGS/se
>     W

##Odd Ratio
>     OR1<-exp(Ob1)*1
>     OR1

