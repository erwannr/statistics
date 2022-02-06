library(tseries)
library(moments)
# Import data
x<-scan("C:\\Documents\\work\\Teaching\\Fa06_TS\\tutorial\\d-hwp3dx8099.dat")
dim(x)
x<-matrix(x,ncol=4,byrow=T)
dim(x)
ex1<-list()
ex1$ts<-ts(x[,3],deltat=1/250,start=c(1980,1))
ex1$len<-length(ex1$ts)

#what is ts?
?ts
is(ex1$ts)
methods(ts)

ex1$sr<-exp(ex1$ts/250)
ex1$pr<-ex1$sr
ex1$pr[1]<-100
for(i in 2:ex1$len) ex1$pr[i]<-ex1$sr[i]*ex1$pr[i-1]

# Stationarity 
ts.plot(ex1$ts,main="lr",xlab="time",ylab="")
ts.plot(ex1$pr,main="pr",ylab="")

oct87<-1987+10/12+19/250
oct97<-1997+10/12+27/250
aug98<-1998+8/12+3/250

lines(c(oct87,oct87),c(140,200))
lines(c(oct97,oct97),c(500,600))
lines(c(aug98,aug98),c(500,600))

adf.test(ex1$pr) #make clear what is H0 and the outcome
adf.test(diff(ex1$pr,lag=1))
adf.test(ex1$ts)

#Order selection
train<-window(ex1$ts,start=c(1980,1),end=c(1989,250))
par(mfrow=c(2,1))
acf(train,ylim=c(-0.4,0.4),main="lr -- 1980-1989",lag.max=9)
acf(train,ylim=c(-0.4,0.4),main="",type="partial",lag.max=9)

aicTabl<-function(dta,parCnt){
	ar<-array(NA,dim=c(parCnt,parCnt),dimnames=c("p","q"))
	for(p in 1:parCnt){
	for(q in seq(1,parCnt+1-p)){
		tmp<-arima(dta,order=c(p,0,q))
		ar[p,q]<-tmp$aic
	}
	}
ar
}
ex1$aic<-aicTabl(train,4)

#Estimation
ex1$fit<-arima(train,c(3,0,2))
ex1$fit
names(ex1$fit)
ex1$fit$var.coef #let's compute the t-values
ex1$fit$tval<-ex1$fit$coef/sqrt(diag(ex1$fit$var.coef))
ex1$fit$tval

#Model checking
tsdiag(ex1$fit)

deflt<-par(no.readonly=TRUE)
par(mfcol=c(1,2))
hist(ex1$fit$resid,nclass="Scott",main="residuals",xlab="",xlim=c(-5,5))
tmp<-seq(-5,5,length=100)
lines(tmp,500*dnorm(tmp,mean=mean(ex1$fit$resid),sd=sd(ex1$fit$resid)))
text(-2.8,490,labels=paste("xsK=",round(kurtosis(ex1$fit$resid)-3,1)))
text(-2.8,450,labels=paste("skew=",round(skewness(ex1$fit$resid,1),1)))
par(pty="s")
qqnorm(ex1$fit$resid, main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles",
            ylab = "Sample Quantiles",xlim=c(-10,10),ylim=c(-10,10))
lines(c(-10,10),c(-10,10))

sq<-ex1$fit$resid^2;
sd<-sqrt(filter(sq,rep(1/10,10))-filter(ex1$fit$resid,rep(1/10,10))^2)
ts.plot(sd,ylab="sample sd of residuals",main=paste("filter size=",c(10)))
acf(ex1$fit$resid^2,main="squared residuals")

#Forecasting
ex1$pred<-predict(ex1$fit,n.ahead=10)
tmp<-ts.plot(window(ex1$ts,start=c(1989,240),end=c(1990,10)),ylim=c(-3,3),main="forecast",ylab="")
lines(ex1$pred$pred)
lines(ex1$pred$pred-ex1$pred$se*2,lty=2)
lines(ex1$pred$pred+ex1$pred$se*2,lty=2)
grid(20,10)

#Wrapping up
ls()
names(ex1)
write(ex1$aic,"C:\\Documents\\work\\Teaching\\Fa06_TS\\tutorial\\aic_table",ncolumns=4)

#Garch
library(fSeries)
ex2<-list()
ex2$fit<-garchFit(train,formula.mean=~arma(0,0),formula.var=~garch(1,1),trace=TRUE,iter=100,init.rec="uev")
ex2$e<-ex2$fit@residuals
ex2$h<-ex2$fit@fit$series$h
ex2$par<-ex2$fit@fit$params$params
ts.plot(ts(ex2$h,deltat=1/250,start=c(1980,1)),ylab="h (cond var)")
par(mfcol=c(2,1))
hist(ex2$e,nclass="Scott",main="",xlim=c(-10,10),xlab="e_t")
text(-9,300,labels=paste("xsK=",round(kurtosis(ex2$e)-3,1)))
text(-9,250,labels=paste("skew=",round(skewness(ex2$e,1),1)))
hist(ex2$e/sqrt(ex2$h),nclass="Scott",main="",xlab="z_t=e_t/sqrt(h_t)",xlim=c(-10,10))
text(-9,250,labels=paste("xsK=",round(kurtosis(ex2$e/sqrt(ex2$h))-3,1)))
text(-9,225,labels=paste("skew=",round(skewness(ex2$e/sqrt(ex2$h),1),1)))

