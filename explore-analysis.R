###
### exploring different analysis options!
### 
### 
source("fun/auxfun.R")
source("fun/data.R")
source("fun/analysis.R")

### You need to run analysis1.R first
load("data/mort.rda") # includes object morts, preddata, country

library("splines")

### Stepwise regression
### For proper analysis, we need to split the basis functions into separate variables
for(i in c("base",names(morts))) {
    if(i=="base") data<-preddata
    else          data<-morts[[i]]
    for(d in -2:3) {
        if(d==-2) lmat<-wlsb (data$temp,by=2,n=10,ref= 20 )
        if(d==-1) lmat<-wlsb (data$temp,by=2,n=10,ref=-20 )
        if(d== 0) lmat<-wlsb (data$temp,by=2,n=10         )
        if(d>  0) lmat<-wlsbs(data$temp,by=2,n=10,degree=d)
        for(j in dimnames(lmat)[[2]]) {
            if(d==-2) nam<-paste("temp.p.",    j,sep="")
            if(d==-1) nam<-paste("temp.m.",    j,sep="")
            if(d== 0) nam<-paste("temp.",      j,sep="")
            if(d>  0) nam<-paste("temp.",d,".",j,sep="")
            data[[nam]]<-na.0(lmat[,j])
        }
    }
    if(i=="base") preddata  <-data
    else          morts[[i]]<-data
}
### Did we dooorite?
names(morts[[countries[1]]])
summary(morts[[countries[1]]])

par(mfcol=c(1,1))
with(preddata,matplot(temp,cbind(temp.left9,temp.p.left9,temp.m.left9,temp.1.b1),type="l"))
with(preddata,matplot(temp,cbind(wlsb(temp,n=2,by=10,ref=0),wlsb(temp,n=2,by=10,ref=-10)),type="l",col=rep(1:2,c(4,4))))

bases<-list()
bases$base<-with(preddata,wlsb (temp,n=5,by=5))
for(i in paste(1:3)) bases[[paste("bs",i,sep="")]]<-with(preddata,wlsbs(temp,n=5,by=5,degree=i))
par(mfcol=c(4,1))
for(i in names(bases)) with(preddata,matplot(temp,bases[[i]],main=i,type="l",lty=1,col=1:2,ylim=c(0,1)))
args(cor)
sapply(bases,cor,use="pair")

fit0a<-list()
for(i in names(morts)) {
    fit0a[[i]]<-glm(formula(paste("n~wlsb(temp)",sep="")),data=morts[[i]],family=poisson)
}
fit0b<-list()
for(i in names(morts)) {
    fit0b[[i]]<-glm(formula(paste("n~wlsbs(temp)",sep="")),data=morts[[i]],family=poisson)
}
matplot(sapply(fit0a,predict),sapply(fit0b,predict))

names(morts[[country]])

mods<-list(base0=paste("n~",paste("temp.left"  ,9:0,sep="",collapse="+"),"+",paste("temp.right"  ,0:9,sep="",collapse="+")),
           basep=paste("n~",paste("temp.p.left",9:0,sep="",collapse="+"),"+",paste("temp.p.right",0:9,sep="",collapse="+")),
           basem=paste("n~",paste("temp.m.left",9:0,sep="",collapse="+"),"+",paste("temp.m.right",0:9,sep="",collapse="+")),
           bs1=paste("n~",paste("temp.1.b",1:20,sep="",collapse="+")),
           bs2=paste("n~",paste("temp.2.b",1:22,sep="",collapse="+")),
           bs3=paste("n~",paste("temp.3.b",1:24,sep="",collapse="+")))
### Basis fits
all.cols<-unique(unlist(sapply(morts,names)))
common.cols<-all.cols[apply(sapply(morts,function(a) all.cols%in%names(a)),1,all)]
mortdata<-do.call("rbind",lapply(names(morts),function(a) cbind(morts[[a]][common.cols],country=a)))

fit1<-list()
fit1$base<-list()
for(i in names(morts)) {
    fit1$base[[i]]<-list()
    for(j in names(mods)) {
        fit1$base[[i]][[j]]<-eval(substitute(glm(formula(mods[[j]]),data=morts[[i]],family=poisson)))
    }
}
names(fit1$base)
(function(a) a[order(a[,4]),])(coef(summary(fit1[[1]])))

library("MASS") # for some utilities

### We will penalize the model complexity with relatively high number
### to get minimal number of covariates
log(sapply(morts,nrow)) # 5.7 for Finland
sapply(morts,function(a) mean(a$n))

### Backward selection
fit1$drops<-list()
for(i in names(morts)) {
    fit1$drops[[i]]<-list()
    for(j in names(mods)) {
        fit1$drops[[i]][[j]]<-step(fit1$base[[i]][[j]],test="Chi",k=6,direction="backward",keep=function(a,b) predfun(a,preddata))
    }
}
drops1[[1]]$anova

### Forward selection 
fit1$adds<-list()
for(i in names(morts)) {
    fit1$adds[[i]]<-list()
    for(j in names(mods)) {
        fit1$adds[[i]][[j]]<-step(update(fit1$base[[i]][[j]],n~1),
                                  list(upper=formula(fit1$base[[i]][[j]]),
                                       lower=n~1),
                                  test="Chi",k=2,direction="forward",keep=function(a,b) predfun(a,preddata))
    }
}

### Both forward and backward selection 
fit1$both<-list()
for(i in names(morts)) {
    fit1$both[[i]]<-list()
    for(j in names(mods)) {
        fit1$both[[i]][[j]]<-step(update(fit1$base[[i]][[j]],n~1),
                                  list(upper=formula(fit1$base[[i]][[j]]),
                                       lower=n~1),
                                  test="Chi",k=6,direction="both",keep=function(a,b) predfun(a,preddata))
    }
}

par(mfrow=c(length(mods),3*length(morts)),mar=c(0,0,2,0),oma=c(0,0,0,0))
for(model in names(mods)) {
    for(co in names(morts)) {
        for(j in names(fit1)[-1]) {
            with(morts[[co]],plot(temp,n,xaxt="n",yaxt="n",main=paste(co,j,model)))
            for(i in 1:ncol(fit1[[j]][[co]][[model]]$keep)) plot.pred.frame(fit1[[j]][[co]][[model]]$keep[,i],"temp",col=i,alpha=0.1,lwd=3)
        }
    }
    cat("\n")
}

################################################################################
### LASSO: idea is to penalize the regression coefficients to find a minimal model
### LASSO penalizes the absolute values of the regression coefficients. This is equivalent to
###  assuming Laplacian (or double exponential) prior. In practise, what it does, it
###  forces some of  the coefficients to zero. In some sense it is a bit like backward selection
library("glmnet")
lasso1 <-lapply(morts,function(a) with(a,   glmnet(cbind(1,wlsb(temp,by=1,n=20)),n,family="poisson")))
lasso1c<-lapply(morts,function(a) with(a,cv.glmnet(cbind(1,wlsb(temp,by=1,n=20)),n,family="poisson")))
par(mfcol=c(length(morts),1),mar=c(2,2,2,2))
### numbers top of the picture are effective numbers of parameters
for(i in names(lasso1)) plot(lasso1c[[i]],main=i)
for(i in names(lasso1)) plot(lasso1c[[i]]$glmnet.fit)
for(i in names(lasso1)) plot(lasso1c[[i]]$glmnet.fit,"lam")

for(i in morts) with(i, plot(temp,n,xlim=c(-30,30)))

### Predictions (no SE:s available)
lasso1pred <-exp(predict(lasso1c$glmnet.fit,with(preddata,cbind(1,wlsb(temp,by=1,n=20)))))
lasso1cpred<-exp(predict(lasso1c           ,with(preddata,cbind(1,wlsb(temp,by=1,n=20)))))

### Pics or it didn't happen
with(morts[[country]],plot(temp,n,pch=16))
matlines(preddata$temp,lasso1pred ,type="l",lty=1,col=rainbow(100))
matlines(preddata$temp,lasso1cpred,type="l",lty=1,col="black",lwd=3)
###
### The problem with LASSO is that it is considered biased. Most methods for calculating 
### standard errors (and thus CIs) assume that the estimate is unbiased so that the
### variances can be estimated from either Fisher or observed information. The
### bias makes this unusable

################################################################################
###
### Other implementation of LASSO, that includes also Ridge Regression (RR)
### Does not seem to be very stable. Needs a bit of rework. Nevertheless,
### RR puts penalty to the sum of squares of the parameters. It puts equal 
###  penalty to all of the coefficients, bringing them closer to zero.
### Package penalized implements both penalties simultaneously. It would
### be nice to be able to use either. Does not seem to work.
library("penalized")
pen1<-with(morts[[country]],penalized(n~wlsb(temp,by=1,n=20),lambda1=100,lambda2=100,model="poisson",steps=1)) # sloooooooooooow
plot(pen1)
coef(pen1,"all")
with(morts[[country]],plot(temp,n))
with(preddata,lines(temp,predict(pen1,wlsb(temp,by=1,n=20))))

################################################################################
### Fullblown GAM
### These are just silly trials.
library("mgcv")
gams<-list()
for(i in names(morts)) {
    gams[[i]]<-list()
    for(k in c(10,20)) {
        for(adj in c("","+momosin(date,2)")) {
            for(va in c("","temp","temp.adj1")) {
                for(bs in c("ad")) {
                    if(va!="")
                        fo<-paste("n~1",adj,"+s(",va,",k=",k,",bs=\"",bs,"\")",sep="")
                    else
                        fo<-paste("n~1",adj,                                   sep="")
                    print(fo)
                    na<-paste(adj,bs,k,va)
                    gams[[i]][[na]]<-gam(formula(fo),data=morts[[i]],family=poisson)
                }
            }
        }
    }
}
pen.edf(gams[[1]][[1]])
sapply(gams,sapply,function(a) min(pen.edf(a)))

names(gams[[1]])
par(mfcol=c(length(gams[[1]])/2,length(gams)*2))
for(i in names(gams)) {
    for(m in names(gams[[i]])) {
        print(paste(i,m))
        plot(gams[[i]][[m]],residuals=FALSE,all.terms=FALSE,shade=TRUE,seWithMean=TRUE,trans=exp,main=paste(i,m))
    }
}

### there is a correlation between date and temp so that the default value for date is missleading
### unfortunately this correlation is different for each country
summary(preddata$date) 

### effect predictions
par(mfcol=c(6,4))
for(i in names(morts)) {
    for(j in grep("ad 10",names(gams[[i]]),value=TRUE)) {
        vari<-ifelse(grepl("temp.adj",j),"temp.adj1","temp")
        with(morts[[i]],plot(get(vari),n,main=paste(i,j),xlim=list("temp"=c(-25,25),"temp.adj1"=c(-15,10))[[vari]]))
        with(predict(gams[[i]][[j]],newdata=preddata,se=TRUE,type="response"),{
            lines(preddata$temp,fit,type="l")
            polygon(c(preddata$temp,rev(preddata$temp)),c(fit-2*se.fit,rev(fit+2*se.fit)),border=NA,col=rgb(0,0,0,.2))
        })
             
    }
}

gamsf<-lapply(gams,function(a) sapply(a[grep("ad 10",names(a))],fitted))
head(gamsf[[1]])
par(mfcol=c(length(gamsf),1))
for(i in names(morts)) {
    with(morts[[i]],{
        matplot(temp,n-cbind(gamsf[[i]]),main=i)
    })
}
pairs(morts[[1]]$n-gamsf[[1]])
    
sapply(gamsf,matplot)

### clearly the effect is modified by the presence of the seasonal pattern
mcols<-rainbow(12)
mshad<-rainbow(12,alpha=.2)
Mcols<-rainbow( 4)
Mshad<-rainbow( 4,alpha=.2)
names(Mcols)<-names(Mshad)<-names(gams[[1]])
par(mfcol=c(12,length(morts)),mar=c(0,0,2,0))
for(i in names(morts)) {
    print(yl<-range(morts[[i]]$n))
    print(xl<-range(morts[[i]]$temp.adj1))
    ##with(morts[[i]],plot(temp.adj1,n,main=paste(i),xlim=xl,ylim=yl,yaxt="n",xaxt="n",col=mcols[as.numeric(format(date,"%m"))]))
    for(m in 1:12) {
        with(subset(morts[[i]],as.numeric(format(date,"%m"))==m),plot(temp.adj1,n,main=paste(i,m),xlim=xl,ylim=yl,yaxt="n",xaxt="n"))
        tmppred<-preddata
        tmppred$date<-as.Date(paste("2014-",m,"-15",sep=""))
        for(j in c(" ad 10 temp "," ad 20 temp ")) {
            pred<-predict(gams[[i]][[j]],newdata=tmppred,se=TRUE,type="response")
            with(pred,{
                lines(preddata$temp,fit,type="l",col=Mcols[j])
                polygon(c(preddata$temp,rev(preddata$temp)),c(fit-2*se.fit,rev(fit+2*se.fit)),border=NA,col=Mshad[j])
            })
        }
    }
}
sapply(morts,function(a) table(format(a$date,"%m")))

zscore<-function(fit,y) {
    with(predict(fit,se=TRUE,type="response"),{
        tse<-fit+se.fit^2
        predse<-sqrt((4/9)*fit^(-2/3)*tse)
        (y^(2/3)-fit^(2/3))/predse
    })
}
?predict.gam
### zscores by model
gamz<-lapply(names(gams),function(a) sapply(gams[[a]],zscore,morts[[a]]$n))
### raw residual as % of observed
gamRO<-lapply(names(gams),function(a) 100*(1-sapply(gams[[a]],predict,type="response")/morts[[a]]$n))
### raw residual as % of expected
gamRE<-lapply(names(gams),function(a) 100*(-1+morts[[a]]$n/sapply(gams[[a]],predict,type="response")))
### raw residual as % of mean
gamRM<-lapply(names(gams),function(a) 100*(morts[[a]]$n-sapply(gams[[a]],predict,type="response"))/mean(morts[[a]]$n))
###
names(gamRO)<-names(gamRE)<-names(gamRM)<-names(gamz)<-names(gams)
head(gamRM[[1]])
varz<-c(" ad 10 ","+momosin(date,2) ad 10 ","+momosin(date,2) ad 10 temp.adj1"," ad 10 temp")
varz%in%names(gams[[1]])
par(mfrow=c(length(morts),length(varz)))
for(i in 1:length(gamz)) for(j in 1:length(varz)) {
    if(j<length(varz))
    matplot(gamz[[i]][,varz[j:length(varz)]]-gamz[[i]][,varz[j]],type="l",lty=1,ylim=c(-10,10),col=j:length(varz))
    else
    matplot(gamz[[i]][,varz[length(varz)]]-gamz[[i]][,varz[2]],type="l",lty=1,ylim=c(-10,10),col=j:length(varz))
}

do.call("rbind",lapply(gams,function(a) anova(a[[varz[1]]],a[[varz[2]]],a[[varz[3]]],a[[varz[4]]],test="Ch")))
sapply(gams,function(a) sapply(a[varz],AIC)) # zmaller iz better

par(mfrow=c(length(morts),length(varz)))
for(i in names(gamz)) for(j in 1:length(varz)) {
    x<-morts[[i]]$temp
    if(j>=1)        y<-gamRM[[i]][,varz[j]]
    #if(j%in%c(2,3)) y<-gamz[[i]][,varz[j]]-gamz[[i]][,varz[j-1]]
    #if(j%in%c(4  )) y<-gamz[[i]][,varz[4]]-gamz[[i]][,varz[  2]]
    y<-y[order(x)]
    x<-x[order(x)]
    plot(x,y,main=paste(i,j),ylim=c(-30,30),
         xlim=c(-30,30),xaxs="i")
    abline(h=0)
    xx<--30:30
    with(predict(gam(y~s(x)),newdata=data.frame(x=xx),se=TRUE),{
        lines(xx,fit)
        polygon(c(xx,rev(xx)),c(fit-2*se.fit,rev(fit+2*se.fit)),col=rgb(0,0,0,.3),border=NA)
    })
}

with(predict(gam(y~x),se=TRUE),cbind(fit,se.fit))

### not working, yet
modeffects(gams[[1]][[2]],morts[[1]])

################################################################################
###
### gamlss has a wide selection of smoothers also
require("gamlss")

gamlsss<-list()
for(co in sample(names(morts))) {
    gamlsss[[co]]<-list()
    for(ef in c("pb","cs","fp","lo")) {
        for(va in c("temp","temp.adj1")) {
            for(ad in c("","momosin(date,3)+")) {
            
                print(fo<-paste("n~",ad,"",ef,"(",ifelse(ef=="lo","~",""),va,")",sep=""))
                res<-try(gamlss(as.formula(fo),family=PO(),data=morts[[co]][,c("n",va,"date")]))
                                        #if(!inherits(res,"try-error"))
                gamlsss[[co]][[paste(ef,va,ad)]]<-res
            }
        }
    }
}
sapply(gamlsss,names)
sapply(gamlsss,sapply,data.class)
sapply(gamlsss,sapply,summary)
par(mfcol=c(length(morts),1))
lapply(lapply(gamlsss,sapply,predict),matplot,type="l",lty=1)
par(mfcol=c(3*2*2,3))
for(co in names(morts)) {
    for(ef in c("pb","cs","fp")) {
        for(va in c("temp","temp.adj1")) {
            for(ad in c("","momosin(date,3)+")) {
                try(term.plot(gamlsss[[co]][[paste(ef,va,ad)]],partial.resid=TRUE,main=paste(co,ef,va,ad),what="mu",terms=ifelse(ad=="",1,2)))
            }
        }
    }
}
par(mfcol=c(3*2*1,3))
for(co in names(morts)) {
    for(ef in c("pb","cs","fp")) {
        for(va in c("temp","temp.adj1")) {
            for(ad in c("momosin(date,3)+")) {
                try(term.plot(gamlsss[[co]][[paste(ef,va,ad)]],partial.resid=TRUE,main=paste(co,ef,va,ad),what="mu",terms=1))
            }
        }
    }
}


summary(morts[[2]])

table(floor(((1:12)%%12)/3),1:12)
################################################################################
###
### analysis by season
#ggg<-gam(n~momosin(date,4)+s(time,bs="cr",k=4)+s(temp,by=monthf,k=1),data=morts$fi,family=poisson)
tpred<-mkpreddata(weeks,morts$fi)
xyplot(temp.adj1~seas,groups=temp,data=tpred,type="l")
ggg<-gam(n~momosin(date,4)+s(time,bs="cr",k=4)+s(temp.adj1,by=factor(floor((month%%12)/3)),k=1,fx=FALSE,bs="cr"),data=morts$fi,family=poisson)
ggG<-gam(n~momosin(date,4)+s(time,bs="cr",k=4)+s(temp.adj1,seas,bs="tp"),data=morts$fi,family=poisson)
gGg<-gam(n~momosin(date,4)+s(time,bs="cr",k=4)+s(temp     ,by=factor(floor((month%%12)/3)),k=1,fx=FALSE,bs="cr"),data=morts$fi,family=poisson)
gGG<-gam(n~momosin(date,4)+s(time,bs="cr",k=4)+s(temp     ,seas,bs="tp"),data=morts$fi,family=poisson)
tpredy<-predict(ggg,type="response",newdata=tpred)
tpredY<-predict(ggG,type="response",newdata=tpred)
tpreDy<-predict(gGg,type="response",newdata=tpred)
tpreDY<-predict(gGG,type="response",newdata=tpred)
predm<-with(tpred,tapply(tpredy,list(seas*366,temp),mean))
predM<-with(tpred,tapply(tpredY,list(seas*366,temp),mean))
preDm<-with(tpred,tapply(tpreDy,list(seas*366,temp),mean))
preDM<-with(tpred,tapply(tpreDY,list(seas*366,temp),mean))
morts$fi$mort<-with(morts$fi,1000*n/pop)
par(mfcol=c(2,2))
contour(weeks,-30:30,predm      ,zlim=c( 800,1200),nlevels=100,col=colorRampPalette(c("red","blue"))(100),labels="")
with(morts$fi,symbols(seas*366,temp,circles=sqrt(100*mort),add=TRUE,inches=FALSE))
contour(weeks,-30:30,      predM,zlim=c( 800,1200),nlevels=100,col=colorRampPalette(c("red","blue"))(100),labels="")
with(morts$fi,symbols(seas*366,temp,circles=sqrt(100*mort),add=TRUE,inches=FALSE))
contour(weeks,-30:30,preDm      ,zlim=c( 800,1200),nlevels=100,col=colorRampPalette(c("red","blue"))(100),labels="")
with(morts$fi,symbols(seas*366,temp,circles=sqrt(100*mort),add=TRUE,inches=FALSE))
contour(weeks,-30:30,      preDM,zlim=c( 800,1200),nlevels=100,col=colorRampPalette(c("red","blue"))(100),labels="")
with(morts$fi,symbols(seas*366,temp,circles=sqrt(100*mort),add=TRUE,inches=FALSE))

contour(weeks,-30:30,predm-predM,zlim=c(-200, 200),nlevels=100,col=colorRampPalette(c("red","blue"))(100),labels="")
with(morts$fi,symbols(seas*366,temp,circles=sqrt(100*mort),add=TRUE,inches=FALSE))
with(morts$fi,summary(sqrt(n/100)))
pen.edf(ggg)
par(mfcol=c(5,2))
plot(gGg,shade=TRUE)
plot(ggg,shade=TRUE)
summary(ggg)
dev.off()
plot(gamfits$fi$MTB)
?contour

par(mfcol=c(4,4))
plot(gamfits$fi$MTB,shade=TRUE)
################################################################################
###
### Search algorithms?
### 
