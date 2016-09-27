###
### first attempt of analysis
### 
### prerequisities
### - modify (or create) defaults-local.txt
###   - check defaults-global.dat for examples
###     - insert xx.momo.path=/path/to/a-momo-outputs as defined in MOMOmaster.do (global WDIR)
###     - OR save your mortality data into download
###       use name mort_xx.dat
###       variables "n" (number of deaths) and "date" (date of the monday of the week)
###       optionally "pop" for population, age groups etc
### - run mort-data.R   produces data/mort-raw.rda
### - run data.R        produces data/weather.rda
### - run indicators.R  produces data/mort.rda
###
source("fun/auxfun.R")
source("fun/data.R")
source("fun/analysis.R")
### load defaults so that no need for country specific code
loadDefaults()

### load the combined weather and mortality data
load("data/mort.rda")

### prediction data
summary(preddata<-data.frame(date=monday(Sys.Date()),temp=seq(-30,30,0.5)))
preddata$temp.adj1<-preddata$temp
preddata$temp.adj2<-preddata$temp

### crude fits - spline
preds<-fits<-list()

fits$crude.spline<-list()
for(i in names(morts)) {
    fits$crude.spline[[i]]<-glm(n~wlsb(temp,n=10,by=3,ref=0),data=morts[[i]],family=poisson)
}
lapply(fits$crude.spline,function(a) coef(summary(a))) # which are 'significant'
preds$crude.spline<-lapply(fits$crude.spline,predfun,preddata)

### crude fits - number
fits$crude.num<-list()
for(i in names(morts)) {
    fits$crude.num[[i]]<-glm(n~wsei(temp,n=10,by=3,ref=0),data=morts[[i]],family=poisson)
}
preds$crude.num<-lapply(fits$crude.num,predfun,preddata)

### plotting

png("out/eu_an_01.png",,width=1200,height=500*length(morts),pointsize=20)
par(mfrow=c(length(morts),1))
for(i in names(morts)) {
    ##png(paste("out/",i,"_an_01.png",sep=""),width=1200,height=800,pointsize=10)       
    with(morts [[i]],plot(temp,n,main=i,xlim=range(preddata$temp)))
    plot(preds$crude.spline[[i]],"temp","blue")
    plot(preds$crude.num   [[i]],"temp","green")
    ##dev.off()
}
dev.off()

### adjusted
###
for(i in names(morts)) { # all datasets
    for(j in c("temp")) { # all variables
        fit<-lm(formula(paste(j,"~momosin(date,2)")),data=morts[[i]])
        pred<-predict(fit,newdata=morts[[i]],se=TRUE)
        morts[[i]][[paste(j,".adj1",sep="")]]<-with(pred,morts[[i]][[j]]-fit)
        morts[[i]][[paste(j,".adj2",sep="")]]<-with(pred,(morts[[i]][[j]]-fit)/se.fit)
    }
}

png("out/eu_an_adj.png",,width=1200,height=500*length(morts),pointsize=20)
par(mfrow=c(length(morts),1))
for(i in names(morts)) {
    ##png(paste("out/",i,"_an_adj.png",sep=""),width=1200,height=800,pointsize=10)       
    with(morts [[i]],matplot(numdate(date),cbind(temp,temp.adj1,temp.adj2),type="l",lty=1,main=i))
    ##dev.off()
}
dev.off()

### adjusted fits - spline
fits$adj.spline<-list()
for(i in names(morts)) {
    fits$adj.spline[[i]]<-glm(n~momosin(date,2)+wlsb(temp.adj1,n=5,by=1,ref=0),data=morts[[i]],family=poisson)
}
lapply(fits$adj.spline,function(a) coef(summary(a)))
preds$adj.spline<-lapply(fits$adj.spline,predfun,preddata)

### crude fits - number
fits$adj.num   <-list()
for(i in names(morts)) {
    fits$adj.num   [[i]]<-glm(n~momosin(date,2)+wsei(temp.adj1,n=5,by=1,ref=0),data=morts[[i]],family=poisson)
}
preds$adj.num<-lapply(fits$adj.num,predfun,preddata)

### plotting
png("out/eu_an_02.png",,width=1200,height=500*length(morts),pointsize=20)
par(mfrow=c(length(morts),1))
for(i in names(morts)) {
    ##png(paste("out/",i,"_an_02.png",sep=""),width=1200,height=800,pointsize=10)       
    with(morts [[i]],plot(temp.adj1,n,main=i,xlim=c(-10,10)))
    plot(preds$adj.spline[[i]],"temp.adj1","blue")
    plot(preds$adj.num   [[i]],"temp.adj1","green")
    ##dev.off()
}
dev.off()

nnames<-function(a) structure(names(a),names=names(a))
pred.date<-lapply(fits,function(a) lapply(nnames(a),function(b) predfun(a[[b]],morts[[b]])))
lapply(pred.date,sapply,data.class)
colors<-rainbow(length(fits))
names(colors)<-names(fits)

png("out/eu_an_03.png",,width=1200,height=500*length(morts),pointsize=20)
par(mfrow=c(length(morts),1))
for(i in names(morts)) {
    ##png(paste("out/",i,"_an_03.png",sep=""),width=1200,height=800,pointsize=10)       
    with(morts [[i]],plot(date,n,main=i))
    for(j in names(pred.date)) plot(pred.date[[j]][[i]],"date",colors[j])
    legend("topright",col=colors,legend=names(colors),lty=1)
    ##dev.off()
}
dev.off()

save(morts,preddata,countries,file="data/mort.rda") # writes over
