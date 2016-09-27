###
### second attempt of analysis
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
library("mgcv")
library("lattice")
### load defaults so that no need for country specific code
loadDefaults()

### load the combined weather and mortality data
load("data/mort.rda")

for(i in names(morts)) morts[[i]]$seas<-seasday(morts[[i]]$date)
for(i in names(morts)) morts[[i]]$time<-numdate(morts[[i]]$date)
for(i in names(morts)) morts[[i]]$month<-as.numeric(format(morts[[i]]$date,"%m"))
for(i in names(morts)) morts[[i]]$monthf<-factor(morts[[i]]$month)
for(i in names(morts)) morts[[i]]$pop   <-rep(max(subset(sites,NUTS0b==toupper(i))$pop0),nrow(morts[[i]]))
for(i in names(morts)) morts[[i]]$country<-rep(i,nrow(morts[[i]]))
for(i in names(morts)) morts[[i]]$temp.adj3<-with(morts[[i]],residuals(gam(temp~s(seas,bs="cc"))))
for(i in names(morts)) morts[[i]]$at.adj3  <-with(morts[[i]],residuals(gam(at  ~s(seas,bs="cc"))))

### four main models
### T0: just the temperature
### TS: temperature, seasonality and trend
### TR: deseasonalized temperature, seasonality and trend
### ST: temperature by season, seasonality and trend
### 
mods<-list(nTn=n~1,
           MTn=n~                               momosin(date,4)+s(time,bs="cs",k=3),
           nRn=NULL,
           MRn=NULL,
           nTS=n~s(temp     ,bs="cs"          ),
           MTS=n~s(temp     ,bs="cs"          )+momosin(date,4)+s(time,bs="cs",k=3),
           nRS=NULL, # does not make sense
           MRS=n~s(temp.adj1,bs="cs"          )+momosin(date,4)+s(time,bs="cs",k=3), 
           nTB=n~s(temp     ,bs="cs",by=monthf)+momosin(date,4)+s(time,bs="cs",k=3),
           MTB=n~s(temp     ,bs="cs",by=monthf)+momosin(date,4)+s(time,bs="cs",k=3),
           nRB=n~s(temp.adj1,bs="cs",by=monthf)+momosin(date,4)+s(time,bs="cs",k=3),
           MRB=n~s(temp.adj1,bs="cs",by=monthf)+momosin(date,4)+s(time,bs="cs",k=3),
           nTI=n~  temp     *           monthf,
           MTI=n~  temp     *           monthf +momosin(date,4)+s(time,bs="cs",k=3),
           nRI=n~  temp.adj1*           monthf,
           MRI=n~  temp.adj1*           monthf +momosin(date,4)+s(time,bs="cs",k=3),
           nTX=n~s(temp     ,seas,bs="tp"     ),
           MTX=n~s(temp     ,seas,bs="tp"     )+momosin(date,4)+s(time,bs="cs",k=3),
           nRX=n~s(temp.adj1,seas,bs="tp"     ),
           MRX=n~s(temp.adj1,seas,bs="tp"     )+momosin(date,4)+s(time,bs="cs",k=3))

### Fit the models
gamfits<-list()
for(i in names(morts)) {
    gamfits[[i]]<-list()
    for(j in names(mods)) {
        cat(i,j)
        if(!is.null(mods[[j]])) {
            tmp<-try(gam(mods[[j]],data=morts[[i]],family=poisson))
            if(!inherits(tmp,"try-error")) {
                gamfits[[i]][[j]]<-tmp
                cat(" OK\n")
            } else {
                cat(" FAIL\n")          
            }
        } else {
            cat(" NULL\n")
        }
    }
}
par(mfcol=c(4,4),mar=c(3,3,1,1))
plot(gamfits[[1]]$MTB,shade=TRUE)
plot(gamfits[[1]]$MRB,shade=TRUE)
### simple predctions to see if any of the models fail spectacularly
png(file=paste("out/eu_an2_predicts.png",sep=""),width=1000,height=250*length(gamfits))
par(mfcol=c(length(gamfits),1),mar=c(2,2,2,1))
for(i in names(morts)) {
  matplot(morts[[i]]$time,sapply(gamfits[[i]],predict,type="response"),type="l",lty=1,xaxt="n",xaxs="i",main=i,
          ylab="", col=colorRampPalette(c("red","green","blue"))(length(mods)),lwd=1)
}
dev.off()
#par(mfcol=c(7,4),mar=c(0,0,2,0))
#for(i in names(gamfits)) for(j in grep("T",names(gamfits[[i]]),value=TRUE)) try(plot(gamfits[[i]][[j]],shade=TRUE,main=paste(i,j),select=1))
#for(i in names(gamfits)) for(j in grep("A",names(gamfits[[i]]),value=TRUE)) try(plot(gamfits[[i]][[j]],shade=TRUE,main=paste(i,j),select=1))
#par(mfcol=c(6,4))
#for(i in names(gamfits)) for(j in grep("T",names(gamfits[[i]]),value=TRUE)) try(plot(gamfits[[i]][[j]],shade=TRUE,main=paste(i,j),select=2))
#for(i in names(gamfits)) for(j in grep("A",names(gamfits[[i]]),value=TRUE)) try(plot(gamfits[[i]][[j]],shade=TRUE,main=paste(i,j),select=2))
#par(mfcol=c(3,4))
#for(i in names(gamfits)) for(j in grep("T",names(gamfits[[i]]),value=TRUE)) try(plot(gamfits[[i]][[j]],shade=TRUE,main=paste(i,j),select=3))
#for(i in names(gamfits)) for(j in grep("A",names(gamfits[[i]]),value=TRUE)) try(plot(gamfits[[i]][[j]],shade=TRUE,main=paste(i,j),select=3))
#lapply(gamfits,lapply,plot,shade=TRUE)
#?smooth.terms
#?plot.gam

sapply(morts,function(a) summary(a$seas*366)) ### chk the scale. not works if different in preddata
sapply(morts,function(a) min(format(a$date,"%j"))) # sumthin pekuliarr
### Create a large matrix of time x temp values. time by day, temp by 1°C
weeks<-seq(2,366,by=7) # need to use daily data, as "seas" may differe for week "03" each year
### For simplicity assume constant calndartime. This is not right as there's a trend in the models
###  and the trend might be correlated with the seasonal terms. In final analyses, this must be fixed
###  somehow.
preds<-list()
for(i in names(gamfits)) {
  preds[[i]]<-list()
  preddata<-mkpreddata(weeks,morts[[i]])
  for(j in names(gamfits[[i]])) {
    cat(i,j,"\n")
    tmp<-try(predict(gamfits[[i]][[j]],type="response",newdata=preddata))
    if(!inherits(tmp,"try-error"))
      preds[[i]][[j]]<-with(preddata,tapply(tmp,list(seas*366,temp),mean))      
  }
}
#with(subset(preddata,temp==2),plot(seas,predtemp2,type="l"))
#(seq(2,366,by=7)%in%(floor((2:366)/7)*7+2))
### observen mean mortality using the same grid. Mostly missing values! (not used)
obs<-lapply(morts,function(a) with(a,(tapply(n,list(factor(round(seas*366/7)*7+2,levels=weeks),
                                                        factor(round(temp),levels=-30:30)),mean))))
### also, put the predictions into the original data files
for(i in names(morts)) for(j in names(gamfits[[i]])) morts[[i]][[paste("pred",j,sep="")]]<-fitted(gamfits[[i]][[j]])
dimnames(preds[[1]][[1]])
### contour plots
for(i in names(preds)) {
  png(file=paste("out/",i,"_an2_contour.png",sep=""),width=1400,height=1000)
  par(mfrow=c(5,4),mar=c(2,2,2,2))
  for(j in names(mods)) {
      cat(i,j,"\n")
      plotpredcont(preds,i,j)
  }
  dev.off()
}
sapply(morts,function(a) range(sqrt(a$n/1000)))
par(mfcol=c(2,2),mar=c(3,3,2,1),mgp=2:0)
for(i in names(preds)) plotpredcont(preds,i,"MTn")
for(i in names(preds)) plotpredcont(preds,i,"nTS")
for(i in names(preds)) plotpredcont(preds,i,"MTS")

par(mfrow=c(2,4))
for(i in names(preds)) plotpredpair(i,"MTn")
for(i in names(preds)) plotpredpair(i,"nTS")
for(i in names(preds)) plotpredpair(i,"MTS")
for(i in names(preds)) plotpredpair(i,"MRS")
for(i in names(preds)) plotpredpair(i,"MTX")
for(i in names(preds)) with(morts[[i]],plot(seas,predMTn,main=i))


#filled.contour(weeks,-30:30,preds$fi$Mnn,zlim=c(800,1200))
#filled.contour(weeks,-30:30,preds$fi$nTS,zlim=c(800,1200))
#filled.contour(weeks,-30:30,preds$fi$MTS,zlim=c(800,1200))
#filled.contour(weeks,-30:30,preds$fi$MTR,zlim=c(800,1200))
#filled.contour(weeks,-30:30,preds$fi$MTB,zlim=c(800,1200))
#filled.contour(weeks,-30:30,preds$fi$MTI,zlim=c(800,1200))
#filled.contour(weeks,-30:30,preds$fi$MTX,zlim=c(800,1200))
### perspective plots
for(i in names(preds)) {
    png(file=paste("out/",i,"_an2_persp.png",sep=""),width=1200,height=1000)
    par(mfcol=c(5,4),mar=c(2,2,2,2))
    for(j in names(mods)) {
        if(!is.null(mods[[j]]))
            persp(weeks,-30:30,preds[[i]][[j]],zlim=range(morts[[i]]$n),main=paste(i,j),theta=40,ylab="",xlab="",zlab="")
        else frame() # empty plot
    }
    dev.off()
}

### different graphs using mid month data
print(midmonth<-paste(round(with(preddata,tapply(seas*366,month,min)))))
(monthvec<-with(preddata,tapply(month,seas*366,max)))
mobs<-lapply(obs,function(a) to.na.0(apply(a,2,tapply,monthvec,sum,na.rm=TRUE)))
mobs[[1]]
#(preds[[i]]$STS)[midmonth,]-0*to.na.0(mobs[[i]])

### fits of sorts
png(file=paste("out/eu_an2_interaction.png",sep=""),width=1000,height=500*length(preds))
par(mfcol=c(length(preds),1),mar=c(0,0,1,0)+.5)
for(i in names(preds)) {
  matplot(-30:30,t(preds[[i]]$MTI[midmonth,]),type="l",lty=1,ylim=range(morts[[i]]$n),lwd=1,
          main=i,xaxt="n",xaxs="i",yaxt="n",
          col=colorRampPalette(c("green","red","blue"))(12))
  matlines(-30:30,t(preds[[i]]$MTI[midmonth,]-0*mobs[[i]]),type="l",lty=1,lwd=3,
          col=colorRampPalette(c("green","red","blue"))(12))
  matlines(-30:30,t(preds[[i]]$MTX[midmonth,]-0*mobs[[i]]),type="l",lty=1,lwd=2,
          col=colorRampPalette(c("green","red","blue"))(12))
}
dev.off()

### matrix of predictions
gampreds<-lapply(gamfits,sapply,function(a) predict(a,type="response"))
gamresid<-lapply(gamfits,sapply,function(a) residuals(a))

for(i in names(morts)) {
  png(file=paste("out/",i,"_an2_tempseas.png",sep=""),width=1000,height=1000)
  par(mfrow=c(4,8),mar=c(2,2,2,0.2))
  for(j in names(mods)) {
    for(k in c("temp","seas")) {
      if(!is.null(mods[[j]]))
        matplot(morts[[i]][[k]],gamresid[[i]][,j],pch=16,ylim=range(gamresid[[i]]),main=paste(i,j,k))
    }
  }
  dev.off()
}

### goodness of fit. As we are comparing different countries, we need to rescale:
### Use null model as as max and STB as min
#cols<-c("dk"="red","fi"="lightblue","pt"="purple","ie"="green",uk="blue",de="orange","sp"="yellow")

aics<-(apply(aics0<-sapply(gamfits,sapply,AIC),2,function(a) (a-a["MTB"])/(a["nTn"]-a["MTB"])))
aicdat<-as.data.frame.table(aics)
names(aicdat)<-c("model","country","relAIC")
aicdat$AIC<-as.vector(aics0)
table(aicdat$seasonality<-factor(substring(as.character(aicdat$model),1,1),levels=c("n","M"),labels=c("None","Sin")))
table(aicdat$adjust     <-factor(substring(as.character(aicdat$model),2,2),levels=c("T","R"),labels=c("Raw","Adj")))
table(aicdat$modeltype  <-factor(substring(as.character(aicdat$model),3,3),levels=c("n","S","I","B","X"),
                                 labels=c("None","Spline","Interaction","Spline IA","Surface")))
head(aicdat)
aicdat<-aicdat[order(aicdat$modeltype),]
png(file=paste("out/eu_an2_aic.png",sep=""),width=1400,height=1000,pointsize=24)
trellis.par.set(fontsize=list(text=24,points=16))
trellis.par.set(superpose.symbol=list(cex=1,pch=16,col=1:4))
trellis.par.set(superpose.line  =list(lwd=1,lty= 1,col=1:4))
dotplot(modeltype~AIC|country,groups=interaction(seasonality,adjust),data=aicdat,
        scales=list(relation="free"),type="o",pch=16,auto.key=list(space="right",lines=TRUE))
dev.off()

