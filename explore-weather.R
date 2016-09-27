################################################################################
###
### explore-weather.R
###
### Explore the weather data downloaded in data.R
###
### Source utility functions
source("fun/auxfun.R")
source("fun/data.R")

### Load the data saved in data.R
loadWeatherData(FALSE) # use all possible data

basecol<-colorRampPalette(1:8)(length(datas))
### define some colors for further use
(colors<-rgb(t(col2rgb(basecol)),maxColorValue=255,alpha=255)) # dark
(shades<-rgb(t(col2rgb(basecol)),maxColorValue=255,alpha=25 )) # transparent
(shadeS<-rgb(t(col2rgb(basecol)),maxColorValue=255,alpha=75 )) # slightly transparent
names(colors)<-names(shades)<-names(shadeS)<-sort(names(datas))
colors
### how many countries, and how to put them on screen cleanly
(ncountries<-length(datas))
(scountries<-ceiling(sqrt(length(datas))))
(Scountries<-ceiling(ncountries/scountries))

### Organize slightly differently for descriptive analysis. This takes time, about 12*5=60s when tested (times number of countries included, 900s for all)
wvars<-c("temp","dewp","slp","stp","visib","wdsp","mxspd","gust","max","min","prcp","sndp")
names(wvars)<-wvars
system.time(mats<-lapply(datas,function(data) lapply(wvars,function(a) {
    cat(a,"...")
    cat(system.time(res<-with(data,tapply(get(a),list(date,site),mean)))[3],"\n")
    res
})))

names(mats) # which countries
names(mats[[1]]) # what variables for each country

sapply(mats$li,function(a) range(a,na.rm=TRUE))
lapply(mats,sapply,dim)

ranges<-cbind(apply(sapply(mats,sapply,min,na.rm=TRUE),1,min),
              apply(sapply(mats,sapply,max,na.rm=TRUE),1,max))
    
ranges
### Simple visualisation, one graph for each country. 
for(j in names(mats)) {
    cat(j,"... ")
    png(paste("out/",j,"01.png",sep=""),width=1200,height=800,pointsize=10)
    par(mfcol=c(4,3),oma=c(0,0,2,0),mar=c(2,2,2,0))
    for(i in names(mats[[j]])) {
        cat(i,"/ ")
        if(any(!is.na(mats[[j]][[i]]))) {
            matplot(numdate(as.Date(dimnames(mats[[j]][[i]])[[1]])),
                    mats[[j]][[i]],
                    ylim=ranges[i,],yaxs="i",xaxs="i",
                    type="l",lty=1,col=rgb(0,0,0,1/ncol(mats[[j]][[i]])),main=toupper(i),ylab="",xlab="")
        } else {
            frame()
        }
    }
    title(main=toupper(j),outer=TRUE)
    dev.off()
    cat("\n")
}

### explorations
sapply(mats[[1]],dim)
lapply(mats[[1]],head,1)
lapply(mats[[1]],apply,2,mean,na.rm=TRUE)

## what is the correlation between instruments within station?
system.time(stcorr<-lapply(mats,function(mat) sapply(1:ncol(mat[[1]]),function(a) cor(sapply(mat,function(b) b[,a]),use="pair")))) # about 2s
for(i in c("mean","min","max","sd")) {
    png(paste("out/eu02",i,".png",sep=""),width=1200,height=800,pointsize=10) # all ccountries in one graph
    par(mfcol=c(Scountries,scountries),pty="s")
    for(j in names(mats)) {
        
        image(1:12,1:12,matrix(apply(stcorr[[j]],1,get(i),na.rm=TRUE),nr=12),breaks=seq(-1,1,length=100),
              col=colorRampPalette(c("red","white","blue"))(99),yaxt="n",xaxt="n",main=j,ylab="",xlab="")
        for(k in 1:2) axis(k,at=1:12,names(mats[[j]]),las=2)
    }
    dev.off() # close the graph
}

## what is the within day variance as compared to total variance?
variances<-lapply(mats,sapply,function(a) {
    m1<-apply(a,1,mean,na.rm=TRUE) # ave over stations
    d1<-apply((a-m1)^2,1,mean,na.rm=TRUE) # var bw stations
    c(mean(m1,na.rm=TRUE),mean(d1,na.rm=TRUE),var(m1,na.rm=TRUE))
})
t(apply(variances$fi,2,function(a) c(a,a[2]+a[3],100*a[2:3]/sum(a[2:3]))))
(icc<-sapply(variances,apply,2,function(a) 100*a[3]/sum(a[2:3])))

png("out/eu03.png",width=1200,height=800,pointsize=10)
par(mfcol=c(3,ncountries),pty="m",mar=c(0,0,2,0))
for(j in sort(names(mats))) {
    cat(j,"...")
    matplot(mats[[j]][[1]],type="l",lty=1,col=rgb(0,0,0,.1),ylim=c(-40,40),ylab="",main=j,xaxt="n",yaxt="n")
    matplot(apply(mats[[j]][[1]],1,mean,na.rm=TRUE),type="l",lty=1,col=rgb(0,0,0,1),ylim=c(-40,40),ylab="",xaxt="n",yaxt="n")
    matplot(t(apply(mats[[j]][[1]],1,function(a) a-mean(a,na.rm=TRUE))),type="l",lty=1,col=rgb(0,0,0,.1),ylim=c(-40,40),ylab="",xaxt="n",yaxt="n")
}
cat("\n")
dev.off() # close the graph

## what is the coefficient of variation?
lmes<-lapply(mats,lapply,apply,1,mean,na.rm=TRUE)
lses<-lapply(mats,lapply,apply,1,sd  ,na.rm=TRUE)
lsms<-lapply(mats,lapply,apply,1,function(a) sum(!is.na(a)))
(lmem<-apply(sapply(lmes,sapply,min,na.rm=TRUE),1,min))
(lmex<-apply(sapply(lmes,sapply,max,na.rm=TRUE),1,max))
(lsex<-apply(sapply(lses,sapply,max,na.rm=TRUE),1,max))
(devs<-pmax(c(lmex-lmem,lsex)))

png(paste("out/eu04.png",sep=""),width=1200,height=800,pointsize=10)
par(mfcol=c(3,4),pty="s")
for(i in names(mats[[1]])) {
    dev<-c(0,1)*devs[i]
    plot(lmes[[1]][[i]],lses[[1]][[i]],type="n",main=paste(i,sep=""),
         ylab="SD",xlab="Mean",
         xlim=lmem[i]+dev,
         ylim=        dev)
    for(j in names(mats)) {
        points(lmes[[j]][[i]],lses[[j]][[i]],pch=16,col=shades[j],cex=1*lsms[[j]][[i]]/max(lsms[[j]][[i]]))
    }
    legend("top",col=colors,pch=15,cex=.6,legend=paste(toupper(names(colors)),"(",round(icc[i,],1),"%)",sep=""),ncol=5,pt.cex=2)
}
dev.off() # close the graph

### what is the seasonal heteroschedasticity (as modelled using the overall mean?)
means<-lapply(mats ,sapply,apply,1,mean,na.rm=TRUE)
nqpoints<-1
quant<-lapply(means,function(a) try(lapply(apply(a,2,tapply,format(as.Date(dimnames(a)[[1]]),"%m"),quantile,seq(0,1,length=2*nqpoints+1),na.rm=TRUE),sapply,function(a) a)))
str(quant)
(qlims<-apply(do.call("rbind",lapply(quant,sapply,range)),2,range,na.rm=TRUE))
quant[[1]][[1]]
for(sc in 0:1) {
    png(paste("out/eu05_",sc,".png",sep=""),width=1200,height=800,pointsize=10)
    par(mfcol=c(3,4),pty="m",mar=c(2,2,2,1))
    for(i in names(quant$fi)) {
        for(j in 1:length(mats)) {
            qmid<-quant[[j]][[i]][nqpoints+1,]
            if(j==1) plot(1:12,qmid-sc*qmid,main=i,ylab="",xlab="",type="n",ylim=qlims[,i]-sc*mean(qmid))
            lines(1:12,qmid-sc*qmid,col=colors[j])
            for(k in 1:nqpoints) {
                polygon(c(1:12,12:1),c(quant[[j]][[i]][k,]-sc*qmid,rev(quant[[j]][[i]][2*nqpoints+2-k,]-sc*qmid)),border=NA,col=shades[j])
            }
        }
    }
    legend("topright",fill=colors,legend=names(mats),border=NA)   
    dev.off() # close the graph
}
### other reports?

### TODO (7.11.2013): all the above is daily vs weekly data. What is the effect of weekly aggregation?

### OK let's calculate the ranges for temp and mean seasonal range
source("fun/analysis.R")
temp.meanmean<-sapply(mats,function(a) mean(a$temp,na.rm=TRUE))
(temp.meanrange<-sapply(mats,function(a) range(a$temp,na.rm=TRUE)))
(temp.tmaxrange<-sapply(mats,function(a) range(a$max ,na.rm=TRUE)))
temp.tminrange<-sapply(mats,function(a) range(a$min ,na.rm=TRUE))
temp.tmodrange<-sapply(mats,function(a) with(as.data.frame.table(a$temp),range(fitted(lm(Freq~momosin(as.Date(as.character(Var1))))))))
(temp.table<-t(rbind(temp.meanrange,temp.tminrange[1,],temp.tmaxrange[2,],temp.tmodrange))[order(temp.meanmean),])
plot(apply(temp.meanrange,2,diff),apply(temp.tmodrange,2,diff),type="n")
text(apply(temp.meanrange,2,diff),apply(temp.tmodrange,2,diff),dimnames(temp.meanrange)[[2]])
matplot(temp.table,nrow(temp.table):1,type="l",col=c(1,1,2,2,3,3),lty=1,yaxt="n",ylab="")
mtext(side=2,at=nrow(temp.table):1,text=dimnames(temp.table)[[1]],las=2)

