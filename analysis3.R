###
### third attempt of analysis
### compare different indicators
### 
source("fun/auxfun.R")
source("fun/data.R")
source("fun/analysis.R")
library("mgcv")
### load defaults so that no need for country specific code
loadDefaults()

### load the combined weather and mortality data
load("data/mort.rda")
for(i in names(morts)) morts[[i]]$seas  <-seasday(morts[[i]]$date)
for(i in names(morts)) morts[[i]]$time  <-numdate(morts[[i]]$date)
for(i in names(morts)) morts[[i]]$month <-as.numeric(format(morts[[i]]$date,"%m"))
for(i in names(morts)) morts[[i]]$monthf<-factor(morts[[i]]$month)
for(i in names(morts)) print(table(morts[[i]]$seasf <-factor(morts[[i]]$month%in%4:9,levels=c(TRUE,FALSE),labels=c("Hot","Cold"))))
for(i in names(morts)) morts[[i]]$pop   <-rep(max(subset(sites,NUTS0b==toupper(i))$pop0),nrow(morts[[i]]))
for(i in names(morts)) morts[[i]]$country<-rep(i,nrow(morts[[i]]))

### above this same as analysis2.R
################################################################################
###
### Let's try different indicators
###
### Fit the models
gam3base<-list()
gam3fits<-list()
for(i in names(morts)) {
    gam3fits[[i]]<-list()
    gam3base[[i]]<-gam(n~s(time,k=3)+momosin(date,2),data=morts[[i]],family=poisson)
    for(l in c("",".adj1",".ind")) {
        for(j in c("","w")) {
            for(k in c("temp","at","wci","humindex")) {
                for(b in c("",",by=seasf")) {
                    cat(i,j,k,l,b)
                    lab<-paste(j,k,l,sep="")
                    nam<-paste(lab,b,sep="")
                    if(l==".ind")
                        mod<-paste("~.+s(",lab,"1",b,")+s(",lab,"2",b,")",sep="")
                    else
                        mod<-paste("~.+s(",lab,b,")",sep="")
                    cat(mod)
                    tmp<-try(update(gam3base[[i]],formula(mod)))
                    if(!inherits(tmp,"try-error")) {
                        gam3fits[[i]][[nam]]<-tmp
                        cat(" OK\n")
                    } else {
                        cat(" FAIL\n")          
                    }
                }
            }
        }
    }
}
###
aictab3<-apply(aictab3a<-sapply(gam3fits,sapply,AIC),2,function(a) (a-min(a))/diff(range(a)))
aicdat3<-as.data.frame.table(aictab3)
names(aicdat3)<-c("lab","country","relAIC")
aicdat3$AIC<-as.vector(aictab3a)
table(aicdat3$weight<-ifelse(grepl("^w[tawh]",as.character(aicdat3$lab)),"weight","nonweight"))
table(aicdat3$adj   <-ifelse(grepl(".adj",as.character(aicdat3$lab)),"adj",
                             ifelse(grepl("[.]ind",as.character(aicdat3$lab)),"ind",
                                    "unadj")))
table(aicdat3$by    <-ifelse(grepl("by=",as.character(aicdat3$lab)),"byseas","noby"))
table(aicdat3$windic <-gsub(".adj1","",gsub("[.]ind","",gsub(",by=seasf","",as.character(aicdat3$lab)))))
table(aicdat3$indic <-with(aicdat3,factor(toupper(substring(windic,1+(weight=="weight"))),levels=c("TEMP","AT","WCI","HUMINDEX"))))
head(aicdat3)
library(lattice)
(trellis.par.get()$fontsize)
png(file="out/eu_an3_aic1.png",width=1000,height=1400,pointsize=48)
trellis.par.set(fontsize=list(text=24,points=16))
trellis.par.set(superpose.line=list(col=c(1,1,2,2,3,3),lty=c(1,2,1,2,1,2)),
                superpose.symbol=list(col=c(1,1,2,2,3,3),pch=rep(c(1,16),c(6,6))))
dotplot(indic~AIC|country,data=aicdat3,groups=interaction(weight,adj,by),type="o",layout=c(1,length(morts)),
        auto.key=list(space="right",lines=TRUE),scales=list(relation="free"))
dev.off()
png(file="out/eu_an3_aic2.png",width=1400,height=1000,pointsize=24)
trellis.par.set(fontsize=list(text=24,points=16))
dotplot(indic~relAIC|adj*weight*by,data=aicdat3,groups=country,type="o",
        layout=c(3,4),
        auto.key=list(space="right",lines=TRUE))
dev.off()
trellis.par.set(fontsize=list(text=12,points=8))

#par(mfrow=c(4*length(morts),4),mar=c(1,1,2,1))
#for(i in names(gam3fits)) for(j in names(gam3fits[[i]])) plot(gam3fits[[i]][[j]],select=2,shade=TRUE,main=paste(i,j))

for(i in names(gam3fits)) morts[[i]][[paste("pred",sep="")]]<-fitted(gam3base[[i]])
for(i in names(gam3fits)) for(j in names(gam3fits[[i]])) morts[[i]][[paste("pred",j,sep="")]]<-fitted(gam3fits[[i]][[j]])

###
### all of a sudden the code moves to visualisations for the Stockholm meeting presentation
### 
get.vars<-function(va) unlist(sapply(morts,function(a) eval(parse(text=va),env=a)))
plot.vars<-function(x,type=c("x","y"),legpos="bottom",postfix=FALSE,weighted=FALSE,adjusted=FALSE,xlab=NULL,
                    ylab=NULL,...) {
    type<-match.arg(type)
    ccc<-rep(names(morts),sapply(morts,nrow))
    cols<-c("dk"="red","fi"="blue","pt"="purple","ie"="green")
    par(mfcol=c(2,2))
    for(i in c("temp","at","wci","humindex")) {
        if(postfix) xx<-gsub("%",i,x) else xx<-x
        if(weighted) ii<-paste("w",i,sep="") else ii<-i
        if(adjusted) ii<-paste(ii,".adj1",sep="")
        if(is.null(ylab)) ylab<-ifelse(type=="x",ii,xx)
        if(is.null(xlab)) xlab<-ifelse(type=="y",ii,xx)
        if(type=="x") plot(get.vars(xx),get.vars(ii),col=cols[ccc],pch=16,main=i,xlab=xlab,ylab=ylab,...)
        if(type=="y") plot(get.vars(ii),get.vars(xx),col=cols[ccc],pch=16,main=i,ylab=ylab,xlab=xlab,...)
        legend(legpos,pch=16,col=cols,legend=names(cols))
    }
}
names(gam3fits[[1]])
### TODO: publish these so they can be included in the presentation
png(file="out/eu_indicators_w_nw.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,3,1),mgp=2:0)
plot.vars("w%","y","topleft",weighted=FALSE,postfix=TRUE)
dev.off()
png(file="out/eu_indicators_a_na.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,3,1),mgp=2:0)
plot.vars("%.adj1","y","topleft",weighted=FALSE,postfix=TRUE)
dev.off()
png(file="out/eu_indicators_vs_temp.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,3,1),mgp=2:0)
plot.vars("temp","x","topleft",weighted=FALSE,ylim=c(-30,40),xlim=c(-30,40))
dev.off()
png(file="out/eu_indicators_vs_seas.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,3,1),mgp=2:0)
plot.vars("seas","x","bottom",weighted=FALSE)
dev.off()
png(file="out/eu_indicators_adj_vs_seas.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,3,1),mgp=2:0)
plot.vars("seas","x","topleft",adjusted=TRUE)
dev.off()
png(file="out/eu_indicators_vs_mort.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,3,1),mgp=2:0)
plot.vars("10000*n/pop","y","bottomleft")
dev.off()
png(file="out/eu_seas_vs_mort.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,1,1),mgp=2:0,mfcol=c(1,1))
cols<-c("dk"="red","fi"="blue","pt"="purple","ie"="green")
plot(get.vars("seas"),get.vars("10000*n/pop"),col=cols[get.vars("country")],pch=16,
     ylab="Mortality",xlab="seas",xaxs="i")
legend("top",pch=16,col=cols,legend=names(cols))
dev.off()
png(file="out/eu_indicators_adj_vs_mort_adj.png",width=1400,height=1000,pointsize=24)
par(mar=c(3,3,3,1),mgp=2:0)
plot.vars("(n-pred)/sqrt(pop*pred)","y","topleft",adjusted=TRUE,postfix=FALSE,ylab="deseas. mortality")
dev.off()
plot.vars("(n-pred%)/sqrt(pop*pred%)","y","topleft",postfix=TRUE,weighted=TRUE)
plot.vars("(n-pred%)/sqrt(pop*pred%)","y","topleft",postfix=TRUE,adjusted=TRUE)
plot.vars("(n-pred%.adj1)/sqrt(pop*pred%.adj1)","y","topleft",postfix=TRUE)
plot(get.vars("seas"),get.vars("(n-predtemp)/sqrt(pop*predtemp)"),col=cols[get.vars("country")],pch=16)
plot(get.vars("seas"),get.vars("(n-pred    )/sqrt(pop*pred    )"),col=cols[get.vars("country")],pch=16)
plot.vars("10000*pred/pop","x","topleft",adjusted=TRUE)
plot.vars("momolag(10000*pred/pop,-1)","y","topleft",)
par(mfcol=c(3,3))
for(i in -4:4) plot(get.vars(paste("momolag(temp,",i,")")),
                    get.vars("(n-predtemp)/sqrt(pop*predtemp)"),col=cols[get.vars("country")],pch=16,main=i)
for(i in 0:8) with(morts$fi,matplot(momolag(temp,i),cbind((n-predtemp)/sqrt(pop*predtemp),(n-pred)/sqrt(pop*pred)),type="h",lty=1,
                                    pch=16,ylim=c(-0.003,0.003)))
for(i in 0:8) with(morts$fi,matplot(momolag(temp.adj1,i),cbind((n-predtemp)/sqrt(pop*predtemp),(n-pred)/sqrt(pop*pred)),type="h",lty=1,
                                    pch=16,ylim=c(-0.003,0.003)))

xcors<-lapply(morts,function(a) with(a,{
  yy<-cbind(temp,temp.adj1,n,desn=(n-pred)/sqrt(pred),detn=(n-predtemp)/sqrt(predtemp),detan=(n-predtemp.adj1)/sqrt(predtemp.adj1))
  cbind(cor(momolag(temp     ,0:150,0),yy),
        cor(momolag(temp.adj1,0:150,0),yy),
         2/sqrt(sum(!is.na(temp) & !is.na(n) & !is.na(pred))),
        -2/sqrt(sum(!is.na(temp) & !is.na(n) & !is.na(pred))))
}))
par(mfcol=c(2,5),mar=c(3,3,1,1),mgp=2:0)
for(i in names(xcors)) {
    for(j in c(0,6)){
        col<-c("black","red","blue","purple","yellow","cyan")
        xx<-xcors[[i]][,1:6+j]
        lim<-xcors[[i]][,13:14]
        xmax<-20
        matplot(0:150,xx,main=paste(i,ifelse(j==0,"raw","adjusted"),"temperature"),
                type="l",pch=16,ylim=c(-1,1),lty=rep(1:2,c(6,6)),lwd=2,
                col=col,ylab="Correlation",xlab="Lag, weeks",xlim=c(0,xmax+5),
                xaxs="i",yaxs="i")
        polygon(c(0:150,150:0),c(lim[,1],rev(lim[,2])),col=rgb(0,0,0,.1),border="lightgrey")
        text(rep(xmax,5),xx[xmax+1,],dimnames(xx)[[2]],col=col,adj=0)
        abline(h=-1:1)
  }
}


cols<-c("de"="orange","dk"="red","fi"="blue","pt"="purple","ie"="green")
labs<-c("temp","deseas.temp","deaths","deseas.deaths","temp adjusted deaths","deseas.temp adj deaths")
png("out/eu_autocorr.png",width=1000,height=1200)
par(mfrow=c(6,2),oma=c(0,0,2,0),mar=c(3,3,2,1),mgp=2:0)
for(i in 1:6)for(j in c(0,6))  {
    matplot(0:150,sapply(xcors,function(a) a[,i+j]),type="l",#ylim=c(-1,1),
            xlim=c(0,20),xaxs="i",yaxs="i",
            ylab="Autocorrelation",xlab="Lag",main=paste(labs[i],"vs",ifelse(j==0,"temp","deseas.temp")),
            lty=1,col=cols[names(xcors)])
    sapply(xcors,function(a) polygon(c(0:150,150:0),c(a[,13],rev(a[,14])),col=rgb(0,0,0,.01),border=NA))
    if(i==1 & j==0) legend("bottomleft",legend=names(xcors),lty=1,col=cols[names(xcors)])
}
dev.off()



par(mfcol=c(1,1))
with(morts$fi,matplot(time,cbind(n,pred,predtemp,predtemp.adj1),type="l",lty=1,lwd=c(1,2,2,2)))
with(morts$fi,matplot(time,cbind(desn=(n-pred)/sqrt(pred),detn=(n-predtemp)/sqrt(predtemp),detan=(n-predtemp.adj1)/sqrt(predtemp.adj1)),
                      type="l",lty=1))
with(morts$fi,plot(acf(cbind(temp,temp.adj1,n,desn=(n-pred)/sqrt(pred),detn=(n-predtemp)/sqrt(predtemp),detan=(n-predtemp.adj1)/sqrt(predtemp.adj1)))))

