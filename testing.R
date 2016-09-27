################################################################################
###
### Misc tests, not for general use
###

library("RCurl")
getForm("https://owncloud.thl.fi/public.php?service=files&t=2254dda583b2bae59635c1b09d39efc2&download&path=//mort_fi.dat",.params=list(password="Elsinore08"))
getURL("https://owncloud.thl.fi/public.php?service=files&t=a3512ac71e9ab593421171843af081c5&download&path=//mort_fi.dat")
ls(pos="package:RCurl")
?postForm


loadWeatherData(FALSE) 


sapply(datas,function(a) table(with(a,tapply(!is.na(temp),site,mean)))) #almost all have all
names(datas[[1]])
names(sites)
sapply(names(datas),function(a) sample(subset(sites,NUTS0i==toupper(a))$usaf,1))
subset(sites,NUTS0i=="FI" & grepl("HELSINKI",station.name))
subset(sites,NUTS3i=="FI1B1")
with(sites,as.Date(tapply(begin,NUTS0i,min,na.rm=TRUE),"1970-1-1"))
with(sites,as.Date(tapply(begin,NUTS0i,max,na.rm=TRUE),"1970-1-1"))
max(subset(sites,NUTS0i=="FI" )$end,na.rm=TRUE)
subset(sites,NUTS0i=="FI" & na.0(end)==max(na.0(end)))
with(sites,tapply(end,NUTS0i,max))
subset(sites,usaf==29720)

with(subset(datas$fi,site%in%c(29780,29340,29720)),xyplot(temp~date,groups=site,type="l"))

plot(subset(x,STAT_LEVL_==0))
with(subset(sites,usaf%in%samplesites),points(lon/1000,lat/1000,col="red",pch=16))

is.character(c("a","b"))
is.character(Sys.Date())
inherits(Sys.Date(),"Date")
inherits(Sys.time(),"Date")
data.class(Sys.time())

loadDefaults()

foo<-list(fi=loadMortData("fi"),
          pt=loadMortData("pt"),
          dk=loadMortData("dk"))
sapply(foo,names)
summary(foo$fi)

?list.files
tmp<-loadMortData("fi",type="a",year=2014)
head(tmp)
table(tmp$group)
with(foo$fi,plot (n/pop~date,type="l",ylim=c(0,.0003)))
with(foo$pt,lines(n/pop~date,type="l",col="red"))
with(subset(tmp,group=="Total" & date>"2010-1-1"),matplot(numdate(date),cbind(n,nc,fit,excess),type="l",lty=1,lwd=2))
summary(tmp)

loadDefaults()
getOption("tempmomo")$fi$momo$root

library("foreign")

dkmort<-read.dta("download/mort_dk.dta")
names(dkmort)
summary(dkmort)

alldates<-sort(unique(c(as.character(monday(mweeks$pt$date)),as.character(monday(morts$pt$date)))))
table(alldates%in%as.character(monday(mweeks$pt$date)),
      alldates%in%as.character(monday(morts$pt$date)))
plot(as.Date(alldates),alldates%in%as.character(mweeks$pt$date)+2*
     alldates%in%as.character(morts$pt$date),type="l")
cbind(alldates,
      alldates%in%as.character(mweeks$pt$date),
      alldates%in%as.character(morts$pt$date))


plot(fgams$seas ,shade=TRUE,all.terms=FALSE,residuals=TRUE,ylim=c(-8,8),scale=-1)
plot(fgams$tseas,shade=TRUE,all.terms=FALSE,residuals=TRUE,ylim=c(-8,8),scale=-1)

(Effects<-list("gbase"=c(0,0,1),"gtrend"=c(1,0,2),"gseas"=c(0,1,2),"gtseas"=c(1,2,3)))
par(mfcol=c(3,4))
for(i in c("gbase","gtrend","gseas","gtseas")) {
    for(ef in 1:3) {
        if(Effects[[i]][ef]==0) {
            cat("Ei ",i,ef,"\n")
            plot(1,1,main=paste(i,ef),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
        } else {
            cat("Joo",i,ef,"\n")
            plot(fgams[[i]],shade=TRUE,select=Effects[[i]][ef],main=paste(i,ef),ylim=c(-0.3,0.3))
        }
    }
}
        
plot(fgams$gbase ,shade=TRUE,all.terms=FALSE,residuals=FALSE,scale=-1)
plot(fgams$gtrend ,shade=TRUE,all.terms=FALSE,residuals=FALSE,scale=-1)
plot(fgams$gseas ,shade=TRUE,all.terms=FALSE,residuals=FALSE,scale=-1)
plot(fgams$gtseas ,shade=TRUE,all.terms=FALSE,residuals=FALSE,scale=-1,select=4)

par(mfcol=c(3,5),mar=c(3,3,1,1))
plot(fgams$ base ,shade=TRUE,all.terms=FALSE,residuals=FALSE)
plot(fgams$gbase ,shade=TRUE,all.terms=FALSE,residuals=FALSE)

sapply(with(mortdata,tapply(n,co,range)),log)
?smooth.terms
?plot.gam

names(gams[[1]])
?predict.gam
loo<-(predict(gams[[1]][[41]],se=TRUE,type="link"))
dim(foo<-(predict(gams[[1]][[41]],se=TRUE,type="lpmatrix")))
dim(voo<-vcov(gams[[1]][[41]]))
laa<-diag((foo)%*%voo%*%t(foo))
plot(loo$se,sqrt(laa))
abline(a=0,b=1)
str(gams[[1]][[41]])
getS3method("attrassign","lm")
attrassign.gam<-function(object,...) attrassign(model.matrix.lm(object),terms(object))
modeffects(gams[[1]][[41]],morts[[1]])
model.matrix(gams[[1]][[41]])
methods("model.matrix")
xxx
xxx<-1:50*2
yyy<-rpois(length(xxx),exp(2+sin(2*pi*xxx/100)))
plot(xxx,yyy)
ggg<-gam(yyy~s(xxx))
plot(ggg,shade=TRUE)
plot(xxx,yyy,xlim=c(-100,200))
lines(xxx    ,predict(ggg,newdata=data.frame(xxx=xxx)),col="red")
lines(-50:150,predict(ggg,newdata=data.frame(xxx=-50:150)),col="blue")
lines( 50: 75,predict(ggg,newdata=data.frame(xxx=50:75)),col="green")

tmp<-morts[[1]][,c("n","date","temp")]

testgam<-gam(n~s(temp),family=poisson,data=tmp)
plot(predict(testgam,newdata=preddata,newdata.guaranteed=TRUE))
?predict.gam
traceback()

summary(datas$fi)

?term.plot

foo<-gamlss(n~momosin(date,3)+cs(temp),family=PO(),data=morts$fi[,c("n","date","temp")])
term.plot(foo,partial.resid=TRUE)

sapply(morts,names)
mortdata<-do.call("rbind",lapply(names(morts),function(a) cbind(morts[[a]][,c("n","date","temp")],co=a)))
mortdata$con<-factor(mortdata$co)
mortdata$time<-numdate(mortdata$date)
table(mortdata$fi<-1.0*(mortdata$co=="fi"))
table(mortdata$dk<-1.0*(mortdata$co=="dk"))
table(mortdata$pt<-1.0*(mortdata$co=="pt"))
summary(mortdata$seas<-as.numeric(format(mortdata$date,"%j"))/365.25)
summary(mortdata)
library("lattice")
xyplot(n~temp|co,data=mortdata)
fgams<-list()
fgams$ base <-gam(n~0+con                                      +s(temp,by=con,bs="ps"),data=mortdata,family=poisson)
fgams$ trend<-gam(n~0+con+s(time,by=con)                       +s(temp,by=con,bs="ps"),data=mortdata,family=poisson)
fgams$ seas <-gam(n~0+con               +s(seas,by=con,bs="cs")+s(temp,by=con,bs="ps"),data=mortdata,family=poisson)
fgams$ tseas<-gam(n~0+con+s(time,by=con)+s(seas,by=con,bs="cs")+s(temp,by=con,bs="ps"),data=mortdata,family=poisson)
fgams$gbase <-gam(n~0+con                                      +s(temp       ,bs="ps"),data=mortdata,family=poisson)
fgams$gtrend<-gam(n~0+con+s(time       )                       +s(temp       ,bs="ps"),data=mortdata,family=poisson)
fgams$gseas <-gam(n~0+con               +s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson)
fgams$gtseas<-gam(n~0+con+s(time       )+s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson)
fgams$fbase <-gam(n~                                           +s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=fi==1)
fgams$ftrend<-gam(n~      s(time       )                       +s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=fi==1)
fgams$fseas <-gam(n~                    +s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=fi==1)
fgams$ftseas<-gam(n~      s(time       )+s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=fi==1)
fgams$dbase <-gam(n~                                           +s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=dk==1)
fgams$dtrend<-gam(n~      s(time       )                       +s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=dk==1)
fgams$dseas <-gam(n~                    +s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=dk==1)
fgams$dtseas<-gam(n~      s(time       )+s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=dk==1)
fgams$pbase <-gam(n~                                           +s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=pt==1)
fgams$ptrend<-gam(n~      s(time       )                       +s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=pt==1)
fgams$pseas <-gam(n~                    +s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=pt==1)
fgams$ptseas<-gam(n~      s(time       )+s(seas       ,bs="cs")+s(temp       ,bs="ps"),data=mortdata,family=poisson,subset=pt==1)

summary(blöö<-gam(n~0+con+
                  s(time,by=dk,bs="ps")+s(time,by=fi,bs="ps")+s(time,by=pt,bs="ps")+
                  s(seas,by=dk,bs="cs")+s(seas,by=fi,bs="cs")+s(seas,by=pt,bs="cs")+
                  s(temp,by=dk,bs="ps")+s(temp,by=fi,bs="ps")+s(temp,by=pt,bs="ps"),
                  data=mortdata,family=poisson))
par(mfrow=c(3,9))
for(i in 1:3) {
    for(co in 1:3) {
        plot(fgams$ tseas,shade=TRUE,scale=-1,select=co+(i-1)*3,ylim=c(-0.3,0.3))
        plot(blöö        ,shade=TRUE,scale=-1,select=co+(i-1)*3,ylim=c(-0.3,0.3))
        plot(fgams[[paste(c("d","f","p")[co],"tseas",sep="")]],shade=TRUE,scale=-1,select=i,ylim=c(-0.3,0.3))
    }
}
?plot.gam
table(mortdata$con)
summary(blöö)
summary(fgams$ftseas)
summary(fgams$ base)
summary(fgams$gbase)
summary(fgams$ trend)
summary(fgams$gtrend)
summary(fgams$ seas)
summary(fgams$gseas)
summary(fgams$ tseas)
summary(fgams$gtseas)
anova(fgams$base ,fgams$trend,test="Chi") 
anova(fgams$base ,fgams$ seas,test="Chi") 
anova(fgams$trend,fgams$tseas,test="Chi")
anova(fgams$seas ,fgams$tseas,test="Chi")

anova(fgams$base ,fgams$gbase ,test="Chi") 
anova(fgams$trend,fgams$gtrend,test="Chi") 
anova(fgams$seas ,fgams$gseas ,test="Chi") 
anova(fgams$tseas,fgams$gtseas,test="Chi") 

par(mfcol=c(3,8),mar=c(3,3,1,1))
plot(fgams$base ,shade=TRUE,all.terms=FALSE,residuals=TRUE,ylim=c(-8,8),scale=-1)
plot(fgams$trend,shade=TRUE,all.terms=FALSE,residuals=TRUE,ylim=c(-8,8),scale=-1)

xxx<-1:50*2
yyy<-rpois(length(xxx),exp(2+sin(2*pi*xxx/100)))
xxn<- -50:150
yyn<-rpois(length(xxn),exp(2+sin(2*pi*xxn/100)))
plot(c(xxx,xxn),c(yyy,yyn),col=rep(1:2,c(length(xxx),length(xxn))),pch=16)

newdat<-data.frame(xxx=xxn,yyy=yyn,date=as.Date("2010-1-1")+7*xxn)

gg1<-glm(yyy~bs(xxx),family=poisson)
ga1<-gam(yyy~ s(xxx),family=poisson)

cleannames(names(attrassign(gg1)))

with(predict(gg1,se=TRUE               ),matplot (xxx,cbind(fit,fit-2*se.fit,fit+2*se.fit),xlim=c(-10,150),type="l",lty=1,col=1,lwd=3))
with(predict(gg1,se=TRUE,newdata=newdat),matlines(xxn,cbind(fit,fit-2*se.fit,fit+2*se.fit),xlim=c(-10,150),type="l",lty=1,col=2))

with(predict(ga1,se=TRUE               ),matplot (xxx,cbind(fit,fit-2*se.fit,fit+2*se.fit),xlim=c(-10,150),type="l",lty=1,col=1,lwd=3))
with(predict(ga1,se=TRUE,newdata=newdat),matlines(xxn,cbind(fit,fit-2*se.fit,fit+2*se.fit),xlim=c(-10,150),type="l",lty=1,col=2))

mem1<-modeffects(gg1,data=newdat,debug=TRUE,response="yyy",outeffs=c("Intercept","xxx"),clean=TRUE)
mea1<-modeffects(ga1,data=newdat,debug=TRUE,response="yyy",outeffs=c("Intercept","xxx"),clean=TRUE)

names(mem1)
matplot(mem1$effects$add)
matplot(mem1$effects$res)
matplot(mea1$effects$add)
matplot(mea1$effects$res)
(mem1$seaseffs[,"xxx","season",,])

modmat(glm   (yyy~bs(xxx),family=poisson),data=data.frame(xxx=xxn))
modmat(gam   (yyy~ s(xxx),family=poisson),data=data.frame(xxx=xxn))
modmat(gamlss(yyy~pb(xxx),family=PO),data=data.frame(xxx=xxn))

attrassign(glm   (yyy~bs(xxx),family=poisson))
attrassign(gam   (yyy~ s(xxx),family=poisson))
gam   (yyy~ s(xxx),family=poisson)$pterms
traceback()
?model.matrix.gamlss
methods("attrassign")
str(model.matrix(gam   (yyy~s(xxx),family=poisson)))
str((gam   (yyy~s(xxx),family=poisson)))
getS3method("attrassign","default")
getS3method("attrassign","coxph")
getS3method("model.matrix","gam")
model.matrix.gam
predict.gam
predict(ga1,type="lpmatrix",newdata=newdat)

gamme<-list()
for(i in names(morts)) gamme[[i]]<-lapply(gams[[i]],modeffects,data=morts[[i]],clean=TRUE)

lapply(lapply(gamme,sapply,function(a) a$seaseffs["add","TEMP","winter","2009",]),t)

modeffects(gams[[1]][[41]],data=morts[[1]],clean=TRUE,debug=TRUE)
lapply(gams,sapply,function(a) cleannames(names(attrassign(a))))

names(sites)


load("data/weather.rda")
load("data/mort.rda")
length(with(morts[[1]],list(temp)))
source("fun/auxfun.R")
source("fun/analysis.R")
source("fun/data.R")
source("fun/tempt.R")
pops<-sites$pop0
names(pops)<-sites$usaf
t1<-with(datas[[1]],tempt(temp,date,site))
t2<-with(datas[[1]],tempt(dewp,date,site))
t3<-tempt("temp","date","site",data=datas[[1]])
dim(t1a<-taggre(t1))
dim(t2a<-taggre(t2))
summary(as.data.frame(t1b <-raggre(t1a,agg="NUTS3",weights=pops,debug=TRUE)))
summary(as.data.frame(t1c1<-raggre(t1a,agg="NUTS0")))
dim(t2c1<-raggre(t2a,agg="NUTS0"))
dim(t1c2<-raggre(t1b,agg="NUTS0"))
dimnames(t1c1)
t12<-combine(t1c1,t2c1)
summary(as.data.frame(t1))
summary(as.data.frame(t1a))
summary(as.data.frame(t1b))
summary(as.data.frame(t1c1))
summary(as.data.frame(t12))
str(as.data.frame.table(t12[,,1,drop=FALSE]))

as.character(with(subset(sites,usaf%in%tregs(t1a)),NUTS0b[match(tregs(t1a),usaf)]))

library(mgcv)
ttt<-raggre(taggre(tempt(c("temp","dewp"),"date","site",data=datas[["fi"]])),agg="NUTS0")
tts<-raggre(taggre(ttransf(tempt(c("temp","dewp"),"date","site",data=datas[["fi"]]),type="spl")),agg="NUTS0")
tss<-raggre(ttransf(taggre(tempt(c("temp","dewp"),"date","site",data=datas[["fi"]])),type="sp"),agg="NUTS0")
sss<-ttransf(raggre(taggre(tempt(c("temp","dewp"),"date","site",data=datas[["fi"]])),agg="NUTS0"),type="sp")
tt2<-tempt(as.data.frame(ttt))

matplot(numdate(tdate(ttt)),ttt[,1,],type="l",col=1,ylim=c(-5,5))
matlines(numdate(tdate(tts)),tts[,1,],col=2)
matlines(numdate(tdate(tss)),tss[,1,],col=3)
matlines(numdate(tdate(sss)),sss[,1,],col=4)
matplot(tts[,1,],sss[,1,])
pairs(cbind(ttt[,1,1],tts[,1,1],tss[,1,1],sss[,1,1]))
apply(t1,2:3,deseas,as.Date(ttime(t1)))

seasday(ttime(t1))
seasday(t(outer(as.Date(c("2014-1-1","2008-12-31")),-5:5,"+")))
deseas(ttt,(ttime(ttt)))
seasday(ttime(ttt))

numdate(as.Date(ttime(ttt)))
all.equal(ttt,tt2)
str(ttt)
str(tt2)

t4<-tempt(as.data.frame(t1))

all.equal(t1,t4)
tname(t4)
names(attributes(t2))
attr(t4,"regtype")

tttt<-raggre(taggre(tempt(c("temp","dewp"),"date","site",data=datas[["fi"]])))

naggre(tttt,at=temp+dewp)

test<-function(data,...) {
  mc<-match.call(expand.dots=TRUE)
  out<-list()
  for(i in names(mc)[-(1:2)]) {
    cat(i,"\n")
    if(!is.null(i)) {
      out[[i]]<-eval(mc[[i]],env=data)
      print(str(mc[[i]]))
    }
  }
  out
}
(test(head(datas[["fi"]]),aa=temp+dewp))


load("data/indicators.rda")

pairs(indicators$ie[c("at1","at2","at3","at4")])

####
install.packages("ISOweek")
library(ISOweek)
help(package="ISOweek")
ISOweek(Sys.Date())

head(sites)

par(mfcol=c(1,1))
with(subset(sites,NUTS0b=="DE"),plot(lon,lat))


################################################################################
### matthias definition of normal weather
tmp.temp<-sapply(morts,function(a) with(a,predict(glm(temp~momosin(date,4)))))
sapply(tmp.temp,range)
effs<-list()
for(i in names(morts)) {
    effs[[i]]<-with(morts[[i]],{
        mod<-glm(temp~momosin(date,4))
        fit<-fitted(mod)
        cbind(pmax(0,temp-max(fit)),-pmax(0,min(fit)-temp))
    })
}
for(i in names(effs)) with(morts[[i]],matplot(effs[[i]],(n-pred)/sqrt(pred),main=i))
effmods<-list()
temmods<-list()
for(i in names(effs)) {
    effmods[[i]]<-gam(n~s(time,k=1)+s(seas,bs="cc",k=14)+s(effs[[i]][,1],k=4)+s(effs[[i]][,2],k=4),family=poisson,data=morts[[i]])
    temmods[[i]]<-gam(n~s(time,k=1)+s(seas,bs="cc",k=14)+s(temp,k=5),family=poisson,data=morts[[i]])
}
par(mfrow=c(4,7))
for(i in names(morts)) {
    plot(effmods[[i]],shade=TRUE,scale=0)
    plot(temmods[[i]],shade=TRUE,scale=0)
}
?plot.gam
sapply(effmods,plot,shade=TRUE)
par(mfrow=c(4,2))
sapply(temmods,plot,shade=TRUE)
dotplot(cbind(sapply(effmods,AIC),
              sapply(temmods,AIC)))


library(xlsx)
ge<-read.xlsx("~/ownCloud/FluMoMo.xlsx",sheetIndex=2)
ge$date<-with(ge,ISOweek2date(sprintf("%04d-W%02d-1",year,week)))
ge$Year<-as.numeric(format(ge$date,"%G"))
ge$Week<-as.numeric(format(ge$date,"%V"))
with(ge,table(Year==year,Week==week))
summary(ge)
subset(ge,Year!=year | Week!=week)
subset(ge,(week%in%c(1,53) | Week%in%c(1,53))& agegroups=="65+")
       
geT<-as.data.frame.table(with(ge,addmargins(tapply(deaths,list(Year,Week,agegroups),sum),3)))
names(geT)<-c("year","week","group","n")
geT$date<-with(geT,ISOweek2date(sprintf("%04d-W%02d-1",as.numeric(as.character(year)),as.numeric(as.character(week)))))
summary(geT)
write.table(geT,file="download/mort_de.txt",sep=",",row.names=FALSE)

subset(ge ,agegroups=="65+" & abs(as.numeric(date-as.Date("2005-12-31")))<14)
subset(geT,   group =="65+" & abs(as.numeric(date-as.Date("2005-12-31")))<14)


tmp<-loadMortData("de",type="txt",extension="txt",datetype="week")
tmp2<-loadMortData("de")
loadDefaults()



with(tmp2,xyplot(n~date|group,type="l",layout=c(1,5)))
with(ge ,xyplot(deaths~week|agegroups,groups=year,type="l"))
with(geT,xyplot(n     ~week|   group ,groups=year,type="l"))
par(mfcol=c(1,1))
matplot(with(subset(geT,group=="Sum"),tapply(n,list(week,year),sum)),type="l")
?ISOweek2date
ls(pos="package:ISOweek")
ISOweek2date("2014-W10-1")

###
### fits?
for(i in names(morts)) morts[[i]]$tempbase<-fitted(glm(temp~momosin(date,4),data=morts[[i]]))
xscl<-function(a) (a-mean(a))/sd(a)
par(mfcol=mkpar(length(morts)))
for(i in names(morts)) with(morts[[i]],matplot(time,cbind(-xscl(tempbase),xscl(pred)),main=i,type="l"))
par(mfcol=mkpar(length(morts)))
for(i in names(morts)) with(morts[[i]],plot(tempbase,pred,main=i,type="l"))
par(mfcol=mkpar(length(morts)))
for(i in names(morts)) with(morts[[i]],plot(temp-tempbase,n-pred,main=i,type="l"))
par(mfcol=mkpar(length(morts)))
for(i in names(morts)) with(morts[[i]],matplot(seas,cbind(-xscl(tempbase),xscl(pred)),main=i,type="p"))
par(mfcol=mkpar(length(morts)))
for(i in names(morts)) with(morts[[i]],matplot(seas,cbind(temp-tempbase,(n-pred)/sqrt(pred)),main=i,type="p"))
par(mfcol=mkpar(length(morts)))
for(i in names(morts)) with(morts[[i]],plot(temp.adj2,(n-pred)/sqrt(pred),main=i,type="p",
                                            pch=c("J","F","M","A","M","J","J","A","S","O","N","D")[as.numeric(monthf)]))
with(morts$fi,table(ceiling(seas*4)))

################################################################################
### export
###
library("foreign")
foreign:::writeForeignStata
source("fun/data.R")
writeStatas(weeks,"wks")

head(datas$be)
names(sites)

table(outfail<-as.numeric(substring(system("grep -c 'no room' daily*.log",TRUE),14)))
table(outobs <-as.numeric(substring(system("grep -h 'observations' daily*.log|sed 's/observations read)//'",TRUE),2)))
(outlines<-as.numeric(substring(system("wc -l out/daily*.dat",TRUE),3,9)))
tapply(outlines[1:length(outfail)],outfail,range)
outobs/(outlines[1:length(outobs)])

################################################################################
###
### file info
str(file.info("testing.R"))

21353/3600

 if(file.exists("data/weather-global.rda")) load("data/weather-global.rda")
names(datas)
for(i in names(datas)) saveRDS(datas[[i]],file=paste("download/",toupper(i),".rda",sep=""))

################################################################################
###
dev.list()
summary(sites)
sort(table(sites$country))
summary(sites$lon)
