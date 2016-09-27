###
### first attempt of indicators
### 
### prerequisities
### - modify (or create) defaults-local.txt
###   - check defaults-global.dat for examples
###     - insert xx.momo.path=/path/to/a-momo-outputs as defined in MOMOmaster.do (global WDIR)
###     - OR save your mortality data into download
###       use name mort_xx.dat
###       variables "n" (number of deaths) and "date" (date of the monday of the week)
###       optionally "pop" for population, age groups etc
### - run mort-data.R produces data/mort-raw.rda
### - run data.R  produces data/weather.rda
###
source("fun/auxfun.R")
source("fun/data.R")
source("fun/analysis.R")
source("fun/tempt.R")
### load defaults so that no need for country specific code
loadDefaults()
### local data only, coz we only have local mortality
loadWeatherData(TRUE) 
names(datas)
## use all countries - new 2014-01-09
print(countries<-names(datas) )


pops<-sites$pop0
names(pops)<-sites$usaf
### Use the tools in tempt.R rather than the more crude weatherAggre
### 
### for each country
###  mean temp mean over stations mean over dow - mean temperature
###  mean temp wgt mean over stations mean over dow - NUTSmet temperature
###  min temp min over stations min over dow    - min temperature
###  max temp max over stations max over dow    - min temperature
###  deseasonalized mean temp with deseasonalization
###   stationwise daily data
###   daily mean over stations
###   weekly mean for station
###   weekly mean over stations
###  mean apparent temperature with apparent temp calculated as deseas above
###  number of stations over 2sd deseasonal daily stationwise data
indicators<-list()
for(i in countries) {
  cat(i,"...")
  ii<-list()
  t1<-tempt("temp","date","site",data=datas[[i]]) # basis for many indicators
  t2<-tempt(c("temp","dewp"),"date","site",data=datas[[i]])
  ii$   temp<-naggre(raggre(taggre(t1,fun=namean),fun=namean,agg="NUTS0"),mtemp=temp)
  ii$  wtemp<-naggre(raggre(taggre(t1,fun=namean),fun=namean,agg="NUTS0",weights=pops),wtemp=temp)
  ii$mintemp<-naggre(raggre(taggre(tempt("min","date","site",data=datas[[i]]),fun=namin),fun=namin),mintemp=min)
  ii$maxtemp<-naggre(raggre(taggre(tempt("max","date","site",data=datas[[i]]),fun=namax),fun=namax),maxtemp=max)
  ii$ s1temp<-naggre(raggre(taggre(ttransf(t1,type="sin"),fun=namean),fun=namean),s1temp=temp)
  ii$ s2temp<-naggre(raggre(ttransf(taggre(t1,fun=namean),type="sin"),fun=namean),s2temp=temp)
  ii$ s3temp<-naggre(taggre(ttransf(raggre(t1,fun=namean),type="sin"),fun=namean),s3temp=temp)
  ii$ s4temp<-naggre(ttransf(raggre(taggre(t1,fun=namean),fun=namean),type="sin"),s4temp=temp)
  ii$ at1   <-raggre(taggre(naggre(t2,at1=0.944*temp+0.0153*dewp^2-2.653)))
  ii$ at2   <-raggre(naggre(taggre(t2),at2=0.944*temp+0.0153*dewp^2-2.653))
  ii$ at3   <-taggre(naggre(raggre(t2),at3=0.944*temp+0.0153*dewp^2-2.653))
  ii$ at4   <-naggre(raggre(taggre(t2)),at4=0.944*temp+0.0153*dewp^2-2.653)
  #i13<-???
  out<-as.data.frame(ii[[1]])
  for(j in names(ii)[-1]) out[[j]]<-as.vector(ii[[j]])
  indicators[[i]]<-out
  cat("\n")
}

### save results
save(indicators,countries,file="data/indicators.rda")
#load("data/indicators.rda")
for(i in names(indicators)) {
  png(file=paste("out/",i,"_ind2_%03d.png"),width=1000,height=1000,pointsize=16)
  pairs(indicators[[i]][,c("mtemp","wtemp","mintemp","maxtemp")],pch=16,col=rgb(0,0,0,.25),cex=1,main=i)
  pairs(indicators[[i]][,paste("s",1:4,"temp",sep="")],pch=16,col=rgb(0,0,0,.25),cex=1,main=i)
  pairs(indicators[[i]][,paste("at",1:4,sep="")],pch=16,col=rgb(0,0,0,.25),cex=1,main=i)
  dev.off()
}
sapply(indicators,function(a) with(a,cor(mtemp,cbind(mintemp,maxtemp,wtemp))))
sapply(indicators,function(a) with(a,cor(s1temp,cbind(s2temp,s3temp,s4temp))))
sapply(indicators,function(a) with(a,cor(at1,cbind(at2,at3,at4))))
### conclusion: not much differences
