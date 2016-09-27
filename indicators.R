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
### load defaults so that no need for country specific code
loadDefaults()
### local data only, coz we usually only have the local mortality
loadWeatherData(TRUE) 
names(datas) ## check what data we haz
## use all countries - new 2014-01-09
print(countries<-names(datas) )

sapply(datas,nrow)
##sapply(datas,names)
### Combine two types of aggregation: 
### a) Mean over all stations
### b) Select one site per country. I (MJV) do not know which station would be best, so 
###    let's use the oldest per country (very bad choice)
### this fails with DE for some reason
###(samplesites<-sapply(names(datas),function(a)
###                     with(subset(sites,NUTS0i==toupper(a) & begin<"2008-01-01"),
###                          usaf[begin==min(begin,na.rm=TRUE)][1])))
(samplesites<-sapply(datas,function(a) with(a,names(rev(sort(tapply(!is.na(temp),site,sum)))[1]))))
##summary(subset(datas$de,site==samplesites["de"]))
### ... But if we happen to work with Finnish data, select Kaisaniemi station
###     as it is closest to the center of Helsinki and has no missing values
if("fi"%in%names(samplesites)) samplesites["fi"]<-29780 # TODO: use default site
### default aggregations
funs<-list("All"=namean, # default
            "mintemp"=namin, "maxtemp"=namax,  # aggregated
           "smintemp"=namin,"smaxtemp"=namax) # one station
weeks<-list()
for(i in names(datas)) {
  cat(i,"...")
  stations<-list(stemp=samplesites[i],smintemp=samplesites[i],smaxtemp=samplesites[i])
  tmp<-try(weatherAggre(datas[[i]],
                             stations=stations,
                             vars=list(
                                  "at"="(0.944*temp+0.0153*dewp^2-2.653)",
                                 "wat"="(0.944*temp+0.0153*dewp^2-2.653)",
                                  "wci"="(13.12 + 0.6215*temp-11.37*wdsp^0.16 + 0.3965*temp*wdsp^0.16)",
                                 "wwci"="(13.12 + 0.6215*temp-11.37*wdsp^0.16 + 0.3965*temp*wdsp^0.16)",
                                 "humindex"="(temp+0.55555*(6.11*exp(5417.7530*(1/273.16-1/(dewp+273.16)))))",
                                 "whumindex"="(temp+0.55555*(6.11*exp(5417.7530*(1/273.16-1/(dewp+273.16)))))",
                                  "temp"="temp", "mintemp"="min", "maxtemp"="max",
                                 "stemp"="temp","smintemp"="min","smaxtemp"="max",
                                 "wtemp"="temp"),
                             sfun=funs,dfun=funs,weights=list(
                                                     wtemp="na.0(pop3)",
                                                       wat="na.0(pop3)",
                                                      wwci="na.0(pop3)",
                                                      whumindex="na.0(pop3)"
                                                     )))
  if(!inherits(tmp,"try-error")) {
    weeks[[i]]<-tmp
    cat("OK\n")
  } else {
    cat("FAIL\n")
  }
}
warnings()
sapply(weeks,names)
sapply(weeks,nrow)
lapply(weeks,head)

nn<-length(weeks)
par(mfcol=mkpar(nn))
for(i in names(weeks)) try(with(weeks[[i]], matplot(temp,cbind(wtemp,stemp,at,wat,wci,wwci),main=i)))
par(mfcol=mkpar(nn))
for(i in names(weeks)) try(with(weeks[[i]], matplot(cbind(temp,at,wci),cbind(wtemp,wat,wwci),main=i)))

### mortality
if(file.exists("data/mort-raw.rda"))
    load("data/mort-raw.rda")
if(!exists("morts")) morts<-list()
for(i in names(weeks)) {
    cat(i,"...")
    if(!i%in%names(morts)) {
        cat("load")
        morts[[i]]<-loadMortData(i)
    }
    if("group"%in%names(morts[[i]])) {
        cat("remove group")
        print(dthtable<-rev(sort(with(morts[[i]],tapply(n,group,sum,na.rm=TRUE)))))
        morts[[i]]<-subset(morts[[i]],group==(names(dthtable)[1]))
    }
    cat("\n")
}
sapply(morts,names)
sapply(morts,nrow)
sapply(morts,function(a) try(table(a$group)))
sapply(morts,function(a) try(with(a,tapply(n,group,sum,na.rm=TRUE))))
### combine
for(i in names(weeks)) {
    cat(i,"...")
    for(j in names(weeks[[i]])[-1]) { # assume first variable is "date"
        morts[[i]][[j]]<-weeks[[i]][[j]][match(monday(morts[[i]]$date),monday(weeks[[i]]$date))]        
    }
    cat(nrow(morts[[i]]),"...")
    morts[[i]]<-subset(morts[[i]],!is.na(n) & !is.na(temp) ) #& !is.na(temp1)) allow missing in temp1
    cat(nrow(morts[[i]]),"\n")
}
lapply(morts,summary)

### adjusted
###
for(i in names(morts)) { # all datasets
    for(j in c("temp","at","wci","wtemp","wat","wwci","humindex","whumindex")) { # all variables
        cat(i,j,"...\n")
        fit<-lm(formula(paste(j,"~momosin(date,2)")),data=morts[[i]])
        pred<-predict(fit,newdata=morts[[i]],se=TRUE)
        morts[[i]][[paste(j,".adj1",sep="")]]<-with(pred,morts[[i]][[j]]-fit)
        morts[[i]][[paste(j,".adj2",sep="")]]<-with(pred,(morts[[i]][[j]]-fit)/se.fit)
        morts[[i]][[paste(j,".ind1",sep="")]]<-i1<-pmax(0,morts[[i]][[j]]-max(pred$fit))
        morts[[i]][[paste(j,".ind2",sep="")]]<-i2<-pmin(0,morts[[i]][[j]]-min(pred$fit))
        morts[[i]][[paste(j,".indX",sep="")]]<-i1+i2
    }
}

png("out/eu_an_adj.png",,width=1200,height=500*length(morts),pointsize=20)
par(mfrow=c(length(morts),1))
for(i in names(morts)) {
    ##png(paste("out/",i,"_ind_adj.png",sep=""),width=1200,height=800,pointsize=10)       
    with(morts [[i]],matplot(numdate(date),cbind(temp,temp.adj1,temp.adj2,temp.ind1,temp.ind2,temp.indX),type="l",lty=1,main=i))
    ##dev.off()
}
dev.off()



rbind(sapply(morts,function(a) with(a,cor(temp    ,cbind(wtemp,at,wci,humindex)))),
      sapply(morts,function(a) with(a,cor(at      ,wat))),
      sapply(morts,function(a) with(a,cor(wci     ,wwci))),
      sapply(morts,function(a) with(a,cor(humindex,whumindex))))

for(i in names(morts)) {
    png(paste("out/",i,"_indi.png",sep=""),width=1000,height=1000,pointsize=16)
    with(morts[[i]],pairs(cbind(temp,wtemp,at,wat,wci,wwci,humindex,whumindex)))
    dev.off()
}

save(morts,countries,sites,file="data/mort.rda")
