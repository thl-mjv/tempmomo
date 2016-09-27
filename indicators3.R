###
### third attempt of indicators: just push through original variables
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
lapply(datas,summary)
### default aggregations
funs<-list("All"=namean)
names(datas[[1]])
weeks<-list()
for(i in names(datas)) {
  cat(i,"...")
  tmp<-try(weatherAggre(datas[[i]],
                             stations=stations,
                             vars=list("temp"="temp","dewp"="dewp","stp"="stp","visib"="visib",
                                 "wdsp"="wdsp","prcp"="prcp","sndp"="sndp"),
                             sfun=funs,dfun=funs))
  if(!inherits(tmp,"try-error")) {
    weeks[[i]]<-tmp
    cat("OK\n")
  } else {
    cat("FAIL\n")
  }
}
warnings()
sapply(weeks,names)
lapply(weeks,head)

### mortality
if(file.exists("data/mort-raw.rda")) load("data/mort-raw.rda")
if(!exists("morts")) morts<-list()
for(i in names(weeks)) {
    if(!i%in%names(morts))
        morts[[i]]<-loadMortData(i)
    if("group"%in%names(morts[[i]])) {
        print(dthtable<-rev(sort(with(morts[[i]],tapply(n,group,sum)))))
        morts[[i]]<-subset(morts[[i]],group==(names(dthtable)[1]))
    }
}
sapply(morts,nrow)
sapply(morts,function(a) try(table(a$group)))
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

save(morts,countries,file="data/mort3.rda")
