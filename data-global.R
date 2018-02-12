################################################################################
###
### data.R
###
### Import data from the NOAA site - ALL EUROPEAN COUNTRIES
###
### Source utility functions
source("fun/auxfun.R")
source("fun/data.R")
require("ISOweek")
### TODO: depend on https://github.com/thl-mjv/euromomo-temp
### Load sites metadata. FALSE means that use the predownloaded and extended sites data
sites<-noaaGetSites(FALSE)

### load the defaults
loadDefaults()
print(getOption("tempmomo")$all)

### DON'T load previously downloaded data, if exists
### Instead use country specific files so that we have OS level assurance of the date of creation
## if(file.exists("data/weather-global.rda")) load("data/weather-global.rda")

### A bit of existentialism: What is Europe?
with(sites,table(Europe,EuropeProper))

### Load data. This will take time, especially first time.
### A version of this data will be made available at 
### "https://owncloud.thl.fi/public.php?service=files&t=a3512ac71e9ab593421171843af081c5"
### No list of all countries
## if(!exists("datas")) datas<-list()
## TODO (20150109) Make sure end of year goes smoothly
print(use.this.year<-tempoption("all","years","this",FALSE)) # DONE: from defaults
print(use.lag      <-tempoption("all","years","lag" ,1    )) # DONE: from defaults getOption("tempmomo")$all$years$lag)  # DONE: from defaults
print(years<-as.numeric(format(Sys.Date(),"%Y"))-8:(!use.this.year))
print(basepath     <-tempoption("all","download","dir","./"))# DONE: from defaults
## Files we already got
files<-list.files(path=paste0(basepath,"/download/"),patt=".rda",full=TRUE)
times<-sapply(files,file.mtime)
names(times)<-gsub(".*/download//","",gsub("[.]rda","",files))
print(times<-(sort((times[match(european.countries,names(times))]-as.numeric(Sys.time()))/(24*3600))))
### All European countries
### TODO: should we just concentrate on those that have not been updated after last Sunday?
###       also, size does matter.
set.seed(Sys.time()) # we don't mind fully random order: this is a hash, not simulation
### Do in the reverse order of need
dothese<-names(times)
#dothese<-"LU"
for(i in dothese) {
    cat(as.character(Sys.time()),":",i,"START --------------------------------------\n")
    savefile<-paste(basepath,"/download/",i,".rda",sep="")
    if(file.exists(savefile)) old<-readRDS(savefile) else old<-NULL
    ii<-tolower(i)
    if(!is.null(old)) {
        cat(as.character(Sys.time()),":",i,"exists!\n")
        last.date<-attr(old,"downloaded")
    } else {
        cat(as.character(Sys.time()),":",i,"does not exist\n")
        last.date<-NULL
    }
    if(is.null(last.date)) last.date<-as.Date("1999-1-1") # long time ago (in a galaxy far far away)     
    ## Use real NUTS area codes based on real maps (but lat/lon coordinates)
    ## Uplag means that we will not update files that are less than week old
    ## ALSO FIXME: SHOULD USE (last day of year<last update) (in noaaGetCountry)
    if(as.numeric(Sys.Date()-last.date)>use.lag) { 
        cat(as.character(Sys.time()),":",i,as.character(last.date),": downloading:\n")
        print(system.time(tmp<-noaaGetCountry(i,force=FALSE,years=years,uplag=use.lag,thisy=use.this.year,countrycode="NUTS0i",basepath=basepath)))
        cur.date<-attr(tmp,"downloaded")
    } else {
        cat(as.character(Sys.time()),":",i,as.character(last.date),": using previous\n")
        tmp<-old
        cur.date<-last.date
    }
    cat("cur",as.character(cur.date),"last",as.character(last.date),"now",as.character(Sys.Date()),"\n")
    ### did anything change?
    if(cur.date>last.date) {
        cat(as.character(Sys.time()),":",i,"Saving new version\n")
        try(wk<-ISOweek(tmp$date))
        try(tmp$year<-as.numeric(substring(wk,1,4)))
        try(tmp$week<-as.numeric(substring(wk,7,9)))
        saveRDS(tmp,file=savefile)
    } else {
        cat(as.character(Sys.time()),":",i,"not saving new version\n")
    }
    ## do we have the stata file? If we do, should we update
    statafile<-paste(basepath,"/out/daily-",ii,".dta",sep="")
    if(cur.date>last.date | !file.exists(statafile)) {
        cat(as.character(Sys.time()),":",i,"Stata\n")
        datas<-list(tmp)
        names(datas)<-ii
        try(writeStatas(datas,"daily",select=ii,basepath=basepath)) # only write the stata files if truly changed!
        ## copy these to output place somewhere determined by defaults-local.txt
        system(paste("cp ",statafile," ",tempoption("all","out","dir","."),sep=""))
    } else {
        cat(as.character(Sys.time()),":",i," Stata already done\n")
    }
    cat(as.character(Sys.time()),":",i,"END ---------------------\n\n")
}
### Moved the rest of the analysis to explore-weather.R
quit("no")
