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
use.this.year<-TRUE # TODO: from defaults
use.lag<-1          # TODO: from defaults
(years<-as.numeric(format(Sys.Date(),"%Y"))-8:(!use.this.year))
## Files we already got
files<-list.files(path="download/",patt=".rda",full=TRUE)
times<-sapply(files,file.mtime)
names(times)<-gsub("download//","",gsub("[.]rda","",files))
(times<-(sort((times[match(european.countries,names(times))]-as.numeric(Sys.time()))/(24*3600))))
### All European countries
### TODO: should we just concentrate on those that have not been updated after last Sunday?
###       also, size does matter.
set.seed(Sys.time()) # we don't mind fully random order: this is a hash, not simulation
### Do in the reverse order of need
for(i in names(times)) {
    cat(as.character(Sys.time()),":",i,"START --------------------------------------\n")
    savefile<-paste("download/",i,".rda",sep="")
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
        print(system.time(tmp<-noaaGetCountry(i,force=FALSE,years=years,uplag=use.lag,thisy=use.this.year,countrycode="NUTS0i")))
        cur.date<-attr(tmp,"downloaded")
    } else {
        cat(as.character(Sys.time()),":",i,as.character(last.date),": using previous\n")
        tmp<-old
        cur.date<-last.date
    }
    ### did anything change?
    if(cur.date>last.date) {
        cat(as.character(Sys.time()),":",i,"Saving new version\n")
        try(wk<-ISOweek(tmp$date))
        try(tmp$year<-as.numeric(substring(wk,1,4)))
        try(tmp$week<-as.numeric(substring(wk,7,9)))
        saveRDS(tmp,file=savefile)
    }
    ## do we have the stata file? If we do, should we update
    statafile<-paste("out/daily-",ii,".dta",sep="")
    if(cur.date>last.date | !file.exists(statafile)) {
        cat(as.character(Sys.time()),":",i,"Stata\n")
        datas<-list(tmp)
        names(datas)<-ii
        try(writeStatas(datas,"daily",select=ii)) # only write the stata files if truly changed!
        ## copy these to output place somewhere determined by defaults-local.txt
        system(paste("cp ",statafile," ",getOption("tempmomo")$all$out$dir,sep=""))
    }
    cat(as.character(Sys.time()),":",i,"END ---------------------\n\n")
}
### Moved the rest of the analysis to explore-weather.R
quit("no")
