###
### Routines for importing data from the NOAA site
### ftp://ftp.ncdc.noaa.gov/pub/data/gsod/
### Data description: ftp://ftp.ncdc.noaa.gov/pub/data/gsod/GSOD_DESC.txt
### Lists of sites  : ftp://ftp.ncdc.noaa.gov/pub/data/gsod/ish-history.csv
###
###
european.countries<-c("AT", "BE", "BG", "CH", "CY", "CZ", "DE", "DK", "EE", "EL", 
                      "ES", "FI", "FR", "HR", "HU", "IE", "IS", "IT", "LI", "LT", "LU", 
                      "LV", "ME", "MK", "MT", "NL", "NO", "PL", "PT", "RO", "SE", "SI", 
                      "SK", "TR", "UK")
#' Load defaults for parameters
#'
#' Loads various defaults files and stores them globally. Simplifies the code as all country specific issues can be resolved in these files.
#' All files with pattern defaults-*.dat will be used in alphapetical order. Thus the latter files override the former. So, using
#' defaults-global.dat and defaults-local.dat the latter will override the former. Also, you can use either .dat or .txt so that .txt overrides .dat
#'
loadDefaults<-function() {
    files1<-list.files(patt="^defaults-.*[.]dat$")
    files2<-list.files(patt="^defaults-.*[.]txt$")
    dats<-unlist(sapply(c(sort(files1),sort(files2)),readLines))
    dats<-dats[!grepl("^#",dats)] # remove comments
    splits<-strsplit(dats,"=")
    labels<-sapply(splits,function(a) strsplit(a[1],"[.]")[[1]])
    values<-sapply(splits,function(a) a[2])
    optmat<-t(rbind(labels,values))
    out<-list()
    for(i in 1:nrow(optmat)) {
        country<-optmat[i,1]
        chapter<-optmat[i,2]
        option <-optmat[i,3]
        if(!country%in%names(out))
            out[[country]]<-list()
        if(!chapter%in%names(out[[country]]))
            out[[country]][[chapter]]<-list()
        out[[country]][[chapter]][[option]]<-optmat[i,4]
    }
    options(tempmomo=out)
    invisible(out)
}
    
#' Download the station specific metadata from the NOAA site
#'
#' Downloads the list of stations, reads it in and cleans the variables
#'
#' @param force should the file be downloaded or use the cached file?
noaaGetSites<-function(force=FALSE) {
    ## todo (28.10.2013): Fix country codes using Ajay's table (done 28.10.2013)
    ## todo (28.10.2013): calculate the NUTS region for each site (done 06.11.2013)
    if(force) {
        sites<-read.table("ftp://ftp@ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.csv",sep=",",header=TRUE)
        names(sites)<-c("usaf","wban","station.name","country","state","call","lat","lon","elev","begin","end")
        sites$begin<-as.Date(as.character(sites$begin),"%Y%m%d")
        sites$end  <-as.Date(as.character(sites$end  ),"%Y%m%d")
        sites$lat  <-with(sites,ifelse(lat == -99999,NA,lat ))
        sites$lon  <-with(sites,ifelse(lon ==-999999,NA,lon ))
        sites$elev <-with(sites,ifelse(elev== -99999,NA,elev))
        sitetrans<-read.table("download/noaa-countries.dat",sep=",",header=TRUE)
        sites$NUTS        <-sitetrans$NUTS   [match(sites$country,sitetrans$NOAA)]
        sites$country.name<-sitetrans$Country[match(sites$country,sitetrans$NOAA)]
        #european.countries<-sort(unique(as.character(subset(sites,!is.na(NUTS))$country))) # this is not useful!
        sites$Europe<-sites$country%in%european.countries
        sites$EuropeProper<-with(sites,(-28<lon) & (lon<35) & (33<lat) & (lat<73)) # this is!
        warning("This sites file does not include proper NUTS region information. Run explore-maps.R")
    } else {
        sites<-read.table("download/sites.dat")
        sites$begin<-as.Date(as.character(sites$begin),"%Y-%m-%d")
        sites$end  <-as.Date(as.character(sites$end  ),"%Y-%m-%d")
    }
    sites
}

#' Match sites to map areas
#'
#' @param sites data frame with lat&lon
#' @param map   ShapeSpatialDataFrame with the regions
#' @param level Level of area
#' @param type  What kind of match should be returned
noaaMapRegions<-function(sites,map,level=0,type=c("inarea","distance","bestquess")) {
    sp<-SpatialPolygons(subset(map,STAT_LEVL_==level)@polygons)
    pcents<-coordinates(sp)
    sites$row<-1:nrow(sites)
    csites<-subset(sites,!is.na(lat) & !is.na(lon))
    csites$x <- csites$lon#/1000
    csites$y <- csites$lat#/1000
    coordinates(csites) <- ~x+y

    e.st<-over(csites,sp,) # primary match
    
    dists <- rdist(pcents,csites@coords)
    dists.min <- apply(dists, 2, which.min)
    
    x.data.country <- subset(x,STAT_LEVL_==level)@data

    inarea<-as.character(x.data.country[e.st,2]) # site residing within the area
    distan<-as.character(x.data.country$NUTS_ID[dists.min]) # by minimun dist
    bestqu<-na.0(inarea,distan) # "best" guess

    ## Shlemiel at work here
    type<-match.arg(type)
    res<-switch(type,
                inarea   =inarea,
                distance =distan,
                bestquess=bestqu)
    res[match(sites$row,csites$row)]
}

##                                                                                          111111111111111111111111111111111111111111111111111
##          11111111112222222222333333333344444444445555555555666666666688888888889999999999000000000011111111112222222222333333333344444444445
## 12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
## STN--- WBAN   YEARMODA    TEMP       DEWP      SLP        STP       VISIB      WDSP     MXSPD   GUST    MAX     MIN   PRCP   SNDP   FRSHTT
## 840010 99999  20130101    74.2  5    68.6  5  1011.9  5  1009.9  5   12.4  5    7.6  5   15.0  999.9    79.5    67.1   0.00I 999.9  000000

#' Downloads data for one station and year
#'
#' The NOAA site has each station's data separated into files. These files are further split into years with one directory for each year.
#' The data is downloaded, read in and transformed to use metric measurement units
#'
#' @param site     the identificator of the station, must be six letter string (i.e. with trailing 0's if shorter)
#' @param year     the year whose data is to be downloaded
#' @param force    should the data be downloaded even if it is found in the cache (default=FALSE)
#' @param thisy    should this years data be downloaded even if force=FALSE
#' @param uplag    how old files should be updated, days
#' @param basepath where the cached files should go
noaaGetSiteYear<-function(site,year=2013,force=FALSE,thisy=TRUE,uplag=0,basepath=NULL) {
    file<-paste(site,"-",year,".op.gz",sep="")
    url <-paste("ftp://ftp.ncdc.noaa.gov/pub/data/gsod/",year,"/",file,sep="")
    if(is.null(basepath)) basepath<-"."
    down<-paste(basepath,"/download/",file,sep="")
    if(file.exists(down)) {
        lastupd<-as.Date(file.info(down)$mtime)
        updsinc<-min(as.numeric(Sys.Date()-lastupd))
        cat("Last updated",as.character(lastupd),",",updsinc,"days ago...")
    } else {
        lastupd<-NA
        updsinc<-uplag+1
        cat("New file, newver updated...")
    }
    ## Need to chech that last date of year is smaller than last update to get through end of year
    if(year==format(Sys.Date(),"%Y") & thisy) force<-(uplag<updsinc)
    if(!file.exists(down) | force) {
        err<-try(download.file(url,down))
        if(inherits(err,"try-err"))
            stop("could not download file")
    }
    data<-try(read.fwf(gzfile(down),header=FALSE,skip=1,as.is=TRUE,
                       widths=diff(c(0,6,12,22,30,33,41,44,52,55,63,66,73,76,83,86,93,100,108,109,116,117,123,124,130,138)),
                       col.names=c("site","wban","yearmoda","temp","tempc","dewp","dewpc","slp","slpc","stp","stpc","visib","visibc",
                                   "wdsp","wdspc","mxspd","gust",
                                   "max","maxf","min","minf","prcp","prcpf","sndp","frshtt")))
    if(!inherits(data,"try-error")) {
        try(data$date <-as.Date(as.character(data$yearmoda),"%Y%m%d"))
        try(data$temp <-to.cels(data$temp))
        try(data$dewp <-to.cels(data$dewp))
        try(data$max  <-to.cels(data$max ))
        try(data$min  <-to.cels(data$min ))
        try(data$slp  <-to.na.0(data$slp,9999.9))
        try(data$stp  <-to.na.0(data$stp,9999.9))
        try(data$visib<-to.na.0(data$visib,999.9)*1.609344)
        try(data$wdsp <-to.na.0(data$wdsp ,999.9)*0.1852)
        try(data$mxspd<-to.na.0(data$mxspd,999.9)*0.1852)
        try(data$gust <-to.na.0(data$gust ,999.9)*0.1852)
        try(data$prcp <-to.na.0(data$prcp ,99.99)*25.4)
        try(data$sndp <-na.0(to.na.0(data$sndp ,999.9)*25.4))
        try(data$frshtt<-sprintf("%06d",data$frshtt))
        try(data$fog  <-as.numeric(substring(data$frshtt,1,1)))
        try(data$rain <-as.numeric(substring(data$frshtt,2,2)))
        try(data$snow <-as.numeric(substring(data$frshtt,3,3)))
        try(data$hail <-as.numeric(substring(data$frshtt,4,4)))
        try(data$thun <-as.numeric(substring(data$frshtt,5,5)))
        try(data$torn <-as.numeric(substring(data$frshtt,6,6)))
        thissite<-subset(sites,sprintf("%06d-%05d",usaf,wban)==site)
        if(nrow(thissite)==0) cat("Site",site,"not found in sites\n")
        try(data$pop0 <-rep(mean(thissite$pop0),nrow(data)))
        try(data$pop3 <-rep(mean(thissite$pop3),nrow(data)))
        try(data$nuts3<-rep(as.character(thissite$NUTS3i[1]),nrow(data)))
        attr(data,"downloaded")<-na.0(lastupd,Sys.Date())
    } else {
        ## empty structure
        data<-structure(list(site = integer(0), wban = integer(0), yearmoda = integer(0), 
                             temp = numeric(0), tempc = integer(0), dewp = numeric(0), 
                             dewpc = integer(0), slp = numeric(0), slpc = integer(0), 
                             stp = numeric(0), stpc = integer(0), visib = numeric(0), 
                             visibc = integer(0), wdsp = numeric(0), wdspc = integer(0), 
                             mxspd = numeric(0), gust = numeric(0), max = numeric(0), 
                             maxf = character(0), min = numeric(0), minf = character(0), 
                             prcp = numeric(0), prcpf = character(0), sndp = numeric(0), 
                             frshtt = character(0), date = structure(numeric(0), class = "Date"), 
                             fog = numeric(0), rain = numeric(0), snow = numeric(0), hail = numeric(0), 
                             thun = numeric(0), torn = numeric(0), pop0 = numeric(0), 
                             pop3 = numeric(0), nuts3 = character(0)),
                        .Names = c("site", 
                                   "wban", "yearmoda", "temp", "tempc", "dewp", "dewpc", "slp", 
                                   "slpc", "stp", "stpc", "visib", "visibc", "wdsp", "wdspc", "mxspd", 
                                   "gust", "max", "maxf", "min", "minf", "prcp", "prcpf", "sndp", 
                                   "frshtt", "date", "fog", "rain", "snow", "hail", "thun", "torn", 
                                   "pop0", "pop3", "nuts3"), downloaded = 17302, row.names = integer(0), class = "data.frame")
    }
    data
}

#' Get all data from all stations in a country
#'
#' Go through all stations in the country by each year and download the data
#'
#' @param usecountry the name of the country as listed in ftp://ftp.ncdc.noaa.gov/pub/data/gsod/country-list.txt
#' @param years vector of years to use
#' @param sitesdata a data frame created using noaaGetSites
#' @param force Should the data be downloaded even if it exists in the cache
#' @param thisy Should this year's potentially incomplete data be downloaded even if it exists in the cache
#' @param basepath Directory where the data is downloaded (in subdirectory download)
#' @param countrycode What variable in the sites -file gives the country name
noaaGetCountry<-function(usecountry="FI",years=2008:2013,sitesdata=sites,force=FALSE,thisy=TRUE,basepath=NULL,uplag=0,
                         countrycode=c("country","NUTS","NUTS0i","NUTS3i")) {
    ## todo (xx.05.2013): error checking for missing sitesdata
    ## todo (28.10.2013): allow selection of country using either fields country, NUTS, country.name
    countrycode<-match.arg(countrycode)
    csites<-try(subset(sitesdata,get(countrycode)==usecountry))
    if(inherits(csites,"try-error")) {
        sitesdata<-noaaGetSites()
        csites<-try(subset(sitesdata,get(countrycode)==usecountry))
    }
    cat(nrow(csites),"stations\n")
    res<-list()
    for(i in years) {
        ## kludge: if either begin or end is missing, use anyways
        ysites<-with(subset(csites,na.0(format(begin,"%Y"),"9999")<=i & na.0(format(end,"%Y"),"0") >= i),sprintf("%06d-%05d",usaf,wban))
        cat(length(ysites),"stations year",i,"\n")
        yres<-list()
        for(j in ysites) {
            cat(usecountry,i,j,which(j==ysites),length(ysites),"...")
            cat(system.time(yres[[j]]<-noaaGetSiteYear(j,i,force=force,thisy=thisy,basepath=basepath,uplag=uplag))[3])
            cat(" seconds to download\n")
        }
        res[[i]]<-do.call("rbind",yres)
    }
    out<-do.call("rbind",res)
    last.date<-try(as.Date(max(as.character(out$date),na.rm=TRUE)))
    print(last.date)
    if(inherits(last.date,"try-error")) last.date<-NULL
    if(is.null(last.date)) last.date<-Sys.Date()
    if(is.na(last.date)) last.date<-Sys.Date()
    print(last.date)
    try(attr(out,"downloaded")<-last.date) # collect the date of last update. Should use max from yres?
    out
}

#' Get already downloaded data
#'
#' The data download takes some time, so we'll save it at data/weather.rda (data.R)  or data/weather-global.rda (data-global.R)
#'
#' @param local should we use global or local data?
loadWeatherData<-function(local=TRUE) {
    if(local) file<-"data/weather.rda" else file<-"data/weather-global.rda"
    if(!file.exists(file)) {
        if(local) stop("Local file not found, run data.R first")
        res<-try(download.file("https://owncloud.thl.fi/public.php?service=files&t=a3512ac71e9ab593421171843af081c5&download&path=//weather-global.rda",down))
        ## TODO (7.11.2013: if local=FALSE fails, try local=TRUE
        if(inherits(res,"try-error"))
            stop("Global file not found, try local")
    } # after this file should exist!
    load(file,.GlobalEnv)
    invisible(NULL)
}
#' Aggregate one weather dataset
#'
#' Currently only mean over 
#'
#' @param data     data to aggregate
#' @param stations which stations to use, default=NULL means all
#' @param vars     which variables to aggregate
#' @param sfun     how to aggregate stations
#' @param dfun     how to aggregate days to weeks
#' @param last     last date to use
#' @param first    first date to use
#' @param weights  weights used in aggregation across stations
weatherAggre<-function(data,stations=NULL,vars=c("temp"),sfun=list("All"=mean),dfun=list("All"=mean),last=NULL,first=NULL,weights=NULL) {
    lstat<-sort(unique(data$site))
    if(is.null(stations)) stations<-list()
    if(is.null(stations$All)) stations$All<-lstat
    if(is.null(names(vars))) names(vars)<-vars
    if(is.null(first)) first<-fullweek(min(as.character(data$date)),"first")
    if(is.null(last )) last <-fullweek(max(as.character(data$date)),"last" )
    out<-data.frame(date=mondays(first,last))
    for(i in names(vars)) {
        cat(i,"...")
        stats<-stations[[i]]
        if(is.null(stats)) stats<-stations$All
        data0<-subset(data,site%in%stats)
        cat(nrow(data0),"...")
        if(nrow(data0)==0) warning("No eligible stations, skipping")
        else {
            ## fixed: instead of get() use eval(parse()) so that you can use expressions
            if(!i%in%names(sfun)) lsfun<-sfun$All else lsfun<-sfun[[i]]
            if(!i%in%names(weights)) wgt<-"1" else wgt<-weights[[i]]
            cat(expr0<-paste("(",wgt,")*(",as.character(vars[[i]]),")",sep=""))
            res<-with(data0,tapply(eval(parse(text=expr0)),date,lsfun)) # aggregate         by stations
            if(wgt!="1")
            wes<-with(data0,tapply(eval(parse(text=wgt  )),date,lsfun)) # aggregate weights by stations
            else wes<-1
            if(!i%in%names(dfun)) ldfun<-dfun$All else ldfun<-dfun[[i]]
            res<-tapply(res/wes,mondayf(names(res),first,last),ldfun)
            out[[i]]<-res
        }
        cat("\n")
    }
    out
}
        
#' Load mortality data
#'
#' Loads mortality data from different sources using subfunctions for each case
#'
#' @param country name of the country
#' @param type    type of the source data
loadMortData<-function(country,type=c("txt","a-momo","stata"),...) {
    if(missing(type))
        type<-getOption("tempmomo")[[country]]$mort$type
    if(is.null(type)) type<-"txt"
    type<-match.arg(type)
    res<-NULL
    if(type=="txt")    res<-loadMortDataTXT (country,...)
    if(type=="a-momo") res<-loadMortDataMOMO(country,...)
    if(type=="stata")  res<-loadMortDataStata(country,...)
    if(is.null(res)) stop("Could not load data") # should never happen as "type" is matched
    if("group"%in%names(res))
        ord<-order(res$group,res$date)
    else
        ord<-order(res$date)
    res[ord,]
}
#' Load mortality data from text file
#'
#' Loads and checks mortality data from a text file
#'
#' @param country   name of the country
#' @param extension file extension for the data
#' @param start     start of the calendar for numeric dates
#' @param vars      vector of variable names
loadMortDataTXT<-function(country,extension=NULL,start=NULL,datetype=NULL,
                               vars=c("n"="n","date"="date","pop"="pop")) {
    defaults<-getOption("tempmomo")
    if(is.null(extension)) extension<-defaults[[country]]$mort$ext
    if(is.null(extension)) extension<-"dat"
    if(is.null(start    )) start    <-defaults[[country]]$mort$start
    if(is.null(start    )) start    <-"1970-1-1"
    if(is.null(datetype )) datetype <-defaults[[country]]$mort$datetype
    if(is.null(datetype )) datetype <-"%Y-%m-%d"
    file<-paste("download/mort_",country,".",extension,sep="")
    if(!file.exists(file)) {
        warning(paste("File",file,"not found"))
        return(NULL)
    }
    data<-read.table(file,sep=",",header=TRUE)
    ## checking.
    remove<-NULL
    ## TODO: what if the data has only year/week?
    if(datetype=="week") {
        if(!require("ISOweek")) stop("You need package ISOweek installed to use A-MOMO inputs")
        starty<-min(data$year,na.rm=TRUE)
        endy  <-max(data$year,na.rm=TRUE)
        ## list all weeks you can have
        calendar<-mondays(paste(starty,"-1-1",sep=""),paste(endy,"-12-31",sep=""))
        calendarweek<-gsub("W","",ISOweek(calendar))
        data$date<-calendar[match(sprintf("%04d-%02d",data$year,data$week),calendarweek)]
        remove<-c(remove,"year","week")
    }
    for(i in c("n","date","pop")) { # for each required field
        if(vars[i]%in%names(data)) { # does it exist in the data?
            if(vars[i]!=i) { # and has different name
                data[[i]]<-data[[vars[i]]] # replace
                remove<-c(remove,vars[i])
            }
        }
        if(i=="date") {
            if(!inherits(data$date,"Date")) {
                if(!is.character(data$date) & !is.numeric(data$date))
                    try(data$date<-as.character(data$date))
                if(is.character(data$date))
                    try(data$date<-as.Date(data$date,format=datetype))
                if(is.numeric(data$date))
                    try(data$date<-as.Date(data$date,start))
            }
        }
    }
    for(i in remove)
        data[i]<-NULL
    data
}
#' Load mortality data from A-MOMO output
#'
#' Loads and checks mortality data from A-MOMO output, removing unneeded variables
#'
#' @param country     Name of the country (short name)
#' @param countryname Name of the country (full name)
#' @param root        Root directory for the A-MOMO outputs
#' @param year        Year number for the data to use, default latest
#' @param week        Week number for the data to use, default latest
loadMortDataMOMO<-function(country,countryname=NULL,root=NULL,year=NULL,week=NULL) {
  if(!require("ISOweek")) stop("You need package ISOweek installed to use A-MOMO inputs")
    defaults<-getOption("tempmomo")
    if(is.null(countryname)) countryname<-defaults[[country]]$momo$country
    if(is.null(countryname)) countryname<-country
    if(is.null(root       )) root       <-defaults[[country]]$momo$root
    if(is.null(root       )) root       <-""
    if(!file.exists(root)) {
        warning("could not find root directory")
        return(NULL)
    }
    if(is.null(year)) year<-format(Sys.Date(),"%Y")
    pattern<-paste("MOMOv4-3-",countryname,"-",year,"-",week,sep="")
    files<-list.files(path=root,full.names=TRUE,include.dirs=TRUE,pattern=pattern)
    rest<-strsplit(substring(files,2+nchar(root)+9+nchar(countryname)+6),"-")
    weekstring<-sapply(rest,function(a) as.numeric(a[1]))
    additional<-sapply(rest,function(a) paste(a[-1],collapse="-"))
    momoroot<-rev(files[order(weekstring,additional)])[1] # take the latest
    if(!file.exists(momoroot)) {
        warning("could not find momoroot directory")
        return(NULL)
    }
    subroot<-list.files(momoroot,patt="EUROMOMO-COMPLETE",include.dirs=TRUE,full.names=TRUE)
    print(datafiles<-list.files(subroot,patt="txt",full.names=TRUE))
    data<-read.table(datafiles[1],sep=",",header=TRUE,na="",as.is=TRUE)
    starty<-min(data$YoDi,na.rm=TRUE)
    endy  <-max(data$YoDi,na.rm=TRUE)
    ## list all weeks you can have
    calendar<-mondays(paste(starty,"-1-1",sep=""),paste(endy,"-12-31",sep=""))
    calendarweek<-gsub("W","",ISOweek(calendar))
    data<-data[,c("group","wk2","excess","Pnb","nbc","nb","zscore")]
    names(data)<-c("group","wk2","excess","fit","nc","n","Z")
    data$date<-calendar[match(data$wk2,calendarweek)]
    data
}
#' Load mortality data from a Stata file
#'
#' Loads and checks mortality data from a Stata file compatible with FluMoMo code distributed by Jens in Spring 2013
#'
#' @param country     Name of the country (short name)
loadMortDataStata<-function(country) {
  if(!require("ISOweek")) stop("You need package ISOweek installed to use A-MOMO inputs")
    file<-paste("download/mort_",country,".dta",sep="")
    if(!file.exists(file)) {
        warning(paste("File",file,"not found"))
        return(NULL)
    }
    require("foreign")
    data<-read.dta(file)
    starty<-min(data$year,na.rm=TRUE)
    endy  <-max(data$year,na.rm=TRUE)
    ## list all weeks you can have
    calendar<-mondays(paste(starty,"-1-1",sep=""),paste(endy,"-12-31",sep=""))
    calendarweek<-gsub("W","",ISOweek(calendar))
    data<-data[,c("week","year","deaths","agegrp")] # remove possible excess variables
    names(data)<-c("week","year","n","group")
    data$date<-calendar[match(sprintf("%04d-%02d",data$year,data$week),calendarweek)]
    data
}
                      
#' Write data to Stata
#'
#' Writes a dataset to Stata with necessary metadata
#'
#' @param data a list of data frames
#' @param name optionally names of the elements to output
#' @param stata path to stata executable
writeStatas<-function(data,name=NULL,select=NULL,stata=NULL,basepath="./") {
    require("foreign")
    if(is.null(name)) name<-names(data)
    if(is.null(select)) select<-names(data)
    if(is.null(stata)) stata<-"/usr/local/stata/stata" # TODO: use options!
    if(!file.exists(stata)) stata<-NULL
    for(i in select) {
        if(i %in% names(data)) {
            codefile<-paste(basepath,"/out/",name,"-",i,".base.do",sep="")
            datafile<-paste(basepath,"/out/",name,"-",i,".dat",sep="")
            dodofile<-paste(basepath,"/out/",name,"-",i,".do",sep="")
            cat(i,codefile, datafile,"... ")
            names(data[[i]])<-gsub("[.]","_",names(data[[i]]))
            cls<-sapply(data[[i]],data.class)
            dts<-names(cls[cls=="Date"])
            for(j in dts) {
                data[[i]][[j]]<-as.numeric(data[[i]][[j]]-as.Date("1960-1-1"))
            }
            write.foreign(data[[i]],
                          datafile=datafile,
                          codefile=codefile,
                          package="Stata")
            cat("set mem 1000m\n",file=dodofile,append=FALSE,sep="")
            cat(readLines(codefile),"\n",
                file=dodofile,append=TRUE,sep="")
            cat("label data \"data for ",i,", uploaded at ",format(Sys.time()),"\"\n",
                file=dodofile,append=TRUE,sep="")
            for(j in dts) {
                cat("format ",j," %td\n",file=dodofile,append=TRUE,sep="")
            }
            ## TODO: variable labels
            cat("save  ",basepath,"/out/",name,"-",i,",replace\n",
                file=dodofile,append=TRUE,sep="")
            if(!is.null(stata)) {
                cmd<-paste(stata," -bq \"do ",dodofile,";exit\"", sep = "")
                cat(cmd,"\n")
                system(cmd)
            } else {
                cat("Stata not found!\n")
            }
        }
    }
}

