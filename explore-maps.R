################################################################################
###
### explore-weather.R
###
### Explore the weather data downloaded in data.R
###
### Source utility functions
source("fun/auxfun.R")
source("fun/data.R")
## Install this package first e.g.
## install.packages("geoR")
library(geoR) 
## install these packages also
library(sp)
library(maptools)
library(fields)

sites<-noaaGetSites(!interactive())
#sites<-noaaGetSites(interactive())
#tmp<-readLines("ftp://ftp@ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.csv")
#head(tmp)
summary(sites)
head(sites)
x <- readShapeSpatial("download/NUTS_RG_03M_2010.shp")
                                        
##country codes
par(mfcol=c(2,2))
plot(x )# OK, these are only european NUTS areas
plot(subset(x,STAT_LEVL_==0)) # countries
plot(subset(x,STAT_LEVL_==2)) 
plot(subset(x,STAT_LEVL_==3)) 

## countries
(countries<-names(with(subset(x@data,STAT_LEVL_==0),table(as.character(NUTS_ID)))))
par(mfcol=c(1,1))
plot(subset(x,STAT_LEVL_==0),xlim=c(-28,35),ylim=c(33,73))
points(coordinates(subset(x,STAT_LEVL_==0)),col="red",cex=3)
text(coordinates(subset(x,STAT_LEVL_==0)),labels=as.character(subset(x,STAT_LEVL_==0)$NUTS_ID),col="red")
with(sites,points(lon,lat,pch=16,cex=1,col="green"))

table(sites$NUTS0i<-noaaMapRegions(sites,x,0,"i"))
table(sites$NUTS0d<-noaaMapRegions(sites,x,0,"d"))
table(sites$NUTS0b<-noaaMapRegions(sites,x,0,"b"))
table(sites$NUTS3i<-noaaMapRegions(sites,x,3,"i"))
table(sites$NUTS3d<-noaaMapRegions(sites,x,3,"d"))
table(sites$NUTS3b<-noaaMapRegions(sites,x,3,"b"))

subset(as.data.frame(with(sites,table(country,NUTS0i))),Freq>0)

(CountryNUTS<-names(table(as.character(sites$NUTS0i))))
(CountryCols1<-apply(col2rgb(sample(rainbow(length(CountryNUTS))))/255,2,function(a) rgb(a[1],a[2],a[3],0.1)))
(CountryCols2<-apply(col2rgb(        CountryCols1                )/255,2,function(a) rgb(a[1],a[2],a[3],1.0)))
names(CountryCols1)<-names(CountryCols2)<-CountryNUTS

### sites in a country by type
par(mfcol=c(1,1))
plot(subset(x,STAT_LEVL_==0),xlim=c(-28,35),ylim=c(33,73))
for(i in CountryNUTS) {
    plot  (subset(x,NUTS_ID==i&STAT_LEVL_==0),col=CountryCols1[i],add=TRUE)
    with(subset(sites ,NUTS0i==i&EuropeProper),points(lon,lat,col=CountryCols2[i],pch=16))
    with(subset(sites ,NUTS0b==i&EuropeProper),points(lon,lat,col=CountryCols2[i],pch= 1))
}

### sites not in a country
par(mfcol=c(1,1))
plot(subset(x,STAT_LEVL_==0),xlim=c(-28,35),ylim=c(33,73))
with(subset(sites,(!NUTS0i%in%CountryNUTS) & EuropeProper),points(lon/1000,lat/1000,col=rgb(1,0,0,.5),,pch=16))
with(subset(sites,( NUTS0i%in%CountryNUTS) & EuropeProper),points(lon/1000,lat/1000,col=rgb(0,0,1,.5),,pch=16))

### All in all we should only use NUTS0i
### Otherwise we mut calculate the distance to the closest point.
### In any case, the North Sea/Atlantic sites will be missassigned
### Some of the offshore sites very close to the shore will be missed because of
###  lat/lon != ETRS89

### What is the distribution of the distances between sites within each country
### Problem: diagonal distances are not right due to the cöordiante systems
(function(a) a[order(a[,"Median"]),])(
    t(sapply(CountryNUTS,function(a) with(subset(sites,NUTS0i==a),summary(dist(cbind(lat,lon)/1000)))))
)
### What is the distribution of the distances between sites across the countries
with(subset(sites,EuropeProper),summary(dist(cbind(lat,lon)/1000)))

### how many sites within radius of 1000?
with(subset(sites,Europe&EuropeProper),tapply(apply(as.matrix(dist(cbind(lat,lon)))<1000,1,sum),NUTS0b,mean))
with(subset(sites,Europe&EuropeProper),cbind(apply(as.matrix(dist(cbind(lat,lon)))<1000,1,sum),as.character(NUTS0b))[order(as.character(NUTS0b)),])

### We need population
pop<-read.table("download/demo_r_gind3_1_Data.csv",sep=",",header=TRUE)
table(pop$country<-substring(as.character(pop$GEO),1,2))
table(pop$stat.level<-nchar(as.character(pop$GEO))-2)
table(pop$is.pop<-substring(as.character(pop$INDIC_DE),1,3)=="Pop")
summary(pop$value<-as.numeric(gsub(":","NA",gsub(" ","",as.character(pop$Value)))))
summary(pop)
with(sites,table(EuropeProper,as.character(NUTS3i)%in%as.character(pop$GEO)))

### do these match with the codes in the map?
plot(subset(x,STAT_LEVL_==0),xlim=c(-28,35),ylim=c(33,73))
with(subset(sites,EuropeProper),points(lon/1000,lat/1000,col=ifelse(!as.character(NUTS3b)%in%as.character(pop$GEO),"red","blue"),pch=16))

### OK, there's mixed NUTS1, NUTS2 and NUTS3
with(subset(pop,country=="AT" & is.pop),tapply(value,stat.level,sum))
with(subset(pop,is.pop & TIME==2012),tapply(value,list(country,stat.level),sum)) # some countries are missing level 3

### do these match with the codes in the map
par(mfcol=c(2,2))
for(i in 0:3) {
    plot(subset(x,STAT_LEVL_==0),xlim=c(-28,35),ylim=c(33,73))
    with(subset(sites,EuropeProper),
         points(lon/1000,lat/1000,col=ifelse(!as.character(NUTS3b)%in%as.character(subset(pop,stat.level==i)$GEO),"red","blue"),pch=16))
}
### matching
summary(sites$pop0<-with(subset(pop,TIME==2012 & stat.level==0 & is.pop),value[match(as.character(sites$NUTS0b),as.character(GEO))]))
summary(sites$pop3<-with(subset(pop,TIME==2012 & stat.level==3 & is.pop),value[match(as.character(sites$NUTS3b),as.character(GEO))]))

with(sites,tapply(pop0,NUTS0b,mean))
with(sites,tapply(pop3,NUTS3b,mean))

head(sites)

### OK, let's save the file
write.table(sites,file="download/sites.dat")

################################################################################
###
### problems with missing sites
###
par(mfcol=c(1,1))
plot(subset(x,STAT_LEVL_==0 & NUTS_ID=="PT"))#,xlim=c(-28,35),ylim=c(33,73))
with(subset(sites,(!NUTS0i%in%"PT") & EuropeProper),points(lon/1000,lat/1000,col=rgb(1,0,0,.5),,pch=16))
with(subset(sites,( NUTS0i%in%"PT") & EuropeProper),points(lon/1000,lat/1000,col=rgb(0,0,1,.5),,pch=16))


loadWeatherData(FALSE) 
#loadWeatherData(TRUE ) 
names(datas)
eumeans<-do.call("c",lapply(datas,function(a) with(a,tapply(temp,site,mean,na.rm=TRUE))))
eumins <-do.call("c",lapply(datas,function(a) with(a,tapply(min ,site,min ,na.rm=TRUE))))
eumaxs <-do.call("c",lapply(datas,function(a) with(a,tapply(max ,site,max ,na.rm=TRUE))))
table(substring(names(eumeans),4)%in%sites$usaf)
summary(sites$mean.temp<-eumeans[match(sites$usaf,substring(names(eumeans),4))])
summary(sites$ max.temp<-eumaxs [match(sites$usaf,substring(names(eumaxs ),4))])
summary(sites$ min.temp<-eumins [match(sites$usaf,substring(names(eumins ),4))])

par(mfcol=c(1,1))
plot(subset(x,STAT_LEVL_==0),xlim=c(-28,35),ylim=c(33,73))
with(subset(sites,EuropeProper),points(lon/1000,lat/1000,col=ifelse(min.temp<0,rgb(0,0,1,.25),rgb(1,0,0,.25)),pch=16,cex=abs(min.temp)/20))
with(subset(sites,EuropeProper),points(lon/1000,lat/1000,col=ifelse(max.temp<0,rgb(0,0,1,.25),rgb(1,0,0,.25)),pch=16,cex=abs(max.temp)/20))
with(subset(sites,EuropeProper),points(lon/1000,lat/1000,col=ifelse(is.na(mean.temp),rgb(1,1,0,.25),rgb(0,0,0,.25)),pch=16))
head(sort(eumeans))
subset(sites,usaf==161550)
subset(datas$it,site==161550)
head(subset(sites,!is.na(mean.temp)))

lapply(2008:2014,function(a) with(subset(sites,format(begin,"%Y")<=a & format(end,"%Y")>=a),
                                  table(NUTS0i=="DE",usaf%in%with(datas$de,sort(unique(site))))))
with(subset(sites,NUTS0i=="DE" & !usaf%in%with(datas$de,sort(unique(site)))),summary(begin))
