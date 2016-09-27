################################################################################
###
### explore-variograms.R
###
### Explore the variograms of the weather data downloaded in data.R
###
### Source utility functions
source("fun/auxfun.R")
source("fun/data.R")

### Load the data saved in data.R
loadWeatherData(FALSE) # FALSE: use global data

basecol<-colorRampPalette(1:8)(length(datas))
### define some colors for further use
(colors<-rgb(t(col2rgb(basecol)),maxColorValue=255,alpha=255)) # dark
(shades<-rgb(t(col2rgb(basecol)),maxColorValue=255,alpha=25 )) # transparent
(shadeS<-rgb(t(col2rgb(basecol)),maxColorValue=255,alpha=75 )) # slightly transparent
names(colors)<-names(shades)<-names(shadeS)<-sort(names(datas))

wvars<-c("temp","dewp","slp","stp","visib","wdsp","mxspd","gust","max","min","prcp","sndp")

## What is the spatio temporal variogram?
## Start with spatial, temporal is less interesting although the seasonality may need to be adjusted for
## Is it meaningful to look into variograms ithin country, or should we calculate the variograms across all sites?
## Also, lat/lon distances need to be calculated differently! (I assume)

## Install this package first e.g.
## install.packages("geoR")
library(geoR) 

for(i in names(datas)) {
    print(summary(datas[[i]]$lat<-with(sites,lat[match(datas[[i]]$site,sites$usaf)]/10)))
    print(summary(datas[[i]]$lon<-with(sites,lon[match(datas[[i]]$site,sites$usaf)]/10)))
}

countries<-names(datas)                         # all countries
#countries<-c("fi","hu","pt","ie","mt","dk")     # some countries
#countries<-countries[countries%in%names(datas)] # make sure all countries have data
#countries<-NULL                                 # no recalculations

if(file.exists("data/variograms.rda")) {
    cat("Variogram data already exists, using it\n")
    load("data/variograms.rda")
    ## countries<-countries[!countries%in%names(variogs)] (there may be missing subcalculations)
} else {
    cat("No variogram data found\n")
    variogs<-list() # DRY!
}

for(i in countries) {
    cat(i,"----------------------")
    if(!i %in% names(variogs)) {
        variogs[[i]]<-list()        
        cat("NEW\n")
    } else cat("OLD\n")
    for(j in wvars) {
        cat(i,j,"----------------")
        if(!j %in% names(variogs[[i]])) {
            variogs[[i]][[j]]<-list()
            cat("NEW\n")
        } else cat("OLD\n")            
        for(m in c(sprintf("%02d",1:12))) { # Could not evaluate full year for some countries
            cat(i,j,m,"----------")
            if(!m %in% names(variogs[[i]][[j]])) {
                cat("NEW\n")
                print(system.time(vg<-try(with(subset(datas[[i]],
                                                      (format(date,"%m")==m | m=="All") & !is.na(get(j))),
                                               variog(coords=cbind(lat,lon),data=get(j),
                                                      breaks=c(0,10,20,40,80,160)))))) # upto 320km
                if(!inherits(vg,"try-error")) variogs[[i]][[j]][[m]]<-vg
            } else cat("OLD\n")            
        }
    }
}

save(variogs,file="data/variograms.rda")

png("out/eu_vg_all.png",width=1200,height=800,pointsize=10)
par(mfcol=c(3,4))
for(i in wvars) {
    with(variogs[[1]][[1]][[1]],plot(u,0*u,type="n",main=i,ylim=c(-1,1)))
    for(j in names(variogs))
        for(k in 1:12)
            try(with(variogs[[j]][[i]][[k]],lines(u,(1-v/var.mark),type="l",lty=1,col=shadeS[j])))
}
dev.off()
png("out/eu_vg_temp.png",width=1200,height=800,pointsize=10)
par(mfrow=c(3,4))
for(k in 1:12) {
    with(variogs[[1]][[1]][[1]],plot(u,0*u,type="n",main=k,ylim=c(-1,1),
                                     xlim=c(0,max(u,na.rm=TRUE)*1.1)))
    for(j in names(variogs))
        try(with(variogs[[j]][["temp"]][[k]],{
            lines(u,(1-v/var.mark),type="l",lty=1,col=colors[j])
            text(u[length(u)],(1-v/var.mark)[length(u)],labels=j,col=colors[j],adj=0)
        }))
}
    
dev.off()

