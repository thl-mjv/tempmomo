################################################################################
###
### data.R
###
### Import data from the NOAA site - local sites
###
### Source utility functions
source("fun/auxfun.R")
source("fun/data.R")
if(file.exists("data/mort-raw.rda")) load("data/mort-raw.rda")
if(!exists("morts")) morts<-list()
### Load sites metadata
sites<-noaaGetSites(FALSE)
### A bit of existentialism: What is Europe

with(sites,table(Europe,EuropeProper))

### Load data. This will take time, when tested about 1500s for the first time and about 120s for reruns with no downloading
###
### Add a row for your own country. Note that the NOAA site uses slightly different country codes. See download/noaa-countries.dat for translations
### 
datas<-list()
use.this.year<-TRUE # should the current year be included and updated. If not, then this file needs to be updated only once a year
(years<-as.numeric(format(Sys.Date(),"%Y"))-6:(!use.this.year))
countries<-names(morts)
if(is.null(countries)) countries<-"fi" # replace with your own country code if needed
### you might want to use
# countries<-"your country" 
### here if you only want to use your own data
for(i in countries) {
    print(system.time(datas[[i]]<-noaaGetCountry(toupper(i),force=FALSE,years=years,
                                                 thisy=use.this.year,countrycode="NUTS0i")) )
}

### Check the results
summary(datas[[1]])
head(datas[[1]])

weather.data.version<-Sys.time()

### Save the results for further use
save(sites,datas,weather.data.version,file="data/weather.rda") # local data

### Moved the rest of the analysis to explore-weather.R
