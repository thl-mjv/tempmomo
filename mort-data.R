###
### load mortality data
### 
source("fun/auxfun.R")
source("fun/data.R")
source("fun/analysis.R")
### load defaults so that no need for country specific code
loadDefaults()
### 
morts<-list()
for(i in tolower(european.countries)) {
    cat(i,"...")
    tmp<-try(loadMortData(i)) # using default type
    for(j in c("txt","a-momo","stata")) 
      if(inherits(tmp,"try-error")) {
        cat(j,"...")
        tmp<-try(loadMortData(i,type=j))
      }
    if(!is.null(tmp) &!inherits(tmp,"try-error")) {
        cat("Save\n")
        morts[[i]]<-tmp
    } else {
        cat("Fail\n")
    }
}
names(morts)

sapply(morts,summary)

png(file="out/eu_mort.png",width=1000,height=250*length(morts))
par(mfcol=c(length(morts),1))
for(i in names(morts)) {
    with(morts[[i]],plot(date,n,type="l",main=i))
}
dev.off()

save(morts,file="data/mort-raw.rda")
