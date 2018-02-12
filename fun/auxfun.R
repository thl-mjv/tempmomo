###
### Auxilliary functions for tempmomo
###
### replace missing value in a with b
na.0<-function(a,b=0) ifelse(is.na(a), b,a)
### replace all a==b with missing value
to.na.0<-function(a,b=0) ifelse(a==b    ,NA,a)
### transform farenheiths to celsius assuming 9999.9 missing
to.cels<-function(a) as.numeric(5*(ifelse(a>99.9,NA,a)-32)/9)
### transform dates to years with fractions
numdate<-function(a) as.numeric(a-as.Date("2000-1-1"))/365.25+2000
### number of day in a year
seasday<-function(a) as.numeric(format(as.Date(a),"%j"))/366
### min and max with NA:s removed
namax <-function(a) max (a,na.rm=TRUE)
namin <-function(a) min (a,na.rm=TRUE)
namean<-function(a) mean(a,na.rm=TRUE)
### calculate Monday from a date
monday<-function(a,...) as.Date(a, ...)-na.0(to.na.0(as.numeric(format(as.Date(a,...),"%w"))), 7) + 1
### create a list of full weeks between dates
mondays<-function(first,last,...) { # optional parameters for as.Date
    first<-monday(first,...)
    last <-monday(last ,...)
    length<-as.numeric(last-first)/7
    first+(0:length)*7
}
### create a factor variable of full series of mondays
mondayf<-function(a,first,last,...) { # optional parameters for as.Date
    expected<-mondays(first,last,...)
    observed<-monday(a,...)
    factor(as.numeric(observed),levels=as.numeric(expected),labels=as.character(expected))
}
### give the monday of first (last) full week
fullweek<-function(a,type=c("first","last"),...) {
    a<-as.Date(a,...)
    m<-monday(a)
    type<-match.arg(type)
    if(type=="first") { # truncate to right
        if(a<m) m<-m+7
    } else { # truncate to left
        if(a>m) m<-m-7
    }
    m
}
### create pertty grid of n frames
mkpar<-function(n,row=TRUE) {
    mn<-ceiling(sqrt(n))
    if(row) c(ceiling(n/mn),mn) else c(mn,ceiling(n/mn))
}
### check options
tempoption<-function(country="all",chapter="momo",option="",default=NA) {
    opt<-getOption("tempmomo")
    val<-NULL
    if(country%in%names(opt))
        if(chapter%in%names(opt[[country]]))
            if(option%in%names(opt[[country]][[chapter]]))
                val<-opt[[country]][[chapter]][[option]]
    if(is.null(val)) val<-default
    if(is.numeric(default)) val<-as.numeric(val)
    if(is.logical(default)) val<-as.numeric(val)!=0
    val
}
