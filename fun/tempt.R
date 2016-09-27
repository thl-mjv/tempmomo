################################################################################
###
### Code for tempt classes
### TODO: 4th dimension (groups such as age)
tempt<-function(x=NULL,date=NULL,region=NULL,name=NULL,data=globalenv(),sitesfile=sites,fun=namean) {
  if(is.null(name)) name<-deparse(substitute(x))
  
  if(is.data.frame(x)) {
        data<-x
        x<-NULL
    }
    
    ## date variable
    if(is.null(date)) date<-attr(data,"date")
    if(is.null(date)) date<-"date"
    if(is.character(date) & length(date)==1) # varname
        date<-data[[date]]
    ## date must be a real date
    date<-as.Date(date)
    ## it must have a frequency (common for all)
    dfreq<-as.numeric(diff(sort(date)))
    if(length(table(dfreq))>1) warning("Irregular times")
    freq<-min(dfreq,na.rm=TRUE)

    ## region variable
    if(is.null(region)) region<-attr(data,"region")
    if(is.null(region)) region<-"region"
    if(is.character(region) & length(region)==1) # varname
        region<-data[[region]]
    ## reg must be a site or other areacode
    regtype<-"UNK"
    if(all(region%in%sitesfile$usaf)) {
        regtype<-"sites"
    }
    if(all(region%in%sitesfile$NUTS0b)) {
        regtype<-"NUTS0"
    }
    if(all(region%in%sitesfile$NUTS3b)) {
        regtype<-"NUTS3"
    }
    if(regtype=="UNK") warning("Unknown region type")

    ## check!
    tt<-table(date,region)
    if(any(tt>1)) warning("Non unique data, using sum to aggregate")

    ### response
    if(is.null(x)) x<-attr(data,"vars")
    if(is.null(x)) {
        x<-names(data)
        x<-x[!x%in%c("region","date")]
    }
    if(is.character(x)) {# varname
        name<-x
        x<-unclass(data[x])
    }
    if(!is.list(x)) x<-list(x)
    names(x)<-name

    dn<-c(dimnames(tt),list(name))
    names(dn)<-NULL # for all.equal
    out<-array(NA,dim=c(dim(tt),length(name)),dimnames=dn)

    for(i in name) {
        out[,,i]<-tapply(x[[i]],list(date,region),fun) #allowing different functions speed the calculations
    }
    ## frequency must be retained in attributes
    attr(out,"freq")<-freq
    ## areacode type must be retained in attributes
    attr(out,"regtype")<-regtype
    class(out)<-"tempt"
    out
} 
tfreq<-function(a) UseMethod("tfreq")
tfreq.tempt<-function(a) attr(a,"freq")
tname<-function(a) UseMethod("tname")
tname.tempt<-function(a) dimnames(a)[[3]]
ttime<-function(a) UseMethod("ttime")
ttime.tempt<-function(a) dimnames(a)[[1]]
tdate<-function(a) as.Date(ttime(a))
tregs<-function(a) UseMethod("tregs")
tregs.tempt<-function(a) dimnames(a)[[2]]
### TODO: arbitrary number of tempts:
combine<-function(t1,t2,replace=c("left","right")) {
  ## todo multiple t:s
  f1<-tfreq(t1)
  f2<-tfreq(t2)
  if(f1!=f2) stop("missmatching frequencies, aggregate first")
  r1<-attr(t1,"regtype")
  r2<-attr(t2,"regtype")
  if(r1!=r2) stop("missmatching regions, aggregate first")
  dts <-sort(unique(c(ttime(t1),ttime(t2))))
  regs<-sort(unique(c(tregs(t1),tregs(t2))))
  nas <-sort(unique(c(tname(t1),tname(t2))))
  out<-array(NA,dim=c(length(dts),length(regs),length(nas)),dimnames=list(dts,regs,nas))
  replace<-match.arg(replace)
  if(replace=="right") {
    out[match(ttime(t1),dts),match(tregs(t1),regs),match(tname(t1),nas)]<-t1
    out[match(ttime(t2),dts),match(tregs(t2),regs),match(tname(t2),nas)]<-t2
  }
  if(replace=="left") {
    out[match(ttime(t2),dts),match(tregs(t2),regs),match(tname(t2),nas)]<-t2
    out[match(ttime(t1),dts),match(tregs(t1),regs),match(tname(t1),nas)]<-t1
  }
  attr(out,"freq")<-f1
  attr(out,"regtype")<-r1
  class(out)<-"tempt"
  out
}
taggre<-function(a,fun=namean,agg=monday,weights=NULL) {
  ntime<-agg(ttime(a))
  if(is.null(weights)) weights<-rep(1,dim(a)[1])
  if(!is.null(names(weights))) {
      weights<-weights[match(ttime(a),names(weights))]
  } else {
      weights<-rep(weights,length=dim(a)[1])
  }
  nfreq<-min(as.numeric(diff(sort(unique(ntime)))))
  out<-apply(a,2:3,function(x) tapply(x*weights,ntime,fun)/fun(weights))
  attr(out,"freq")<-nfreq
  attr(out,"regtype")<-attr(a,"regtype")
  class(out)<-"tempt"
  out
}
ttransf<-function(a,fun=deseas,...) {
  out<-apply(a,2:3,fun,date=as.Date(ttime(a)),...)
  attributes(out)<-attributes(a)
  out
}
deseas<-function(x,date,type=c("spline","sin"),...) {
  type<-match.arg(type)
  out<-NULL
  if(type=="spline")
    fit<-gam(x~s(seasday(date),bs="cc",...))
  if(type=="sin") 
    fit<-lm(x~momosin(date,...))
  pred<-rep(NA,length(x))
  pred[!is.na(x)]<-fitted(fit)
  x-pred
}
raggre<-function(a,fun=namean,agg=c("NUTS0","NUTS3"),weights=NULL,debug=FALSE) {
    ag0<-attr(a,"regtype")
    agg<-match.arg(agg)
    if(debug) cat(ag0,agg,"\n")
    if(ag0==agg) return(a) # nothing NEEDS to be done
    if(ag0=="sites") {
        nreg<-as.character(with(subset(sites,usaf  %in%tregs(a)),get(paste(agg,"b",sep=""))[match(tregs(a),usaf  )]))
    }
    if(ag0=="NUTS3") {
        nreg<-as.character(with(subset(sites,NUTS3b%in%tregs(a)),get(paste(agg,"b",sep=""))[match(tregs(a),NUTS3b)]))
    }
    if(ag0=="NUTS0") {
        return(a) # nothing CAN be done
    }
    names(nreg)<-tregs(a)
    if(debug) print(head(nreg))
    if(is.null(weights)) weights<-rep(1,dim(a)[2])
    if(!is.null(names(weights))) {
        weights<-weights[match(tregs(a),names(weights))]
    } else {
        weights<-rep(weights,length=dim(a)[2])
    }
    if(debug) print(head(weights))
    out<-apply(a,c(1,3),function(x) tapply(x*weights,nreg,fun)/fun(weights))
    if(length(table(nreg))>1) {
        out<-aperm(out,c(2,1,3))
    } else {
        if(debug) print(dim(a))
        if(debug) print(str(out))
        dim(out)<-c(dim(a)[1],1,dim(a)[3])
        dimnames(out)<-list(dimnames(a)[[1]],nreg[1],dimnames(a)[[3]])
    }
    attr(out,"freq")<-tfreq(a)
    attr(out,"regtype")<-agg
    class(out)<-"tempt"
    out
}
as.data.frame.tempt<-function(a) {
    out<-as.data.frame.table(unclass(a[,,1,drop=FALSE]))[1:2]
    nam<-tname(a)
    names(out)<-c("date","region")
    for(i in nam) out[[i]]<-as.vector(a[,,i])
    out$date<-as.Date(as.character(out$date))
    attr(out,"date")<-"date"
    attr(out,"region")<-"region"
    attr(out,"vars")<-nam
    out
}
naggre<-function(a,...) {
  data<-as.data.frame(a)
  mf<-match.call(expand.dots=TRUE)[-(1:2)]
  #if(is.null(names(mf))) names(mf)<-paste("X",1:length(mf),sep="")
  for(name in names(mf))
    data[[name]]<-eval(mf[[name]],env=data)
  tempt(names(mf),"date","region",data=data)
}
