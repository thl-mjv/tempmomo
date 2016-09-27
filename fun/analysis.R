#' Smart calculation of lags
#'
#' Create a matrix of lagged values. Assumes equidistant, presorted values
#'
#' @param a original variable, could be a matrix
#' @param n which lags to use,
#' @param na value to use for the values 
momolag<-function(a,n=1,na=NA) {
  ca<-data.class(a)
  if(is.null(dim(a))) {
    a<-matrix(a,ncol=1)
    dimnames(a)<-list(a,"")
  }
  if(length(dim(a))==1) {
    dim(a)<-c(dim(a),1)
    dimnames(a)<-list(dimnames(a),"")
  }
  nc<-ncol(a)
  an<-nrow(a)
  nn<-dimnames(a)[[2]]
  if(is.null(nn)) nn<-rep("",nc)
  res<-lapply(n,function(n) {
    if(n==0) {
      return(a)
    }
    nas<-matrix(na,nrow=pmin(abs(n),an),ncol=nc)
    if(an<n)
      return(nas)
    dro<-1:pmin(abs(n),an) 
    if(n>0)
      return(rbind(nas,a[dro-an-1,,drop=FALSE]))
    else
      return(rbind(a[-dro,,drop=FALSE],nas))
  })
  res<-do.call("cbind",res)
  ##  dim(res)<-c(na,nc*length(n))
  dimnames(res)<-list(dimnames(a)[[1]],as.vector(t(outer(nn,n,paste,sep=""))))
  if(ca=="Date") {
    at<-attributes(res)
    res<-as.Date(res,origin="1970-1-1")
    attributes(res)<-c(at,class=class(res))
  }
  res
}
#' Sine terms
#'
#' @param a    original variable
#' @param n    number of terms peer year
#' @param freq length of single cycle assuming it can not be determined from the data
#' @param use  which terms to use (default 1:n)
momosin<-function(a,n=1,freq=NA,use=NULL) UseMethod("momosin")
momosin.time<-function(a,n=1,freq=NA,use=NULL) {
  if(is.na(freq)) freq<-frequency(a)
  momosin.default(unclass(a),n=n,freq=freq,use=use)
}
momosin.Date<-function(a,n=1,freq=NA,use=NULL) {
  if(is.na(freq)) freq<-1 
  x<-as.POSIXlt(a)
  momosin.default(x$year+x$yday/365.25,n=n,freq=freq,use=use) # slightly wrong for leap years
}
momosin.character<-function(a,n=1,freq=NA,use=NULL,pad=0) {
  if(is.na(freq)) freq<-frequency(a)
  ## assume that this is of the form YYYYxWW
  yr<-as.numeric(substring(a,1,4))
  wk<-as.numeric(substring(a,6+pad,7+pad))
  momosin.default(yr+wk/52.25,n=n,freq=1,use=use)
}
momosin.default<-function(a,n=1,freq=NA,use=NULL) {
  if(is.na(freq))
    stop("Could not determine frequency")
  if(is.null(use))
    use<-1:n
  n<-length(use)
  res<-matrix(NA,nr=length(a),nc=2*n)
  dimnames(res)<-list(NULL,paste(rep(c("Sin","Cos"),c(n,n)),rep(use,2)))
  for(i in use) res[,paste("Sin",i)]<-sin(2*pi*i*a/(freq))
  for(i in use) res[,paste("Cos",i)]<-cos(2*pi*i*a/(freq))
  res
}
#' Season
#'
#' Influenza seasons
#'
#' @param time A time variable
#' @param start Which week starts the season
momoseason<-function(time,start=26) UseMethod("momoseason")
momoseason.Date<-function(time,start=26) {
  yy<-as.numeric(format(time,"%G")) ## ISO8601
  ww<-as.numeric(format(time,"%V"))
  ifelse(ww<start, yy-1,yy)
}
momoseason.character<-function(time,start=26,pad=0) {
  yy<-as.numeric(substring(a,1,4))
  ww<-as.numeric(substring(a,6+pad,7+pad))
  ifelse(ww<start, yy-1,yy)
}
### output as factor for regression and tables
momoseasonf<-function(a) { 
  res<-momoseason(a)
  lev<-sort(unique(res))
  lnam<-paste(lev,lev+1,sep="-")
  out<-factor(res,levels=lev,labels=lnam)
  co<-try(contr.sum(length(lev)))
  if(!inherits(co,"try-error")) {
    dimnames(co)<-list(lnam,lnam[-length(lnam)])
    contrasts(out,)<-co
  }
  out
}
#' Prediction
#'
#' Something witty about prediction (this is a very crude version but can be used as a placeholder for better implementation)
#'
#' @param fit model fit to use
#' @param newdata new data to use
predfun<-function(fit,newdata) {
    res<-newdata
    pred<-predict(fit,newdata=newdata,type="lin",se=TRUE)
    ## fixme: prediction intervals
    ## fixme: proper transformations
    res$pred  <-exp(with(pred,fit))
    res$predup<-exp(with(pred,fit+2*se.fit)) # CRUDE
    res$predlo<-exp(with(pred,fit-2*se.fit))
    class(res)<-c("pred.frame","data.frame")
    res
}
#' Plot prediction
#'
#' Plots predictions
#'
#' @param pred  pred.frame object from predfun
#' @param var   covariate to use as X
#' @param col   color of prediction
#' @param alpha transparency of the confidence area
#' @param add   if TRUE overlay on previous plot
#' @param ...   other options passed to lines
plot.pred.frame<-function(pred,var,col="black",alpha=.2,add=TRUE,...) {
    col<-col2rgb(col)
    col1<-rgb(col[1],col[2],col[3],255*1    ,maxColorValue=255) # full intensity
    col2<-rgb(col[1],col[2],col[3],255*alpha,maxColorValue=255)
    with(pred,{
        x<-get(var)
        if(!add) plot (x,pred,type="n")
        polygon(c(x,rev(x)),c(predlo,rev(predup)),col=col2,border=NA)
        lines(x,pred,col=col1,...)
    })
    invisible(list(col1,col2))
}

#' Create a linear spline basis
#'
#' For modeling temperature, create a spline basis that estimates the effect of departures from some preset reference value
#'
#' @param a   the variable for the basis
#' @param by  distance between knots
#' @param n   number of knots ot each direction
#' @param ref the reference value
wlsb<-function(a,by=1,n=20,ref=0) {
    out<-matrix(NA,nrow=length(a),ncol=2*n)
    na<-as.character(a)
    a<-a-ref
    for(i in 1:n-1) {
        out[,n-i  ]<-pmax(0,-a+(-i*by)) # left
        out[,n+i+1]<-pmax(0, a-( i*by)) # right
    }
    dimnames(out)<-list(na,c(paste("left",n:1-1,sep=""),paste("right",1:n-1,sep="")))
    out
}
#' Create a general spline basis
#'
#' For modeling temperature, create a spline basis based on same knots as wlsb
#'
#' @param a   the variable for the basis
#' @param by  distance between knots
#' @param n   number of knots ot each direction
#' @param ref the reference value
#' @param degree degree of the spline
wlsbs<-function(a,by=1,n=20,ref=0,degree=1) {
    knots<-(-n:n)*by-ref
    out<-bs(a,knots=knots,degree=degree,Boundary.knots=range(knots))
    dimnames(out)<-list(a,paste("b",1:ncol(out),sep=""))
    out
}
#' Create indicators for excess
#'
#' For modeling temperature, create variables giving numbers of excess weeks
#'
#' @param a   the variable for the basis
#' @param by  distance between knots
#' @param n   number of knots ot each direction
#' @param ref the reference value
wsei<-function(a,by=1,n=20,ref=0) {
    out<-matrix(NA,nrow=length(a),ncol=2*n)
    for(i in 1:n-1) {
        out[,n-i  ]<-1*(a< (-i*by)) # left
        out[,n+i+1]<-1*(a> ( i*by)) # right
    }
    dimnames(out)<-list(a,c(paste("left",n:1-1,sep=""),paste("right",1:n-1,sep="")))
    out
}
#' Create dummy variables
#'
#' Helper function for creating a matrix of 0,1 dummies for each category. 
#'
#' @param a    A grouping variable such as factor
#' @param nlev fixed size of the resulting matrix
dodummy<-function(a,nlev=NULL) {
    le<-names(table(a)) # missing values omitted, also non existent levels included
    mi<-"dodummy missing value"
    if(!is.null(nlev)) le<-c(le,rep(NA,nlev))[1:nlev]
    res<-1*outer(na.0(as.character(a),mi),le,"==")
    dimnames(res)<-list(a,le)
    res
}
#' Replacement strings for cleanmat
cleanmat<-list(c("momoseasonf","Season"),c("mktrend","Trend"),c("date",""),c("temp","TEMP"),
               c("momosin","Sin"),c("tmpdate","Trend"), # can't have both tmpdate and mktrend
               c(".all",""),c("mkseas",""),c("wseas",""),c(".HA",""),c(".adj2",""),c(".adj",""),c("abs",""),c("pop",""),c("bsdeg",""),
               c("bsthin",""),c("momolag",""),c("pm",""),c("[():/<>, =*0-9-]",""),
               c("c[(]",""),c("s[(]",""),c("bs[(]",""),c("I[(]",""))

#' Clean effect names from a formula
#'
#' @param a character vector from formula object
cleannames<-function(a) {
  for(i in rev(cleanmat)) {
    a<-gsub(i[1],i[2],a)
  }  
  a
}
#' Model matrix
#'
#' Get the disegn matrix from a model object for a different data
#'
#' @param mod   model object
#' @param data  new data
modmat<-function(mod,data) {
    UseMethod("modmat")
}
modmat.default<-function(mod,data) {
    cat("modmat.default\n")
    model.matrix(formula(mod),data=model.frame(formula(mod),data=data))#,na.action=na.fail))
}
modmat.gam<-function(mod,data) {
    cat("modmat.gam\n")
    model.matrix(mod,newdata=data)
}
#' Model effects
#'
#' Calculate model effects (predictions, excesses attributable to a factor) from a fit object
#'
#' @param mod   model fit object from glm, gam, gamlss, ...
#' @param data  data used in prediction
#' @param debug should debugging infortmation be printed during the calculations?
#' @param response a string giving the name of the response variable
#' @param date     a string giving the name of the date variabe for defining the seasons
#' @param savecoc  should the variance covariance matrices for the effects be saved in outputs (these are large so avoid if possible)?
#' @param outeffs  names of the effects to be output
modeffects<-function(mod,data,debug=FALSE,response="n",date="date",savecov=FALSE,clean=FALSE,
                     outeffs=c("Intercept","Season","Trend","Sin","TEMP","SLP","DEWP","InfA","InfB","RSV")) {
  require("survival")
  ## decide the link function
  modfamily<-try(family(mod))
  invlink<-NULL
  if(inherits(modfamily,"character")) {invlink<-get(modfamily)()$mu.linkinv;} ###?
  if(inherits(modfamily,"family"   )) {invlink<-modfamily$linkinv;varfun<-modfamily$mu.eta}
  if(inherits(modfamily,"try-error")) stop("Not a model object")
  if(is.null(invlink)) stop("Unknown family")
  ## linear prediction
  ### not working with GAM. Need to do a function for model martices.
  mm<-modmat(mod,data=data)#,na.action=na.fail))
  aa<-attrassign(mod) # from survival
  ## clean names
  if(debug) print(names(aa))
  if(clean)
      names(aa)<-cleannames(names(aa))
  if(debug) print((aa))
  ## OFFSET HERE!
  ##
  ## Also: chect that all is in outeffs 
  aa<-aa[match(outeffs,names(aa))]
  names(aa)<-outeffs
  coe<-na.0(coef(mod))
  vaa<-vcov(mod)
  vaa<-na.0(vaa[match(names(coe),dimnames(vaa)[[1]]),match(names(coe),dimnames(vaa)[[2]])])            
  bb<-NULL
  lps<-matrix(NA,nrow=nrow(mm),ncol=length(aa))
  if(debug) print(dim(lps))
  dimnames(lps)<-list(data[[date]],names(aa))
  effs<-list(add=lps,res=lps)
  ceff<-array(NA,dim=list(nrow(mm),nrow(mm),length(aa)),dimnames=list(data[[date]],data[[date]],names(aa)))
  covs<-list(add=ceff,res=ceff)
  eff0<-na.0(lps[,1]) # starting values for effects
  vf0<-0*mm
  for(i in names(aa)) {
      cat(i,"...")
      bb<-c(bb,aa[[i]]) # cumulatively
      #if(debug) print(bb)
      cat(length(bb),"...")
      lm<-mm
      lm[,-bb]<-0 # set the remaining effects at 0
      #if(debug) print(head(lm))
      cat("lp...")
      lp<-as.vector(lm%*%coe)
      cat("eff...")
      eff<-as.vector(na.0(invlink(lp)))
      cat("vf...")
      vf<-varfun(lp)*lm
      cat("eff$add..")
      effs$add[,i]<-eff-eff0
      cat("eff$res..")
      effs$res[,i]<-data[[response]]-eff
      cat("cov$add..")
      covs$add[,,i]<-(vf-vf0)%*%vaa%*%t(vf-vf0)
      cat("cov$res..")
      covs$res[,,i]<-diag(eff)+vf%*%vaa%*%t(vf)
      eff0<-eff
      vf0<-vf
      cat(summary(eff),"\n")

  }
  ses<-lapply(covs,apply,3,function(a) sqrt(diag(a)))
  dummies<-list()
  ### seasonal
  dummies$season<-dodummy(momoseason(data[[date]]))
  ### winter
  dummies$winter<-dodummy(ifelse( as.numeric(format(data[[date]],"%V"))%in%c(40:53,1:20),momoseason(data[[date]]),NA),nlev=ncol(dummies$season))
  ### summer
  dummies$summer<-dodummy(ifelse(!as.numeric(format(data[[date]],"%V"))%in%c(40:53,1:20),format(data[[date]],"%Y"),NA),nlev=ncol(dummies$season))
  ###
  outmat<-array(NA,dim=c(length(effs),length(aa),length(dummies),ncol(dummies[[1]]),3),
                dimnames=list(names(effs),names(aa),names(dummies),dimnames(dummies[[1]])[[2]],c("Est","(",")")))
  if(debug) print(sapply(dummies,dim))
  nnames<-function(a) structure(names(a),names=names(a))
  for(e in names(effs)) {
      for(a in names(aa)) {
          for(b in names(dummies)) {
              cat(e,a,b,"... ")
              oeff<-as.vector(t(dummies[[b]])%*%effs[[e]][,a])
              oses<-sqrt(diag(t(dummies[[b]])%*%covs[[e]][,,a]%*%dummies[[b]]))
              outmat[e,a,b,,]<-oeff+cbind(0,-2*oses,2*oses)
              cat(round(oeff),"\n")
          }
      }
  }
  if(!savecov) covs<-NULL
  list(effects=effs,
       covariances=covs,
       ses=ses,
       seaseffs=outmat,
       resp=data[[response]],date=data[[date]])
}
#' Create data for contour predictions
#'
#' @param weeks dates to use
#' @param data  data to use as background
#' @param basedate date of the start of the year
mkpreddata<-function(weeks=1:366,data,basedate="2010-1-1") {
    preddata<-expand.grid(seas=weeks/366,temp=-30:30)
    preddata$date  <-as.Date(basedate)+preddata$seas*366
    preddata$time  <-numdate(preddata$date)
    preddata$month <-as.numeric(format(preddata$date,"%m"))
    preddata$monthf<-factor(preddata$month)
    try(preddata$temp.pred1<-predict(gam(temp~s(seas,bs="cc"),data=data),newdata=preddata))
    try(preddata$temp.pred2<-predict(gam(temp~momosin(date,4),data=data),newdata=preddata))
    try(preddata$temp.adj3 <-with(preddata,temp-temp.pred1))
    try(preddata$temp.adj1 <-with(preddata,temp-temp.pred2))
    preddata
}
#' Plot contour predictions
#' 
#' @param preds object with predictions
#' @param country country to plot
#' @param model model to plot
plotpredcont<-function(pred,country,model) {
  if(!is.null(preds[[country]][[model]])) {
    contour(weeks,-30:30,preds[[country]][[model]],zlim=range(morts[[country]]$n),
            nlevels=100,main=paste(country),sub=model,ylim=range(morts[[i]]$temp),
            col=colorRampPalette(c("blue","red"))(100),drawlabels=FALSE)
    with(morts[[country]],symbols(seas*366,temp,circles=sqrt(n/400),add=TRUE,inches=FALSE,
                                  fg=NA,bg=rgb(0,0,0,.75)))
  } else {
    frame()
  }
  invisible(NULL)  
}
plotpredpair<-function(country,model) {
  with(morts[[country]],{
    pred<-get(paste("pred",model,sep=""))
    par(mar=c(3,3,2,0))
    plot(temp,pred,main=i,pch=16,ylim=range(c(pred,n)))
    matlines(rbind(temp,temp),rbind(pred,n),type="l",lty=1,col=2)
    par(mar=c(3,0,2,3))
    plot(seas,pred,main=i,pch=16,ylim=range(c(pred,n)),yaxt="n")
    axis(4)
    matlines(rbind(seas,seas),rbind(pred,n),type="l",lty=1,col=2)
  })
}