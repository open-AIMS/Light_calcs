
calculate.MS8.logger.dat=function(parsed.file.f,
                               absorption.dat=NA,
                               night.correct=T,
                               include.daily=T){
  # scratch: absorption.dat=NA;night.correct=T; include.daily=T
  require(doBy)
  require(plyr)
  require(caTools)

  #-----------------------------------------------------------------------------
  # remove cleaned data
  all.dat.f.min=parsed.file.f[which(is.na(parsed.file.f$Exclude.data)==T),]
  numeric.vars=c("year","month","day","hour","minute",
                 "nm425","nm455","nm485" ,"nm515","nm555","nm615","nm660","nm695",
                 "MS8.temperature","DL3.temperature",
                 "MS8.tilt","depth_M","voltage","NTU",
                 "ntu.tilt","NTU.temperature","deci.j.date")
  for(g in 1:length(numeric.vars)){
      all.dat.f.min[,numeric.vars[g]]=as.numeric(as.character(all.dat.f.min[,numeric.vars[g]]))
  }

  # generate some columns and objects
  all.dat.f.min$j.date=floor(all.dat.f.min$deci.j.date)
  all.dat.f.min$deci.time=all.dat.f.min$deci.j.date-all.dat.f.min$j.date
  all.dat.f.min=all.dat.f.min[order( all.dat.f.min$j.date+ all.dat.f.min$deci.time),]
  spec.colms=c("nm425","nm455","nm485","nm515","nm555","nm615","nm660","nm695")
  spec.dat=do.call("cbind",lapply(spec.colms,FUN=function(x){
        as.numeric(all.dat.f.min[,spec.colms][,x])}))
  colnames(spec.dat)=spec.colms

  #-----------------------------------------------------------------------------
  # re-set data based on min night time values
  if(night.correct==T){
   night.dat=all.dat.f.min[which(all.dat.f.min$deci.time <0.15 |
                                 all.dat.f.min$deci.time >0.85),
                                 c("j.date",spec.colms)]

   correction.f=join(all.dat.f.min[c("j.date","deci.time")],
                     summaryBy(as.formula(paste(paste(spec.colms,collapse="+"),"~j.date")),
                               data=night.dat,FUN=mean,keep.names=T),type="left")
   correction.f=correction.f[order(correction.f$j.date+correction.f$deci.time),spec.colms]
   spec.dat=(spec.dat-correction.f)
  }

  #Convert MS8EN data to W/m^2/nm
  spec.dat=spec.dat/100

  #convert to photons per second by dividing by h*c / lambda
  h=6.63e-34
  c=3.00e+08
  lambda.vec=c(425,455,485,515,555,615,660,695)
  names(lambda.vec)=spec.colms
  phot_per_sec=do.call("cbind",lapply(spec.colms,FUN=function(x){
       tt=spec.dat[,x]/(h*c/(lambda.vec[x]*1e-9))
       #tt[which(tt<0)]=0
       return(tt)}))
  colnames(phot_per_sec)=spec.colms
  #convert to microMoles per second by dividing by Avogadro's number * 1e6
  avo=6.022e17
  micromolespersec=phot_per_sec / avo
  head(micromolespersec)
  #interpolate multispectral micromolespersec to high res 1nm spacing between 400 - 700 nm.
  new.lambda.vec=400:700
  #now interpolate using your favorite function.
  micromolespersec_ipol=t(apply(micromolespersec,MARGIN=1,FUN=function(x){
         out=rep(NA,length(new.lambda.vec))
         if(length(na.omit(x))>5){ # if at least 5 nm readings
          out=approx(lambda.vec,x,xout=new.lambda.vec,rule=2)$y}
         return(out)}))
  #now calculate the integral of micromolespersec_ipol.
  all.dat.f.min$PAR.I=apply(micromolespersec_ipol,MARGIN=1,FUN=function(x){trapz(new.lambda.vec,x)})
  # Just a sum works? Implied by the excel spreadsheet
  all.dat.f.min$PAR=rowSums(micromolespersec_ipol)

  # Add PUR based on an action spectrum
  if(length(absorption.dat)==1){
   absorption.dat=read.csv("U:/1 TOWNSVILLE PORT - CLEVELAND BAY/1  REBECCA FISHER ANALYSES/RScripts/absorption_dat.csv")
  }
  absorption.dat=merge(data.frame(nm=new.lambda.vec),absorption.dat,all=T)
  absorption.vals=approx(absorption.dat$nm,
                                          absorption.dat$Absorption,
                                          xout=new.lambda.vec,rule=2)$y
  micromolespersec_ipol.pur=t(apply(micromolespersec_ipol,MARGIN=1,FUN=function(x){
        x*absorption.vals}))
  all.dat.f.min$PUR.I=apply(micromolespersec_ipol.pur,MARGIN=1,FUN=function(x){
        trapz(new.lambda.vec,x)})
  # Just a sum works? Implied by the excel spreadsheet
  all.dat.f.min$PUR=rowSums(micromolespersec_ipol.pur)

  # remove light data below zero
  all.dat.f.min$PAR.I[which(all.dat.f.min$PAR.I<0)]=0
  all.dat.f.min$PAR[which(all.dat.f.min$PAR<0)]=0
  all.dat.f.min$PUR.I[which(all.dat.f.min$PUR.I<0)]=0
  all.dat.f.min$PUR[which(all.dat.f.min$PUR<0)]=0

  ### now calcualte daily summaries, including DLI based PAR and PUR
  # first make daily dataset
  all.dat.f.daily=list()
  if(include.daily==T){
   all.dat.f.daily=summaryBy(as.formula(paste(
                                     paste(c(spec.colms,"MS8.temperature","DL3.temperature",
                                             "MS8.tilt","depth_M","voltage",
                                             "PAR.I","PAR","PUR.I","PUR",
                                             "NTU","ntu.tilt","NTU.temperature"),
                                           collapse="+"),"~",
                                     paste(c("file.f","sensor","machine_number",
                                             "Start.Date","Start.time",
                                             "j.date","Date"),
                                           collapse="+"))),
                           FUN=mean,data=all.dat.f.min,keep.names=F,na.rm=T)
  daily.max.dat.f=summaryBy(as.formula(paste(
                                     paste(c(spec.colms,"MS8.temperature","DL3.temperature",
                                             "MS8.tilt","depth_M","voltage",
                                             "PAR.I","PAR","PUR.I","PUR",
                                             "NTU","ntu.tilt","NTU.temperature"),
                                           collapse="+"),"~",
                                     paste(c("file.f","sensor","machine_number",
                                             "Start.Date","Start.time",
                                             "j.date","Date"),
                                           collapse="+"))),
                           FUN=max,data=all.dat.f.min,keep.names=F,na.rm=T)
   all.dat.f.daily=join(all.dat.f.daily,daily.max.dat.f)
   all.dat.f.daily=all.dat.f.daily[order(all.dat.f.daily$j.date),]

   # calculate DLI
   light.cols.vec=c("PAR.I","PAR","PUR.I","PUR")
   # make a vector of every second of the day
   new.dat=data.frame(deci.time=seq(from=0,to=1,by=((1/60)/60)/24))

   dli.out=lapply(all.dat.f.daily$j.date,FUN=function(x){
     x.dat=na.omit(all.dat.f.min[which(all.dat.f.min$j.date==x),
                          c("deci.time",light.cols.vec)])
     DLI.PAR.I=NA; DLI.PAR=NA; DLI.PUR.I=NA; DLI.PUR=NA
     if(nrow(x.dat)>6*4){
      DLI.PAR.I=sum(approx(x.dat$deci.time,x.dat$PAR.I,xout=new.dat$deci.time,rule=2)$y)/(10^6)
      DLI.PAR=sum(approx(x.dat$deci.time,x.dat$PAR,xout=new.dat$deci.time,rule=2)$y)/(10^6)
      DLI.PUR.I=sum(approx(x.dat$deci.time,x.dat$PUR.I,xout=new.dat$deci.time,rule=2)$y)/(10^6)
      DLI.PUR=sum(approx(x.dat$deci.time,x.dat$PUR,xout=new.dat$deci.time,rule=2)$y)/(10^6)
     }
     out=c(DLI.PAR.I,DLI.PAR,DLI.PUR.I,DLI.PUR)
     names(out)=c("DLI.PAR.I","DLI.PAR","DLI.PUR.I","DLI.PUR")
     return(out)})
   all.dat.f.daily=cbind(all.dat.f.daily,do.call("rbind",dli.out))
  }

### Now return all of the datasets ---------------------------------------------
  return(list(
    all.dat.minute=all.dat.f.min,
    all.dat.daily=all.dat.f.daily)

)
}
### end.function ---------------------------------------------------------------











