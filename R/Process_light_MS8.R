#    Copyright 2020 Australian Institute of Marine Science
#
#    Licenced under a Creative Commons Attribution 4.0 International licence.
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       hhttps://creativecommons.org/licenses/by/4.0/
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

#' process.MS8.logger.dat
#'
#' Calculates PAR and PUR for parsed MS8 logger data
#'
#' @param  file.f The parsed MS8 logger file to be processed (including path, if required)
#'
#' @param absorption.dat A data.frame containing an absorption spectrum to use for calculating PUR. Must contain columns fo nm and normalized absoption/action (values between zero and one)
#' 
#' @param night.correct Should a night correction be applied to the calculations?
#' 
#' @param include.daily Should daily data be included in the output?
#'
#' @export
#' @return A list containing the processed MS8 logger data.

proccess.MS8.logger.dat=function(file.f,
                                 absorption.dat,
                                 night.correct=T,
                                 include.daily=T){
  # scratch: absorption.dat=NA;night.correct=T; include.daily=T
  require(doBy)
  require(plyr)
  require(caTools)
  
  if(missing(absorption.dat)){
    absorption.dat  <-  Light.calcs:::data_list$absorptionDat
  }

  tt=scan(file.f,what="character")
  tt=unlist(strsplit(tt,split=","))
  xx=unlist(lapply(tt,FUN=nchar))
  tt=tt[which(xx>0)]

  log.dat.f=unlist(strsplit(c(paste(tt[1:3],collapse=" "),paste(tt[4:5],collapse=" "),
                    tt[7],tt[11:12],tt[14]),split=","))
  names(log.dat.f)=c("source","device","firmware","session_start_date",
                     "session_start_time","log_file")

  port.dat.f=list()
  port.indices.f=grep("Port",tt)
  for(p in 1:(length(port.indices.f)-1)){
       port.dat.f=c(port.dat.f,list(tt[port.indices.f[p]:((port.indices.f[p+1])-1)]))
  }
  port.dat.f=c(port.dat.f,list(tt[port.indices.f[length(port.indices.f)]:length(tt)]))

  #DL3 data
  DL3.index=unlist(lapply(port.dat.f,FUN=function(x){length(grep("DL3",x))>0}))
  DL3.dat.raw=port.dat.f[DL3.index]
  DL3.length=table(unlist(lapply(DL3.dat.raw,FUN=length)))
  dominant.length=as.numeric(names(DL3.length)[which.max(DL3.length)])
  DL3.dat.f=do.call("rbind",DL3.dat.raw[unlist(lapply(DL3.dat.raw,FUN=function(x){length(x)==dominant.length}))])
  if(length(DL3.dat.f)>0){
  colnames(DL3.dat.f)=c("port","date","time","sensor","device_number","voltage_raw_counts",
                        "depth_raw_counts","temp_raw_counts","voltage","depth_M","Temp")}
  head(DL3.dat.f)

  # MS8 data
  MS8.index=unlist(lapply(port.dat.f,FUN=function(x){length(grep("MS8",x))>0}))
  MS8.dat.raw=port.dat.f[MS8.index]
  MS8.length=table(unlist(lapply(MS8.dat.raw,FUN=length)))
  dominant.length=as.numeric(names(MS8.length)[which.max(MS8.length)])
  MS8.dat.f=do.call("rbind",MS8.dat.raw[unlist(lapply(MS8.dat.raw,FUN=function(x){length(x)==dominant.length}))])
  if(length(MS8.dat.f)>0){
  colnames(MS8.dat.f)=c("port","Date_1",
   "time_1",
   "sensor",
   "machine_number",
   "Date",
   "time",
   "nm425", "nm455", "nm485", "nm515", "nm555", "nm615", "nm660", "nm695",
   "tilt","Temp")}
  head(MS8.dat.f)

  # NTU data if available
  NTU.index=unlist(lapply(port.dat.f,FUN=function(x){length(grep("NTU",x))>0}))
  NTU.dat.raw=port.dat.f[NTU.index]
  NTU.length=table(unlist(lapply(NTU.dat.raw,FUN=length)))
  dominant.length=as.numeric(names(NTU.length)[which.max(NTU.length)])
  NTU.dat.f=do.call("rbind",NTU.dat.raw[unlist(lapply(NTU.dat.raw,FUN=function(x){length(x)==dominant.length}))])
  if(length(NTU.dat.f)>0){
  colnames(NTU.dat.f)=c("port","Date_1","time_1","sensor","machine_number",
                        "","","dark.counts","meas.counts","NTU","ntu.tilt","ntu.temp")}
  head(NTU.dat.f)

  ### process DL3 data -----------------------------------------------------------
  head(DL3.dat.f)
  DL3.dat.f.out=data.frame(DL3.dat.f[,c("port","sensor","device_number")])
  DL3.dat.f.out$Date=as.Date(DL3.dat.f[,"date"],format="%d/%m/%Y")
  DL3.dat.f.out$Time=DL3.dat.f[,"time"]
  Start.Date=as.Date(log.dat.f["session_start_date"],format="%d/%m/%Y")
  DL3.dat.f.out$Start.Date=Start.Date
  DL3.dat.f.out$Start.time=log.dat.f["session_start_time"]
  DL3.dat.f.out$j.date=julian(DL3.dat.f.out$Date,origin=Start.Date)
  ff=strsplit(DL3.dat.f[,"time"],split=":")
  length.ff=unlist(lapply(ff,FUN=length))
  ff[which(length.ff!=3)]=rep(NA,3)
  time.matrix=do.call("rbind",ff)
  DL3.dat.f.out$hour=as.numeric(time.matrix[,1])
  DL3.dat.f.out$minute=as.numeric(time.matrix[,2])
  DL3.dat.f.out$second=as.numeric(time.matrix[,3])
  DL3.dat.f.out$deci.time=((DL3.dat.f.out$hour*60*60)+
                           (DL3.dat.f.out$minute*60)+
                            DL3.dat.f.out$second)/(24*60*60)
  DL3.dat.f.out$deci.j.date=DL3.dat.f.out$j.date + DL3.dat.f.out$deci.time
  vec.rr=c("voltage_raw_counts","depth_raw_counts","temp_raw_counts","voltage","depth_M")
  dl3.rr=do.call("cbind",lapply(vec.rr,FUN=function(x){as.numeric(DL3.dat.f[,x])}))
  colnames(dl3.rr)=vec.rr
  DL3.dat.f.out=cbind(DL3.dat.f.out,dl3.rr)
  DL3.dat.f.out$DL3.temperature=as.numeric(gsub(">","",DL3.dat.f[,"Temp"]))

  ### process NTU data -----------------------------------------------------------
  head(NTU.dat.f)
  NTU.dat.f.out=list()
  if(length(NTU.dat.f)>0){
  NTU.dat.f.out=data.frame(NTU.dat.f[,c("port","sensor","machine_number")])
  NTU.dat.f.out$Date=as.Date(NTU.dat.f[,"Date_1"],format="%d/%m/%Y")
  NTU.dat.f.out$Time=NTU.dat.f[,"time_1"]
  Start.Date=as.Date(log.dat.f["session_start_date"],format="%d/%m/%Y")
  NTU.dat.f.out$Start.Date=Start.Date
  NTU.dat.f.out$Start.time=log.dat.f["session_start_time"]
  NTU.dat.f.out$j.date=julian(NTU.dat.f.out$Date,origin=Start.Date)
  ff=strsplit(NTU.dat.f[,"time_1"],split=":")
  length.ff=unlist(lapply(ff,FUN=length))
  table(length.ff)
  ff[which(length.ff!=3)]=rep(NA,3)
  time.matrix=do.call("rbind",ff)
  NTU.dat.f.out$hour=as.numeric(time.matrix[,1])
  NTU.dat.f.out$minute=as.numeric(time.matrix[,2])
  NTU.dat.f.out$second=as.numeric(time.matrix[,3])
  NTU.dat.f.out$deci.time=((NTU.dat.f.out$hour*60*60)+
                           (NTU.dat.f.out$minute*60)+
                            NTU.dat.f.out$second)/(24*60*60)
  NTU.dat.f.out$deci.j.date=NTU.dat.f.out$j.date + NTU.dat.f.out$deci.time
  vec.rr=c("dark.counts","meas.counts","NTU","ntu.tilt")
  ntu.rr=do.call("cbind",lapply(vec.rr,FUN=function(x){as.numeric(NTU.dat.f[,x])}))
  colnames(ntu.rr)=vec.rr
  NTU.dat.f.out=cbind(NTU.dat.f.out,ntu.rr)
  #NTU.dat.f.out$NTU.temperature=as.numeric(gsub(">","",NTU.dat.f[,"Temp"]))
  NTU.dat.f.out$NTU.temperature=as.numeric(do.call("rbind",
            strsplit(NTU.dat.f[,"ntu.temp"],split="*",fixed=T))[,1])
 }

  ### process MS8 data -----------------------------------------------------------
  head(MS8.dat.f)
  MS8.dat.f.out=data.frame(MS8.dat.f[,c("port","sensor","machine_number")])
  MS8.dat.f.out$Date=as.Date(MS8.dat.f[,"Date_1"],format="%d/%m/%Y")
  MS8.dat.f.out$Time=MS8.dat.f[,"time_1"]
  Start.Date=as.Date(log.dat.f["session_start_date"],format="%d/%m/%Y")
  MS8.dat.f.out$Start.Date=Start.Date
  MS8.dat.f.out$Start.time=log.dat.f["session_start_time"]
  MS8.dat.f.out$j.date=julian(MS8.dat.f.out$Date,origin=Start.Date)
  ff=strsplit(MS8.dat.f[,"time_1"],split=":")
  length.ff=unlist(lapply(ff,FUN=length))
  ff[which(length.ff!=3)]=rep(NA,3)
  time.matrix=do.call("rbind",ff)
  MS8.dat.f.out$hour=as.numeric(time.matrix[,1])
  MS8.dat.f.out$minute=as.numeric(time.matrix[,2])
  MS8.dat.f.out$second=as.numeric(time.matrix[,3])
  MS8.dat.f.out$deci.time=((MS8.dat.f.out$hour*60*60)+
                           (MS8.dat.f.out$minute*60)+
                            MS8.dat.f.out$second)/(24*60*60)
  MS8.dat.f.out$deci.j.date=MS8.dat.f.out$j.date + MS8.dat.f.out$deci.time

  # extract the light sensor data for each wavelength
  spec.colms=c("nm425","nm455","nm485","nm515","nm555","nm615","nm660","nm695")
  spec.dat=do.call("cbind",lapply(spec.colms,FUN=function(x){
        as.numeric(MS8.dat.f[,x])}))
  colnames(spec.dat)=spec.colms
  MS8.dat.f.out=cbind(file.f,MS8.dat.f.out,spec.dat)

  # re-set data based on min night time values
  if(night.correct==T){
  night.dat=MS8.dat.f.out[which(MS8.dat.f.out$deci.time <0.15 |
                                MS8.dat.f.out$deci.time >0.85),
                                c("j.date",spec.colms)]

  correction.f=join(MS8.dat.f.out[c("j.date","deci.time")],
                     summaryBy(as.formula(paste(paste(spec.colms,collapse="+"),"~j.date")),
                               data=night.dat,FUN=mean,keep.names=T),type="left")[,spec.colms]

  spec.dat=(spec.dat-correction.f)}

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
  MS8.dat.f.out$PAR.I=apply(micromolespersec_ipol,MARGIN=1,FUN=function(x){trapz(new.lambda.vec,x)})
  # Just a sum works? Implied by the excel spreadsheet
  MS8.dat.f.out$PAR=rowSums(micromolespersec_ipol)
  MS8.dat.f.out$tilt=as.numeric(MS8.dat.f[,"tilt"])
  MS8.dat.f.out$MS8.temperature=as.numeric(gsub(">","",MS8.dat.f[,"Temp"]))

  # Add PUR based on an action spectrum
  absorption.dat=merge(data.frame(nm=new.lambda.vec),absorption.dat,all=T)
  absorption.vals=approx(absorption.dat$nm,
                                          absorption.dat$Absorption,
                                          xout=new.lambda.vec,rule=2)$y
  micromolespersec_ipol.pur=t(apply(micromolespersec_ipol,MARGIN=1,FUN=function(x){
        x*absorption.vals}))
  MS8.dat.f.out$PUR.I=apply(micromolespersec_ipol.pur,MARGIN=1,FUN=function(x){
        trapz(new.lambda.vec,x)})
  # Just a sum works? Implied by the excel spreadsheet
  MS8.dat.f.out$PUR=rowSums(micromolespersec_ipol.pur)

  ### now merge the depth from the DL3 and NTU sensors  -------------------------
  MS8.key=paste(MS8.dat.f.out$j.date,MS8.dat.f.out$hour,MS8.dat.f.out$minute,round(MS8.dat.f.out$second))
  MS8.dat.f.out$key=MS8.key
  DL3.key=paste(DL3.dat.f.out$j.date,DL3.dat.f.out$hour,DL3.dat.f.out$minute,round(DL3.dat.f.out$second))
  DL3.dat.f.out$key=DL3.key
  # now merge
  all.dat.f=merge(MS8.dat.f.out,DL3.dat.f.out[,c("key","depth_M","DL3.temperature","voltage")],all.x=T)
  all.dat.f$NTU=NA
  all.dat.f$ntu.tilt=NA
  all.dat.f$NTU.temperature=NA

  if(length(NTU.dat.f)>0){
  NTU.key=paste(NTU.dat.f.out$j.date,NTU.dat.f.out$hour,NTU.dat.f.out$minute,round(NTU.dat.f.out$second))
  NTU.dat.f.out$key=NTU.key
  all.dat.f=merge(all.dat.f,NTU.dat.f.out[,c("key","NTU","ntu.tilt","NTU.temperature")],all.x=T)
  }


  all.dat.f.min=summaryBy(as.formula(paste(
                                     paste(c("second","deci.time","deci.j.date",
                                             spec.colms,"MS8.temperature","DL3.temperature",
                                             "tilt","depth_M","voltage",
                                             "PAR.I","PAR","PUR.I","PUR",
                                             "NTU","ntu.tilt","NTU.temperature"),
                                           collapse="+"),"~",
                                     paste(c("file.f","port","sensor","machine_number","Start.Date","Start.time",
                                             "j.date","hour","minute" ,"Date"),
                                           collapse="+"))),
                           FUN=mean,data=all.dat.f,keep.names=T,na.rm=T)
  all.dat.f.min=all.dat.f.min[order(all.dat.f.min$deci.j.date),]
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
                                             "tilt","depth_M","voltage",
                                             "PAR.I","PAR","PUR.I","PUR",
                                             "NTU","ntu.tilt","NTU.temperature"),
                                           collapse="+"),"~",
                                     paste(c("file.f","sensor","machine_number",
                                             "Start.Date","Start.time",
                                             "j.date","Date"),
                                           collapse="+"))),
                           FUN=mean,data=all.dat.f,keep.names=F,na.rm=T)
  daily.max.dat.f=summaryBy(as.formula(paste(
                                     paste(c(spec.colms,"MS8.temperature","DL3.temperature",
                                             "tilt","depth_M","voltage",
                                             "PAR.I","PAR","PUR.I","PUR",
                                             "NTU","ntu.tilt","NTU.temperature"),
                                           collapse="+"),"~",
                                     paste(c("file.f","sensor","machine_number",
                                             "Start.Date","Start.time",
                                             "j.date","Date"),
                                           collapse="+"))),
                           FUN=max,data=all.dat.f,keep.names=F,na.rm=T)
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
     if(nrow(x.dat)>20){
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
   raw.dat=list(MS8=MS8.dat.f.out,
                DL3=DL3.dat.f.out,
                NTU=NTU.dat.f.out),
   all.dat.minute=all.dat.f.min,
   all.dat.daily=all.dat.f.daily)
)
}
### end.function ---------------------------------------------------------------












