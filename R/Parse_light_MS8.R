
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

#' parse.MS8.logger.dat
#'
#' Reads in MS8 logger data
#'
#' @param  file.f The MS8 logger file to be processed (including path, if required)
#'
#' @export
#' @return A list containing the processed MS8 logger data, including a combined dataset as well as the individual channels present.
#' 
parse.MS8.logger.dat=function(file.f){
  # scratch: absorption.dat=NA;night.correct=T; include.daily=T
  require(doBy)
  #require(plyr)
  #require(caTools)

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
  MS8.dat.f.out$MS8.tilt=as.numeric(MS8.dat.f[,"tilt"])
  MS8.dat.f.out$MS8.temperature=as.numeric(gsub(">","",MS8.dat.f[,"Temp"]))

  # extract the light sensor data for each wavelength
  spec.colms=c("nm425","nm455","nm485","nm515","nm555","nm615","nm660","nm695")
  spec.dat=do.call("cbind",lapply(spec.colms,FUN=function(x){
        as.numeric(MS8.dat.f[,x])}))
  colnames(spec.dat)=spec.colms
  MS8.dat.f.out=cbind(file.f,MS8.dat.f.out,spec.dat)

  ### now merge the depth from the DL3 and NTU sensors  -------------------------
  MS8.key=paste(MS8.dat.f.out$j.date,MS8.dat.f.out$hour,MS8.dat.f.out$minute,round(MS8.dat.f.out$second))
  MS8.dat.f.out$key=MS8.key
  DL3.key=paste(DL3.dat.f.out$j.date,DL3.dat.f.out$hour,DL3.dat.f.out$minute,round(DL3.dat.f.out$second))
  DL3.dat.f.out$key=DL3.key

  # now merge
  all.dat.f=merge(MS8.dat.f.out,DL3.dat.f.out[,c("key","depth_M","DL3.temperature","voltage")],all.x=T)

  if(length(NTU.dat.f)==0){
   all.dat.f$NTU=NA
   all.dat.f$ntu.tilt=NA
   all.dat.f$NTU.temperature=NA}else{
   NTU.key=paste(NTU.dat.f.out$j.date,NTU.dat.f.out$hour,NTU.dat.f.out$minute,round(NTU.dat.f.out$second))
   NTU.dat.f.out$key=NTU.key
   all.dat.f=merge(all.dat.f,NTU.dat.f.out[,c("key","NTU","ntu.tilt","NTU.temperature")],all.x=T)
  }

  all.dat.f$day=format(all.dat.f$Date, format="%d")
  all.dat.f$month=format(all.dat.f$Date, format="%m")
  all.dat.f$year=format(all.dat.f$Date, format="%Y")
  all.dat.f.min=summaryBy(as.formula(paste(
                                     paste(c(spec.colms,"MS8.temperature","DL3.temperature",
                                             "MS8.tilt","depth_M","voltage",
                                             "NTU","ntu.tilt","NTU.temperature","deci.j.date"),
                                           collapse="+"),"~",
                                     paste(c("file.f","port","sensor","machine_number","Date",
                                           "year","month","day","hour","minute"),
                                           collapse="+"))),
                           FUN=mean,data=all.dat.f,keep.names=T,na.rm=T)
  all.dat.f.min=all.dat.f.min[order(all.dat.f.min$deci.j.date),]


all.dat.f.min$Exclude.data=NA
all.dat.f.min$Notes=NA

### Now return all of the datasets ---------------------------------------------
  return(list(
   raw.dat=list(MS8=MS8.dat.f.out,
                DL3=DL3.dat.f.out,
                NTU=NTU.dat.f.out),
   all.dat.merged=all.dat.f.min)
)
}
### end.function ---------------------------------------------------------------












