### Demonstration of light manipulation scripts to run the light model
# developed by IMSO and to calculate PAR and PUR from MS8 logger files.
# written: 31-aug-2018
# last edited:
# author: r.fisher@aims.gov.au
# with contributions from Weijech, Simon Spagnol and Ross Jones

# requires 
require(doBy) 
require(caTools)
require(plyr)
require(NISTunits)


# NOTE 2: The Process_light_MS8.R function does not do any data cleaning, aside from a night
# based linear offset shift (assumption is that night readings = zero light), and removal of PAR data
# below zero in the minute and daily summary datasets

# please be respectful. If you are going to share outside our group or use this code in any way, let me know who
# with and what they are doing with it :)

source("R/IMO_Light_model.R")
source("R/Process_light_MS8.R")
source("R/Parse_light_MS8.R")
source("R/Calculate_light_MS8.R")
### Light experiment cheat sheert

calc.PAR()

No_Tur=list(
TSS=0.5,
ND.adj=c(1,0.78926,0.078626))

Tur=list(
TSS=10,
ND.adj=c(1,29.35557909,2.935557909))

exp1.dat=rbind(
   data.frame(No_Tur),
   data.frame(Tur))

apply(exp1.dat,MARGIN=1,FUN=function(x){
        calc.PAR(TSS=x["TSS"],ND.adj=x["ND.adj"])})

### Process an MS8 sensor data file
log.files="Data/Florence Bay 29 May to 27 June 2017 deployment2.txt"
#log.dat=list()
#for(f in 1:length(log.files)){
f=1
file.f=log.files[f]

#### this is all you need to run the function
out.dat=proccess.MS8.logger.dat(file.f)
   # note, don't worry too much about these warning messages. This example file has some issues
   # that I don't think are typical. PLEASE NEVER OPEN RAW FILES IN EXCEL, EVER.
   # make a copy first.
names(out.dat)


### Example of new parsing code that allows a data cleaning step ---------------
log.files="Data/Florence Bay 29 May to 27 June 2017 deployment2.txt"
f=1
file.f=log.files[f]
parsed.dat=parse.MS8.logger.dat(file.f)
 # clean the output
 
# now calculate light values from the cleaned data file.
clean.dat=calculate.MS8.logger.dat(tt)


### example light model using a different absorption profile.
#The default absorption data is:
absorption.dat=read.csv("Data/absorption_dat.csv")

# go get the spectral dataset.
AS.dat=read.csv("Data/ActionSpetrum_dat.csv")
ES.dat=read.csv("Data/Experimental_spectrum.csv")

head(AS.dat)
dim(AS.dat)

# currently available options
spectra=unique(AS.dat$study.spectra)
par(mfrow=c(1,1))
plot(NA,xlim=c(400,700),ylim=c(0,1))
for(x in 1:length(spectra)){
 pd=AS.dat[which(AS.dat$study.spectra==spectra[x]),]
 lines(pd$nm,pd$Absorption,col=x)
 }
legend("bottomleft",legend=spectra,lty=1,col=1:length(spectra))
# curernt default = Chlorophyl  a
lines(absorption.dat$nm,absorption.dat$Absorption,col="grey",lty=3,lwd=2)

# use Wangpraseurt and compare to Chlorophyl a
wang.abs.dat=AS.dat[which(AS.dat$study.spectra=="Wangpraseurt"),]
out.dat.wang=proccess.MS8.logger.dat(file.f,absorption.dat=wang.abs.dat)
# plot DLI values
Wang.dat.daily=out.dat.wang$all.dat.daily
Wang.dat.daily=Wang.dat.daily[which(Wang.dat.daily$j.date>0),]
Wang.dat.daily=Wang.dat.daily[which(Wang.dat.daily$PAR.I.max<2500),]  # clearly not possible.

head(Wang.dat.daily)


# use a spectra averaged between Kuhl and Wangpraseurt
Kuhl.abs.dat=AS.dat[which(AS.dat$study.spectra=="Kuhl"),]

new.lambda.vec=400:700
Wang.vals=approx(wang.abs.dat$nm,
                 wang.abs.dat$Absorption,
                 xout=new.lambda.vec,rule=2)$y
Kuhl.vals=approx(Kuhl.abs.dat$nm,
                 Kuhl.abs.dat$Absorption,
                 xout=new.lambda.vec,rule=2)$y
ave.abs.vals=(Wang.vals+Kuhl.vals)/2
ave.abs.dat=data.frame(nm=new.lambda.vec,Absorption=ave.abs.vals)

out.dat.aveSpec=proccess.MS8.logger.dat(file.f,absorption.dat=ave.abs.dat)
ave.dat.daily=out.dat.aveSpec$all.dat.daily
ave.dat.daily=ave.dat.daily[which(ave.dat.daily$j.date>0),]
ave.dat.daily=ave.dat.daily[which(ave.dat.daily$PAR.I.max<2500),]  # clearly not possible.

### weighted average of all spectra
spectra.list=list()
spectra.weights=rep(1/length(spectra),length(spectra))  # equally weighted spectra
for(x in 1:length(spectra)){
 pd=AS.dat[which(AS.dat$study.spectra==spectra[x]),]
 spectra.list=c(spectra.list,
                list(
                approx(pd$nm,pd$Absorption,xout=new.lambda.vec,rule=2)$y*spectra.weights[x]))
 }

weighted.ave.abs.spectra=rowSums(data.frame(spectra.list))

par(mfrow=c(1,1))
plot(NA,xlim=c(400,700),ylim=c(0,1))
for(x in 1:length(spectra)){
 pd=AS.dat[which(AS.dat$study.spectra==spectra[x]),]
 lines(pd$nm,pd$Absorption,col=x)
 }
legend("bottomleft",legend=spectra,lty=1,col=1:length(spectra))
# curernt default = Chlorophyl  a
lines(absorption.dat$nm,absorption.dat$Absorption,col="grey",lty=3,lwd=2)
lines(new.lambda.vec,weighted.ave.abs.spectra,col="orange",lwd=3)


























