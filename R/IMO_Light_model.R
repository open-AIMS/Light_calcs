


calc.PAR=function(SZA=0,
                  TSS=0,
                  Depth=0,
                  solar.zenith.dat=NA,
                  model.dat=NA,
                  absorption.dat=NA,
                  ND.adj=1,
                  nm.range=c(400,700)){
  #SZA=0;TSS=0.5;Depth=5;solar.zenith.dat=NA;model.dat=NA;absorption.dat=NA;nm.range=c(400,700)
  #install.packages("NISTunits", dependencies = TRUE)
  require(NISTunits)
  #NISTdegTOradian(180)
  #NISTradianTOdeg(pi)

  # read in supporting data
  if(length(solar.zenith.dat)==1){
   solar.zenith.dat=read.csv("U:/1 TOWNSVILLE PORT - CLEVELAND BAY/1  REBECCA FISHER ANALYSES/RScripts/solar_zenith_dat.csv")
  }
  if(length(model.dat)==1){
   model.dat=read.csv("U:/1 TOWNSVILLE PORT - CLEVELAND BAY/1  REBECCA FISHER ANALYSES/RScripts/model_dat.csv")
  }
  if(length(absorption.dat)==1){
   absorption.dat=read.csv("U:/1 TOWNSVILLE PORT - CLEVELAND BAY/1  REBECCA FISHER ANALYSES/RScripts/absorption_dat.csv")
  }

  # extract variables
  Solar.Zenith.Angle=SZA

  sza=solar.zenith.dat$sza  #P9
  frac=solar.zenith.dat$frac   #Q9

  Es_C0=model.dat$Es_C0
  Es_C1=model.dat$Es_C1
  Es_C2=model.dat$Es_C2
  Es_C3=model.dat$Es_C3

  EDpve=(Es_C0+Es_C1*Solar.Zenith.Angle+Es_C2*Solar.Zenith.Angle^2+
     Es_C3*Solar.Zenith.Angle^3)*cos(NISTdegTOradian(Solar.Zenith.Angle))

  fracval=frac[which(sza==Solar.Zenith.Angle)]#sum(R9:R99)#fracval

  EDnve=EDpve*fracval

  Kd_M=model.dat$Kd_M
  Kd_C=model.dat$Kd_C

  NTU_C=0.89848512
  NTU_M=0.70246627
  NTU=NTU_M*TSS+NTU_C

  Kd=Kd_M*NTU+Kd_C

  microW.cm2.nm=(10^2)*EDnve*exp(-Kd*Depth)*ND.adj

  h=6.63E-34
  c.val=3.00E+08
  avo=6.02E+17
  lambda=model.dat$nm

  absorption.dat=merge(data.frame(nm=lambda),absorption.dat,all=T)
  absorption.vals=approx(absorption.dat$nm,
                                          absorption.dat$Absorption,
                                          xout=lambda,rule=2)$y
  W.m2.mn=microW.cm2.nm/100

  microM.sec=(W.m2.mn/((h*c.val)/(lambda*0.000000001)))/avo

  # sum.index
  sum.index=match(nm.range[1]:nm.range[2],lambda)
  PAR=sum(microM.sec[sum.index])
  PUR=sum((microM.sec*absorption.vals)[sum.index])

  return(data.frame(PAR=PAR,PUR=PUR))
}


calc.PAR()


