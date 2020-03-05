#    Copyright 2020 Australian Institute of Marine Science and In-situ Marine Optics
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

#' calc.PAR
#'
#' Calculates PAR and PUR for a given depth and TSS level, based on solar zenthith angle, observed model data collected in Cleveland Bay and an absportion spectra
#'
#' @param  SZA Solar Zenthis Angle
#' 
#' @param TSS Total suspended solids level in the water column (generally derived as conversion factor from NTU readings)
#' 
#' @param Depth water depth at which the light reading should be estimated
#' 
#' @param solar.zenith.dat Fraction of illumination given a solar zenith angle (SZA)
#' 
#' @param model.dat Model data derived from IMO light measurements
#'
#' @param absorption.dat A data.frame containing an absorption spectrum to use for calculating PUR. Must contain columns fo nm and normalized absoption/action (values between zero and one)
#' 
#' @param ND.adj xxx
#' 
#' @param nm.range Light eavelength range over which PUR and PAR will be calculated
#'
#' @export
#' @return Calculated PUR and PAR light data.


calc.PAR=function(SZA=0,
                  TSS=0,
                  Depth=0,
                  solar.zenith.dat=NA,
                  model.dat=NA,
                  absorption.dat=NA,
                  ND.adj=1,
                  nm.range=c(400,700)){
  #SZA=0;TSS=0.5;Depth=5;solar.zenithDat=NA;modelDat=NA;absorptionDat=NA;nm.range=c(400,700)
  #install.packages("NISTunits", dependencies = TRUE)
  require(NISTunits)
  #NISTdegTOradian(180)
  #NISTradianTOdeg(pi)

  # read in supporting data

  if(length(model.dat)==1){
   model.dat=Light.calcs::modelDat
  }
  if(length(absorption.dat)==1){
   absorption.dat=Light.calcs::absorptionDat
  }
  if(length(solar.zenith.dat)==1){
    solar.zenith.dat=Light.calcs::solarZenithDat
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


