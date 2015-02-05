###########################################################################
#### Competition (crowding) functions
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-29-2015

###FROM CHENGJIN'S CODE
wriG=function(r){
  return(2*pi*integrate(function(z) z*exp(-alphaG*(z^2))*Cr(z-r),r,r+r.U)$value+
           pi*Ctot*exp(-alphaG*((r+r.U)^2))/alphaG)
}
WriG=Vectorize(wriG,vectorize.args="r") #Wri, and wri

wriS=function(r){
  return(2*pi*integrate(function(z) z*exp(-alphaS*(z^2))*Cr(z-r),r,r+r.U)$value+
           pi*Ctot*exp(-alphaS*((r+r.U)^2))/alphaS)
}
WriS=Vectorize(wriS,vectorize.args="r") #Wri, and wri



####
#### SINGLE SPECIES VERSIONS
####

#' Function to calculate size-dependent crowding for growth, assuming no overlap
#' @aliases get_crowd_growth_ss
#' @param r
#' @param i
#' @param j
#' @param r.U
#' @param Cr
#' @param Ctot
wrijG_ss=function(r,i,r.U,Cr,Ctot){
  return(2*pi*integrate(function(z) z*exp(-alphaG[i]*(z^2))*Cr[[i]](z-r),r,r+r.U[i])$value+
           pi*Ctot[i]*exp(-alphaG[i]*((r+r.U[i])^2))/alphaG[i]);   
}
WrijG_ss=Vectorize(wrijG_ss,vectorize.args="r")

#' Function to calculate size-dependent crowding for survival, assuming no overlap
#' @aliases get_crowd_survival_ss
#' @param r
#' @param i
#' @param j
#' @param r.U
#' @param Cr
#' @param Ctot
wrijS_ss=function(r,i,r.U,Cr,Ctot){
  return(2*pi*integrate(function(z) z*exp(-alphaS[i]*(z^2))*Cr[[i]](z-r),r,r+r.U[i])$value+
           pi*Ctot[i]*exp(-alphaS[i]*((r+r.U[i])^2))/alphaS[i]);   
}
WrijS_ss=Vectorize(wrijS,vectorize.args="r")