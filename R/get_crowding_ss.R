###########################################################################
#### Competition (crowding) functions
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-29-2015

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
wriG_ss=function(r,i,r.U,Cr,Ctot){
  return(2*pi*integrate(function(z) z*exp(-alphaG[i]*(z^2))*Cr[[i]](z-r),r,r+r.U[i])$value+
           pi*Ctot[i]*exp(-alphaG[i]*((r+r.U[i])^2))/alphaG[i]);   
}
WriG_ss=Vectorize(wriG_ss,vectorize.args="r")

#' Function to calculate size-dependent crowding for survival, assuming no overlap
#' @aliases get_crowd_survival_ss
#' @param r
#' @param i
#' @param j
#' @param r.U
#' @param Cr
#' @param Ctot
wriS_ss=function(r,i,r.U,Cr,Ctot){
  return(2*pi*integrate(function(z) z*exp(-alphaS[i]*(z^2))*Cr[[i]](z-r),r,r+r.U[i])$value+
           pi*Ctot[i]*exp(-alphaS[i]*((r+r.U[i])^2))/alphaS[i]);   
}
WriS_ss=Vectorize(wriS_ss,vectorize.args="r")