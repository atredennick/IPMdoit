###########################################################################
#### Competition (crowding) functions
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-29-2015

#' Function to calculate size-dependent crowding for growth, assuming no overlap
#' @aliases get_crowd_growth
#' @param r
#' @param i
#' @param j
#' @param r.U
#' @param Cr
#' @param Ctot
wrijG=function(r,i,j,r.U,Cr,Ctot){
  return(2*pi*integrate(function(z) z*exp(-alphaG[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
           pi*Ctot[j]*exp(-alphaG[i,j]*((r+r.U[j])^2))/alphaG[i,j]);   
}
WrijG=Vectorize(wrijG,vectorize.args="r")

#' Function to calculate size-dependent crowding for survival, assuming no overlap
#' @aliases get_crowd_survival
#' @param r
#' @param i
#' @param j
#' @param r.U
#' @param Cr
#' @param Ctot
wrijS=function(r,i,j,r.U,Cr,Ctot){
  return(2*pi*integrate(function(z) z*exp(-alphaS[i,j]*(z^2))*Cr[[j]](z-r),r,r+r.U[j])$value+
           pi*Ctot[j]*exp(-alphaS[i,j]*((r+r.U[j])^2))/alphaS[i,j]);   
}
WrijS=Vectorize(wrijS,vectorize.args="r")