###########################################################
#### Functions for calculating mean field genet crowding;
#### Can be with or without genet overlap.
####
#### Andrew Tredennick
#### 1-30-2015

#' Calculate crowding with genet overlap allowed
#' @author Andrew Tredennick
#' @aliases crowd_overlap
#' @param A Area of the quadrat.
#' @param N Current population density vector.
#' @param v
#' @param nt Current population vector.
#' @param h
#' @param alphaG Alpha parameters (vector with n=n_spp) for growth.
#' @param alphaS Alpha parameters (vector with n=n_spp) for survival.
#' @param WmatG
#' @param WmatS
#' @param n_spp Number of species.
#' @param Ctot
#' @param Cr
#' @param b.r
#' @param expv
#' @param r.U
#' @param v.r
crowd_overlap <- function(A, N, vt, h, alphaG, alphaS, WmatG, WmatS,
                          n_spp, Ctot, Cr, b.r, expv, r.U, v.r){
  for(ii in 1:n_spp){ 
    # first do all overlap W's
    Xbar=cover*A/N       # multiply by A to get cover back in cm^2
    varX=varN(inits$v,nt,inits$h,Xbar,N) 
    
    muWG = pi*Xbar*N/(A*alphaG[ii,])
    muWS = pi*Xbar*N/(A*alphaS[ii,])
    
    muWG[is.na(muWG)]=0
    muWS[is.na(muWS)]=0
    
    inits$WmatG[[ii]]=matrix(muWG,nrow=length(inits$v[[ii]]),ncol=n_spp,byrow=T)
    inits$WmatS[[ii]]=matrix(muWS,nrow=length(inits$v[[ii]]),ncol=n_spp,byrow=T)
    
    # now do conspecific no overlap W
    inits$Ctot[ii]=inits$h[ii]*sum(inits$expv[[ii]]*nt[[ii]]) 
    inits$Cr[[ii]]=splinefun(inits$b.r[[ii]],inits$h[ii]*c(0,cumsum(inits$expv[[ii]]*nt[[ii]])),method="natural")
    
    inits$WmatG[[ii]][,ii]=WrijG(inits$v.r[[ii]],ii,ii,inits$r.U,inits$Cr,inits$Ctot)/A
    inits$WmatS[[ii]][,ii]=WrijS(inits$v.r[[ii]],ii,ii,inits$r.U,inits$Cr,inits$Ctot)/A
  }
  return(list(WmatG=WmatG, WmatS=WmatS))
}



#' Calculate crowding with genet overlap not allowed
#' @author Andrew Tredennick
#' @aliases crowd_no_overlap
#' @param A Area of the quadrat.
#' @param nt Current population vector.
#' @param h
#' @param alphaG Alpha parameters (vector with n=n_spp) for growth.
#' @param alphaS Alpha parameters (vector with n=n_spp) for survival.
#' @param WmatG
#' @param WmatS
#' @param n_spp Number of species.
#' @param Ctot
#' @param Cr
#' @param b.r
#' @param expv
#' @param r.U
#' @param v.r
crowd_no_overlap <- function(A, vt, h, alphaG, alphaS, WmatG, WmatS,
                             n_spp, Ctot, Cr, b.r, expv, r.U, v.r){
  for(ii in 1:n_spp){
    Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
    Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural") 
  }
  for(jj in 1:n_spp){ 
    
    WfunG=splinefun(size.range,WrijG(size.range,jj,jj,r.U,Cr,Ctot))
    WfunS=splinefun(size.range,WrijS(size.range,jj,jj,r.U,Cr,Ctot))
    
    for(ii in 1:n_spp) { 
      WmatG[[ii]][,jj]=WfunG(v.r[[ii]])/A 
      WmatS[[ii]][,jj]=WfunS(v.r[[ii]])/A 
    }
  }
  return(list(WmatG=WmatG, WmatS=WmatS))
} # end NoOverlap if