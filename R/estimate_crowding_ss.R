###########################################################
#### Functions for calculating mean field genet crowding;
#### Can be with or without genet overlap.
####
#### Andrew Tredennick
#### 1-30-2015

####
#### SINGLE SPECIES VERSIONS
####

#' Calculate crowding with genet overlap allowed
#' @author Andrew Tredennick
#' @aliases crowd_overlap_ss
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
#' @param v
crowd_overlap_ss <- function(A, N, vt, h, alphaG, alphaS, WmatG, WmatS,
                             n_spp, Ctot, Cr, b.r, expv, r.U, v.r, v){
  for(ii in 1:n_spp){ 
    # first do all overlap W's
    Xbar=cover*A/N       # multiply by A to get cover back in cm^2
    varX=varN(v,nt,h,Xbar,N) 
    
    muWG = pi*Xbar*N/(A*alphaG[ii])
    muWS = pi*Xbar*N/(A*alphaS[ii])
    
    muWG[is.na(muWG)]=0
    muWS[is.na(muWS)]=0
    
    WmatG[[ii]]=matrix(muWG,nrow=length(v[[ii]]),ncol=n_spp,byrow=T)
    WmatS[[ii]]=matrix(muWS,nrow=length(v[[ii]]),ncol=n_spp,byrow=T)
    
    # now do conspecific no overlap W
    Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
    Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural")
    
    WmatG[[ii]]=WriG_ss(v.r[[ii]],ii,r.U,Cr,Ctot)/A
    WmatS[[ii]]=WriS_ss(v.r[[ii]],ii,r.U,Cr,Ctot)/A
  }
  return(list(WmatG=WmatG, WmatS=WmatS))
}



#' Calculate crowding with genet overlap not allowed
#' @author Andrew Tredennick
#' @aliases crowd_no_overlap_ss
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
#' @param size.range
crowd_no_overlap_ss <- function(A, vt, h, alphaG, alphaS, WmatG, WmatS,
                                n_spp, Ctot, Cr, b.r, expv, r.U, v.r, size.range){
  for(ii in 1:n_spp){
    Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
    Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural") 
  }
  for(jj in 1:n_spp){ 
    
    WfunG=splinefun(size.range,WriG_ss(size.range,jj,r.U,Cr,Ctot))
    WfunS=splinefun(size.range,WriS_ss(size.range,jj,r.U,Cr,Ctot))
    
    for(ii in 1:n_spp) { 
      WmatG[[ii]]=WfunG(v.r[[ii]])/A 
      WmatS[[ii]]=WfunS(v.r[[ii]])/A 
    }
  }
  return(list(WmatG=WmatG, WmatS=WmatS))
}


