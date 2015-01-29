##################################################
#### A bunch of utility functions.
####
#### Andrew Tredennick
#### 1-29-2015


#' Function to sum total cover of each species
#' @author Andrew Tredennick
#' @aliases sum_cover
#' @param v
#' @param nt The population vector.
#' @param h
#' @param A Area of the quadrat (scalar)
sum_cover=function(v,nt,h,A){
  out=lapply(1:n_spp,function(i,v,nt,h,A) h[i]*sum(nt[[i]]*exp(v[[i]]))/A,v=v,nt=nt,h=h,A=A)
  return(unlist(out))
} 


#' Function to sum total density of each species
#' @author Andrew Tredennick
#' @aliases sum_N
#' @param nt The population vector.
#' @param h
sum_N=function(nt,h){
  out=lapply(1:n_spp,function(i,nt,h) h[i]*sum(nt[[i]]),nt=nt,h=h)
  return(unlist(out))
}

#' Function to calculate size variance of each species
#' @author Andrew Tredennick
#' @aliases var_N
#' @param v
#' @param nt The population vector.
#' @param h
#' @param Xbar
#' @param N Density vector for each size class.
varN=function(v,nt,h,Xbar,N){
  out=lapply(1:n_spp,function(i,v,nt,h,Xbar,N) h[i]*sum((exp(v[[i]]-Xbar[i])^2)*nt[[i]])/N[i],v=v,nt=nt,h=h,Xbar=Xbar,N=N)
  return(unlist(out))
}  
