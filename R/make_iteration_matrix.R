###########################################################################
#### Function that, when called, creates key variables and storage matrices
#### and vectors for simulating an IPM
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-29-2015

####
#### Combined kernel ------------------------------------
####
#' Make iteration matrix
#' 
#' @aliases make_K_matrix
#' @author Andrew Tredennick
#' @param u
#' @param v
#' @param muWG
#' @param muWS
#' @param rec_params
#' @param recs_per_area
#' @param growth_params
#' @param surv_params
#' @param do_year
#' @param do_spp
make_K_values=function(v,u,muWG,muWS, #state variables
                       rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp){  #growth arguments
  f(v,u,rec_params,recs_per_area,do_spp)+S(u,muWS,surv_params,do_year,do_spp)*G(v,u,muWG,growth_params,do_year,do_spp) 
}

####
#### Function to format the W matrix for the outer product ---------
####
#' Expand crowding matrix
#' 
#' @aliases expand_W_matrix
#' @author Andrew Tredennick
#' @param u
#' @param v
#' @param W
expand_W_matrix=function(v,u,W){
  if(dim(W)[1]!=length(u)) stop("Check size of W")
  n_spp=dim(W)[2]
  W=as.vector(W)
  W=matrix(W,length(W),ncol=length(v))
  W=as.vector(t(W))
  W=matrix(W,nrow=length(u)*length(v),ncol=n_spp)
  return(W)
}

####
#### Function to make iteration matrix based only on mean crowding ----------
####
#' Make iteration matrix based only on mean crowding
#' 
#' @aliases make_K_matrix
#' @author Andrew Tredennick
#' @param v
#' @param muWG
#' @param muWS
#' @param rec_params
#' @param recs_per_area
#' @param growth_params
#' @param surv_params
#' @param do_year
#' @param do_spp
#' @param h
make_K_matrix=function(v,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp,h) {
  muWG=expand_W_matrix(v,v,muWG)
  muWS=expand_W_matrix(v,v,muWS)
  
  K.matrix=outer(v,v,make_K_values,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp)
  return(h[do_spp]*K.matrix)
}



####
#### Function to make iteration matrix based only on mean crowding ----------
####
#' Make iteration matrix based only on mean crowding: single species
#' 
#' @aliases make_K_matrix_ss
#' @author Andrew Tredennick
#' @param v
#' @param muWG
#' @param muWS
#' @param rec_params
#' @param recs_per_area
#' @param growth_params
#' @param surv_params
#' @param do_year
#' @param do_spp
#' @param h
make_K_matrix_ss=function(v,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp,h) {  
  muWG=rep(muWG,length(v))
  muWS=rep(muWS,length(v))
  K.matrix=outer(v,v,make_K_values,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp)
  return(h[do_spp]*K.matrix)
}
