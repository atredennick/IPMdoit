###########################################################################
#### Function that, when called, creates key variables and storage matrices
#### and vectors for simulating an IPM
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-29-2015

####
#### Combined kernel ------------------------------------
####
#' Make iteration values
#' 
#' @aliases make_K_values_clim
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
make_K_values_clim=function(v,u,muWG,muWS, #state variables
                       rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp,climate){  #growth arguments
  f(v,u,rec_params,recs_per_area,do_spp)+S(u,muWS,surv_params,do_year,do_spp,climate)*G(v,u,muWG,growth_params,do_year,do_spp,climate) 
}

####
#### Function to make iteration matrix based only on mean crowding ----------
####
#' Make iteration matrix based on mean crowding and climate
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
make_K_matrix_clim=function(v,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp,h,climate) {
  muWG=expand_W_matrix(v,v,muWG)
  muWS=expand_W_matrix(v,v,muWS)
  
  K.matrix=outer(v,v,make_K_values,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp,climate)
  return(h[do_spp]*K.matrix)
}



####
#### Function to make iteration matrix based only on mean crowding ----------
####
#' Make iteration matrix based on mean crowding and climate: single species
#' 
#' @aliases make_K_matrix_ss_clim
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
make_K_matrix_ss_clim=function(v,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp,h,climate) {  
  muWG=rep(muWG,length(v))
  muWS=rep(muWS,length(v))
  K.matrix=outer(v,v,make_K_values_clim,muWG,muWS,rec_params,recs_per_area,growth_params,surv_params,do_year,do_spp,climate)
  return(h[do_spp]*K.matrix)
}
