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
#' @param demographic_stochasticity


make_K_values=function(v,u,muWG,muWS, #state variables
                       rec_params,recs_per_area,growth_params,surv_params,
                       do_year,do_spp){  #growth arguments
  fecundity <- f(v,u,rec_params,recs_per_area,do_spp)
  survival <- S(u,muWS,surv_params,do_year,do_spp)
  growth <- G(v,u,muWG,growth_params,do_year,do_spp)
  fecundity+survival*growth
}



####
####  FXNS for demographic stochasticity
####
####
####  Pairwise multiplication of the population vector
####
#' Get pairs
#' 
#' @aliases get_pairs
#' @author Andrew Tredennick
#' @param X 
#' @param pop_vector
get_pairs <- function(X, pop_vector){
  pairs <- expand.grid(X, X)
  #   pairs$tag <- pairs[,1] - pairs[,2]
  pairs$multi <- pairs[,1]*pairs[,2]*pop_vector
  return(pairs$multi)
}  

####
####  Calculate population vector covariance matrix
####
#' Get nt covariance
#' 
#' @aliases get_pairs
#' @author Andrew Tredennick
#' @param K Combined survival*growth iteration matrix
get_cov <- function(K){
  test <- apply(K, MARGIN = 2, FUN = "get_pairs", 
                pop_vector=(popv))
  mat_dim <- sqrt(dim(test)[1])
  test <- as.data.frame(test)
  test$tag <- rep(c(1:mat_dim), each=mat_dim)
  cov_str <- matrix(ncol=mat_dim, nrow=mat_dim)
  for(do_i in 1:mat_dim){
    tmp <- subset(test, tag==do_i) #subset out the focal i
    rmtmp <- which(colnames(tmp)=="tag") #get rid of id column
    # Sum over k columns
    cov_str[do_i,] <- (-h^2) * apply(tmp[,-rmtmp], MARGIN = 2, FUN = "sum")
  }
  diag(cov_str) <- 1
  return(cov_str)
}


####
#### Generate Multivariate Poisson vector
####
#' Generate poisson vector
#' 
#' @aliases GenerateMultivariatePoisson
#' @author Andrew Tredennick
#' @param pD Length of the population vector
#' @param samples Number of random samples to generate at each vector element (defaults to 1)
#' @param R Population vector covariance matrix from 'get_cov'
#' @param lambda Surival*growth contribution (iteration matrix * population vector)
GenerateMultivariatePoisson<-function(pD, samples=1, R, lambda){
  normal_mu=rep(0, pD)
  normal = mvrnorm(samples, normal_mu, R)
  pois = normal
  p=pnorm(normal)
  for (s in 1:pD){pois[s]=qpois(p[s], lambda[s])}
  return(pois)
}

make.R.values=function(v,u, #state variables
                       Rpars,rpa,doYear,doSpp){
  f(v,u,Rpars,rpa,doSpp)
}

make.P.values <- function(v,u,muWG,muWS, #state variables
                          Gpars,Spars,doYear,doSpp){  #growth arguments
  S(u,muWS,Spars,doYear,doSpp)*G(v,u,muWG,Gpars,doYear,doSpp) 
}

make.P.matrix <- function(v,muWG,muWS,Gpars,Spars,doYear,doSpp,h) {
  muWG=expand_W_matrix(v,v,muWG)
  muWS=expand_W_matrix(v,v,muWS)
  
  P.matrix=outer(v,v,make.P.values,muWG,muWS,Gpars,Spars,doYear,doSpp)
  return(h*P.matrix)
} 

make.R.matrix=function(v,Rpars,rpa,doYear,doSpp,h) {
  R.matrix=outer(v,v,make.R.values,Rpars,rpa,doYear,doSpp)
  return(h*R.matrix)
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
make_K_matrix=function(v,muWG,muWS,rec_params,recs_per_area,growth_params,
                       surv_params,do_year,do_spp,h,demo_stoch=FALSE,popv) {
  if(demo_stoch==FALSE){
    muWG=expand_W_matrix(v,v,muWG)
    muWS=expand_W_matrix(v,v,muWS)
    K.matrix=outer(v,v,make_K_values,muWG,muWS,rec_params,recs_per_area,
                   growth_params,surv_params,do_year,do_spp)
  }
  
  if(demo_stoch==TRUE){
    P.matrix <- make.P.matrix(v,muWG,muWS,Gpars,Spars,doYear,doSpp,h=h)  
    R.matrix <- make.R.matrix(v,Rpars,recs_per_area,doYear,doSpp,h=h)  
    covmat <- get_cov(K=P.matrix)
    pCont <- GenerateMultivariatePoisson(pD = length(popv),
                                         samples = 1,
                                         R = covmat,
                                         lambda = P.matrix%*%popv)
    rCont <- rpois(length(popv),R.matrix%*%popv)
    K.matrix=list(pCont, rCont)
  }
  
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
