################################################################
####  Functions for demographic stochasticity approximation ####
####                                                        ####
#### Andrew Tredennick: atredenn@gmail.com                  ####
#### Date: 7-15-2015                                        ####
################################################################


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
                pop_vector=(nt[[doSpp]]))
  mat_dim <- sqrt(dim(test)[1])
  test <- as.data.frame(test)
  test$tag <- rep(c(1:mat_dim), each=mat_dim)
  cov_str <- matrix(ncol=mat_dim, nrow=mat_dim)
  for(do_i in 1:mat_dim){
    tmp <- subset(test, tag==do_i) #subset out the focal i
    rmtmp <- which(colnames(tmp)=="tag") #get rid of id column
    # Sum over k columns
    cov_str[do_i,] <- (-h[doSpp]^2) * apply(tmp[,-rmtmp], MARGIN = 2, FUN = "sum")
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

