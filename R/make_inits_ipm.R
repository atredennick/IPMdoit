###########################################################################
#### Function that, when called, creates key variables and storage matrices
#### and vectors for simulating an IPM
####
#### Andrew Tredennick: atredenn@gmail.com
#### Date: 1-29-2015

###########################################################################
#### There are several necessary global variables that must be set
#### and used as inputs for this function.

#' Function to create iteration matrix size and initial vectors for IPM
#' 
#' @author Andrew Tredennick
#' @name make_inits_ipm
#' @aliases make_inits_ipm
#' @param n_spp A scalar for the number of species to simulate.
#' @param iter_matrix_dims A vector whose length is equal to 'n_spp' and represents
#'                        the dimensions of iteration matrix. 
#' @param max_size A vector whose length is equal to 'n_spp' and represents the
#'                 maximum size allowable for a genet from each species.
#' @return Returns a list of matrices and vectors for IPM simulating.

make_inits_ipm <- function(n_spp, iter_matrix_dims, max_size){
  v=v.r=b.r=expv=Cr=WmatG=WmatS=list(n_spp)
  h=r.L=r.U=Ctot=numeric(n_spp)
  for(i in 1:n_spp){
    # minimum (0.9*minimum size from data) and maximum sizes (1.1*maximum size from data)
    L=log(0.2)
    U=log(max_size[i])*1.1     
    
    # boundary points b and mesh points y. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
    b = L+c(0:iter_matrix_dims[i])*(U-L)/iter_matrix_dims[i] 
    
    # v calculates the middle of each n-equal-sized portion.
    v[[i]] = 0.5*(b[1:iter_matrix_dims[i]]+b[2:(iter_matrix_dims[i]+1)])
    
    # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
    h[i] = v[[i]][2]-v[[i]][1]  
    
    # variables for Wr approximation
    b.r[[i]]=sqrt(exp(b)/pi)
    v.r[[i]]=sqrt(exp(v[[i]])/pi)
    expv[[i]]=exp(v[[i]])
    r.L[i] = sqrt(exp(L)/pi)
    r.U[i] = sqrt(exp(U)/pi)
    WmatG[[i]]=matrix(NA,length(v.r[[i]]),n_spp)  # storage of size-specific W values for each focal species
    WmatS[[i]]=matrix(NA,length(v.r[[i]]),n_spp)
    
  } # next species
  rm(i)
  tmp=range(v.r)
  size_range=seq(tmp[1],tmp[2],length=50) # range across all possible sizes
  return(list(v=v, L=L, U=U, h=h, WmatG=WmatG, WmatS=WmatS, expv=expv, 
              b=b, size_range=size_range, v.r=v.r, b.r=b.r, Ctot=Ctot, Cr=Cr,
              r.U=r.U))
}


