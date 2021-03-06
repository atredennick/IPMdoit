###############################################################
#### Script to run a multi-species IPM using Idaho data and
#### the 'IPMdoit' package.
####
#### Andrew Tredennick
#### 1-29-2015

# Clear the workspace
rm(list=ls(all=TRUE))

####
#### Load required packages -----------------------
####
library(IPMdoit); library(ggplot2)
library(boot); library(mvtnorm)
library(msm); library(statmod)
library(plyr); library(reshape2)
library(gridExtra)


####
#### Set up global parameters ---------------------
####
A <- 10000 #Area of 100cm x 100cm quadrat
tlimit <- 100  ## number of years to simulate
burn_in <- 50    # years to cut before calculations
spp_list <- c("ARTR","HECO","POSE","PSSP")
iter_matrix_dims <- c(75,75,50,50)     #Set matrix dimension for each species
max_size <- c(3000,202,260,225)    # in cm^2: PSSP=225 HECO=202  POSE=260  ARTR=3000  # minSize=0.2  cm^2
Nyrs <- 22
doGroup <- NA  # NA for spatial avg., values 1-6 for a specific group
n_spp <- length(spp_list)
NoOverlap_Inter=F

####
#### Source import scripts to bring in vital rate functions and parameters ------------- 
####
source("import2ipm_Growth.R")
source("import2ipm_Survival.R")
source("import2ipm_Recruitment.R")

####
#### Get initial vectors and matrices built --------------
####
inits <- make_inits_ipm(n_spp = n_spp, 
                        iter_matrix_dims = iter_matrix_dims, 
                        max_size = max_size)


####
#### Run simulation -----------------
####
# initial population density vector
nt <- inits$v

# loop through species really quick to set low initial densities
for(i in 1:n_spp) nt[[i]][]=0.001
new.nt <- nt #set initial density vector to be fed into IPM

# set up matrix to record cover
covSave <- matrix(NA,tlimit,n_spp)
covSave[1,] <- sum_cover(inits$v,nt,inits$h,A)

# set up list to store size distributions
sizeSave <- list(NULL)
for(i in 1:n_spp){
  sizeSave[[i]] <- matrix(NA,length(inits$v[[i]]),(tlimit))
  sizeSave[[i]][,1] <- nt[[i]]/sum(nt[[i]])
}

# initial densities 
Nsave <- matrix(NA,tlimit,n_spp)
Nsave[1,] <- sum_N(nt,inits$h)

yrSave <- rep(NA,tlimit)

# Loop through simulation times and iterate population
pb <- txtProgressBar(min=2, max=tlimit, char="+", style=3, width=65)
for (t in 2:(tlimit)){
  #draw from observed year effects
  allYrs <- c(1:Nyrs)
  doYear <- sample(allYrs,1)
  yrSave[t] <- doYear
  
  #get recruits per area
  cover <- covSave[t-1,]
  N <- Nsave[t-1,]
  recs_per_area <- get_rpa(Rpars,cover,doYear)
  
  #calculate size-specific crowding
  alphaG <- Gpars$alpha 
  alphaS <- Spars$alpha 
  if(NoOverlap_Inter==F){#T: heterospecific genets cannot overlap; F: overlap allowed
    crowd_list <- crowd_overlap(A, N, inits$vt, inits$h, alphaG, alphaS, inits$WmatG, inits$WmatS,
                                n_spp, inits$Ctot, inits$Cr, inits$b.r, inits$expv, inits$r.U, inits$v.r, inits$v)
  }else{
    crowd_list <- crowd_no_overlap(A, inits$vt, inits$h, alphaG, alphaS, inits$WmatG, inits$WmatS,
                                   n_spp, inits$Ctot, inits$Cr, inits$b.r, inits$expv, inits$r.U, inits$v.r,
                                   inits$size_range)
  } # end NoOverlap if
  
  for(doSpp in 1:n_spp){  
    if(cover[doSpp]>0){    
      # make kernels and project
      K_matrix=make_K_matrix(inits$v[[doSpp]],crowd_list$WmatG[[doSpp]],crowd_list$WmatS[[doSpp]],
                             Rpars,recs_per_area,Gpars,Spars,doYear,doSpp,inits$h,demo_stoch=FALSE)	
      new.nt[[doSpp]]=K_matrix%*%nt[[doSpp]] 
      sizeSave[[doSpp]][,t]=new.nt[[doSpp]]/sum(new.nt[[doSpp]])  
    }    
  } # next species
  
  nt=new.nt 
  covSave[t,]=sum_cover(inits$v,nt,inits$h,A)  # store the cover as cm^2/cm^2
  Nsave[t,]=sum_N(nt,inits$h)
  
  setTxtProgressBar(pb, t)
  flush.console()
  if(sum(is.na(nt))>0) browser()  
} # next time step


####
#### Some example plots to visually check results ----------------------------
####
cover <- as.data.frame(100*covSave[(burn_in+1):tlimit,])
colnames(cover) <- spp_list
cover$year <- seq((burn_in+1),tlimit)
cover_df <- melt(cover, id.vars = "year")
colnames(cover_df) <- c("Sim_Year", "Species", "Cover")

g_cover_box <- ggplot(cover_df)+
  geom_boxplot(aes(x=Species, y=Cover, fill=Species))+
  guides(fill=FALSE)
g_cover_ts <- ggplot(cover_df)+
  geom_line(aes(x=Sim_Year, y=Cover, color=Species, group=Species))
grid.arrange(g_cover_box, g_cover_ts, nrow=1, ncol=2)

