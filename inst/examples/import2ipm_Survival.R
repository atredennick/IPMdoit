# import & format survival parameters
# then define survival function

# survival parameters
# survival parameters
Spars=list(intcpt=rep(NA,n_spp),intcpt.yr=matrix(0,Nyrs,n_spp),
           slope=rep(NA,n_spp),slope.yr=matrix(0,Nyrs,n_spp),
           nb=matrix(0,n_spp,n_spp),
           alpha=matrix(NA,n_spp,n_spp))

# nb.yr=array(0,dim=c(n_spp,Nyrs,n_spp)),

for(i in 1:n_spp){
  infile=infile <- paste("../extdata/Surv_params_",spp_list[i],".csv",sep="")
  Sdata=read.csv(infile)
  
  Spars$intcpt[i]=Sdata$Intercept[1]
  
  tmp=which(names(Sdata)=="Group")
  if(length(tmp)>0) Spars$intcpt.gr[,i]=Sdata$Group[!is.na(Sdata$Group)] # get spatial average
  
  tmp=which(names(Sdata)=="Intercept.yr")
  if(length(tmp)>0) Spars$intcpt.yr[,i]=Sdata$Intercept.yr
  
  Spars$slope[i]=Sdata$logarea[1]
  
  # random effects on slope
  tmp=which(names(Sdata)=="logarea.yr")
  if(length(tmp)>0) Spars$slope.yr[,i]=Sdata[,tmp]
  
  # get competition coefficients
#   tmp=paste("crowd",1:length(spp_list),sep="")
#   tmp=which(is.element(names(Sdata),tmp))
#   if(length(tmp)>0) Spars$nb[i,]=as.numeric(Sdata[1,tmp])
  # get crowding coefficients
#     Spars$nb[i]=as.numeric(Sdata$crowd)[1]
    tmp=paste("crowd",1:length(spp_list),sep="")
    tmp=which(is.element(names(Sdata),tmp))
    if(length(tmp)>0) Spars$nb[i,]=as.numeric(Sdata[1,tmp])
  # get crowd X size interaction
#     Spars$slopeXnb[i]=as.numeric(Sdata$logarea.crowd)[1]
  
#   # get competition X size interactions coefficients
#   tmp=paste("logarea.W",1:length(spp_list),sep="")
#   tmp=which(is.element(names(Sdata),tmp))
#   if(length(tmp)>0) Spars$slopeXnb[i,]=as.numeric(Sdata[1,tmp])  
  
#   # get yr random effects on competition
#   tmp=paste("W",1:length(spp_list),".yr",sep="")
#   tmp=which(is.element(names(Sdata),tmp))
#   if(length(tmp)>0) Spars$nb.yr[i,,]=as.matrix(Sdata[,tmp])
  
  Spars$alpha[i,]=Sdata$alpha[1:length(spp_list)]
} # next i
yrList=Sdata$year
rm(Sdata)

## survival function: probability an individual of size u survives  (u is on log scale)

##crowding (w) is based on the discretized size points...(mid-points)
S=function(u,W,Spars,doYear,doSpp){
  mu=Spars$intcpt[doSpp]+Spars$intcpt.yr[doYear,doSpp]+
     (Spars$slope[doSpp]+Spars$slope.yr[doYear,doSpp])*u+
     W%*%(Spars$nb[doSpp,])
#   mu=Spars$intcpt[doSpp]+Spars$intcpt.yr[doYear,doSpp]+
#     (Spars$slope[doSpp]+Spars$slope.yr[doYear,doSpp])*u+
#     W%*%(Spars$nb[doSpp,])+
#     (W*u)%*%Spars$slopeXnb[doSpp]
  
  return(inv.logit(mu))
}

