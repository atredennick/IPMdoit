# Import and format growth parameters
# then define growth function

# growth parameters
Gpars=list(intcpt=rep(NA,n_spp),intcpt.yr=matrix(0,Nyrs,n_spp), intcpt.gr=matrix(0,6,n_spp),
  slope=rep(NA,n_spp),slope.yr=matrix(0,Nyrs,n_spp),
  nb=matrix(0,n_spp,n_spp),alpha=matrix(NA,n_spp,n_spp),
  sigma2.a=rep(NA,n_spp),sigma2.b=rep(NA,n_spp))
for(i in 1:n_spp){
  infile <- paste("../extdata/Growth_params_",spp_list[i],".csv",sep="")
  Gdata=read.csv(infile)
  Gpars$intcpt[i]=Gdata$Intercept[1]
  tmp=which(names(Gdata)=="Group")
  if(length(tmp)>0) Gpars$intcpt.gr[,i]=Gdata$Group[!is.na(Gdata$Group)] 
  Gpars$intcpt.yr[,i]=Gdata$Intercept.yr
  Gpars$slope[i]=Gdata$logarea.t0[1]
  # random effects on slope
  tmp=which(names(Gdata)=="logarea.t0.yr")
  if(length(tmp)>0) Gpars$slope.yr[,i]=Gdata[,tmp]
  # get competition coefficients
  tmp=paste("crowd",1:length(spp_list),sep="")
  tmp=which(is.element(names(Gdata),tmp))
  if(length(tmp)>0) Gpars$nb[i,]=as.numeric(Gdata[1,tmp])
  Gpars$alpha[i,]=Gdata$alpha[1:length(spp_list)]
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth function
G=function(v,u,W,Gpars,doYear,doSpp){
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+(Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    W%*%(Gpars$nb[doSpp,])
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=dnorm(v,mu,sqrt(sigma2))
  out
}

