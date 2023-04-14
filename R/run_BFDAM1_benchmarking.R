
# x = A list which has 3 data frames (DF). 
# The 1st DF (betavals) should contain the methylation data where rows are CpG sites and the columns are the methylation values from Group1 and Group2 samples. Note: The 1st column of this DF should be a column called 'windows' that should contain the information about each CpG site and to which window number it was assigned.
# The 2nd DF (cpglocationdf) should contain the annotations of the CpGs obtained using the Illumina manifest file. 
# The 3rd DF (windowinfo) should contain the information about the span (region start - end) of each genomic window created by the user.
# iter = Number of MCMC iterations to be run
# burn = Number of burn-ins
# seed = A random number for reproducibilty of results
# grplab1 = Label for group 1 in methylation data frame/matrix
# grplab2 = Label for group 2 in methylation data frame/matrix

fit.BFDAM1 = function(x, iter = 20000,
                      burn = 1000, seed = 1234,
                      grplab1 = 'Normal', grplab2 = 'Tumor'){
  set.seed(seed)
  control.lab = grplab1
  treat.lab = grplab2
  ptm.start = proc.time()
  bvals = x$betavals
  location = x$cpglocationdf
  clust.locations = x$windowinfo
  nwin<- length(unique(bvals$windows))
  results<-NULL
  for(w in 1:length(unique(bvals$windows))){
    print(w)
    data<- bvals
    datsub=data[data$windows==w,]
    cols=ncol(datsub)
    n=nrow(datsub)
    x1=c(1:n)
    yall=datsub[,2:cols]
    
    #### Log-Transforming Methylation data ####
    
    yall<- as.matrix(logit2(yall))
    
    #### Separating Normal and Cancer groups data ####
    
    ynorm=yall[,grep(control.lab,names(yall))]
    ycanc=yall[,grep(treat.lab,names(yall))]
    
    #### Combined group mean fitting ####
    
    gm_all<-rep(0,n)
    gm_all<- apply(yall,1,mean)
    fitg <- ncs(x1,gm_all,burn,iter)
    overall_mean_fit<-fitg$fhat
    overall_sigma_fit<-fitg$sig2
    overall_tau_fit<-fitg$tau
    lambda_overall<-cbind(overall_tau_fit, overall_sigma_fit)
    
    #### Normal group mean fitting ####
    
    gm_norm<-rep(0,n)
    gm_norm<- apply(ynorm,1,mean)
    fitgnorm <- ncs(x1,gm_norm,burn,iter)
    norm_mean_fit<-fitgnorm$fhat
    norm_sigma_fit<-fitgnorm$sig2
    norm_tau_fit<-fitg$tau
    lambda_norm<-cbind(norm_tau_fit, norm_sigma_fit)
    
    #### Cancer group mean fitting ####
    
    gm_canc<-rep(0,n)
    gm_canc<- apply(ycanc,1,mean)
    fitgcanc <- ncs(x1,gm_canc,burn,iter)
    canc_mean_fit<-fitgcanc$fhat
    canc_sigma_fit<-fitgcanc$sig2
    canc_tau_fit<-fitg$tau
    lambda_canc<-cbind(canc_tau_fit, canc_sigma_fit)
    
    #### Creating Sigma Matrix Combined group ####
    
    sig2.mat_mod<- matrix(0,iter,dim(yall)[2])
    sig_over_mat<- matrix(0,iter,dim(yall)[2]*n)
    sdj<- apply(yall, 2, var)
    for(i in 1:iter)
    {
      for(j in 1:dim(yall)[2])
      {
        sig2.mat_mod[i,j]<- sdj[j] + overall_sigma_fit[i]
      }
    }
    for(i in 1:iter)
    {
      sig_over_mat[i,]=rep(sig2.mat_mod[i,],each=n)
    }
    
    #### Creating Sigma Matrix Normal group ####
    
    sig2.matnorm_mod<- matrix(0,iter,dim(ynorm)[2])
    sig_norm_mat<- matrix(0,iter,dim(ynorm)[2]*n)
    sdj.norm<- apply(ynorm,2,var)
    for(i in 1:iter)
    {
      for(j in 1:dim(ynorm)[2])
      {
        sig2.matnorm_mod[i,j]<- sdj.norm[j] + norm_sigma_fit[i]
      }
    }
    for(i in 1:iter)
    {
      sig_norm_mat[i,]=rep(sig2.matnorm_mod[i,],each=n)
    }
    
    #### Creating Sigma Matrix Cancer group ####
    
    sig2.matcanc_mod<- matrix(0,iter,dim(ycanc)[2])
    sig_canc_mat<- matrix(0,iter,dim(ycanc)[2]*n)
    sdj.canc<- apply(ycanc,2,var)
    for(i in 1:iter)
    {
      for(j in 1:dim(ycanc)[2])
      {
        sig2.matcanc_mod[i,j]<- sdj.canc[j] + canc_sigma_fit[i]
      }
    }
    for(i in 1:iter)
    {
      sig_canc_mat[i,]=rep(sig2.matcanc_mod[i,],each=n)
    }
    
    #### Calculating Marginal Likelihood Overall ####
    
    mean.matall<- t(apply(overall_mean_fit,1,function(x) rep(x,dim(yall)[2])))
    lik_overall<-matrix(0,iter,dim(yall)[2]*n)
    sum_lik_overall<-rep(0,iter)
    for(l in 1:iter){
      for(k in 1:n){
        m<- mean.matall[l,k]
        s<- sqrt(sig_over_mat[l,k])
        lik_overall[l,k]<- log(dnorm(gm_all[k],m,s))
      }
      sum_lik_overall[l]<-sum(lik_overall[l,])
    }
    marginal_lik_overall<-mean(sum_lik_overall)
    
    #### Calculating Marginal Likelihood Normal group ####
    
    mean.matnorm<- t(apply(norm_mean_fit,1,function(x) rep(x,dim(ynorm)[2])))
    lik_norm<-matrix(0,iter,dim(ynorm)[2]*n)
    sum_lik_norm<-rep(0,iter)
    for(l in 1:iter){
      for(k in 1:n){
        m<- mean.matnorm[l,k]
        s<- sqrt(sig_norm_mat[l,k])
        lik_norm[l,k]<- log(dnorm(gm_norm[k],m,s))
      }
      sum_lik_norm[l]<-sum(lik_norm[l,])
    }
    marginal_lik_norm<- mean(sum_lik_norm)
    
    #### Calculating Marginal Likelihood Cancer group ####
    
    mean.matcanc<- t(apply(canc_mean_fit,1,function(x) rep(x,dim(ycanc)[2])))
    lik_canc<-matrix(0,iter,dim(ycanc)[2]*n)
    sum_lik_canc<-rep(0,iter)
    for(l in 1:iter){
      for(k in 1:n){
        m<-mean.matcanc[l,k]
        s<-sqrt(sig_canc_mat[l,k])
        lik_canc[l,k]<-log(dnorm(gm_canc[k],m,s))
      }
      sum_lik_canc[l]<-sum(lik_canc[l,])
    }
    marginal_lik_canc<- mean(sum_lik_canc)
    
    #### Calculating Log Bayes Factors ####
    
    M1_bayes<-marginal_lik_overall
    M2_bayes<-marginal_lik_canc+marginal_lik_norm
    logbayesfactor<-M1_bayes-M2_bayes
    ptm.end<- proc.time()
    time.taken<- (ptm.end - ptm.start)[3]/60
    
    #### Summarizing output ####
    
    tmp<-data.frame(Window=w,
                    M1_bayes=M1_bayes,
                    M2_bayes=M2_bayes,
                    logbayesfactor=logbayesfactor,
                    time.min= time.taken)
    results<-data.frame(rbind(results,tmp))
  }
  
  #### Calling DMRs and Non-DMRs ####
  
  if(is.null(clust.locations$Truth)){
    bayes.cutoff<- -5
    Pred=ifelse(results$logbayesfactor<bayes.cutoff,'Dmr','Non-Dmr')
    len.agg=data.frame(ID=paste('win',1:nwin,sep=''),
                       Length=dim(datsub)[1])
    output=data.frame(Method="BFDAM1",
                      len.agg,
                      M1_Bayes = M1_bayes,
                      M2_Bayes = M2_bayes,
                      logBayes = results$logbayesfactor,
                      predicted=Pred,
                      time.min=results$time.min)
  }else{
    truth<- clust.locations$Truth
    bayes.cutoff<- -5
    Pred=ifelse(results$logbayesfactor<bayes.cutoff,'Dmr','Non-Dmr')
    len.agg=data.frame(ID=paste('win',1:nwin,sep=''),
                       Length=dim(datsub)[1])
    output=data.frame(Method="BFDAM1",
                      len.agg,
                      Truth=truth,
                      M1_Bayes = M1_bayes,
                      M2_Bayes = M2_bayes,
                      logBayes = results$logbayesfactor,
                      predicted=Pred,
                      time.min=results$time.min)
  }
  return(output)
}