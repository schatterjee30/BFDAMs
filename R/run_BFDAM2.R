
# iter = Number of MCMC iterations to be run
# burn = Number of burn-ins
# seed = A random number for reproducibilty of results
# grplab1 = Label for group 1 in methylation data frame(beta.df)
# grplab2 = Label for group 2 in methylation data frame(beta.df)

#### Required Libraries ####

library(readr)

#### The BFDAM2 Function ####

fit.BFDAM2 = function(beta.df, chr.pos, nCpGs = 100, iter = 20000, burn = 1000, 
                      seed = 1234, grplab1 = 'Normal', grplab2 = 'Tumor'){
  
  set.seed(seed)
  ptm.start<-proc.time()
  
  #### Organizing Input Data ####
  
  beta.df<-beta.df[!duplicated(rownames(beta.df)),]
  chr.pos<-chr.pos[!duplicated(rownames(chr.pos)),]
  beta.df<-beta.df[rownames(beta.df) %in% rownames(chr.pos),]
  if(is.character(chr.pos$chr)){
    chr.pos$chr<-readr::parse_number(chr.pos$chr)
  }
  if(is.character(chr.pos$pos)){
    chr.pos$pos<-as.numeric(chr.pos$pos)
  }
  
  #### Assigning Equally Spaced and Sized Windows ####
  
  split.chrpos<-split(chr.pos, chr.pos$chr)
  chrpos.win<-lapply(split.chrpos, function(x){
    nCpGs<-100
    x$windows<-rep(seq_len(nrow(x)%/%nCpGs+1L),each=nCpGs,len=nrow(x))
    return(x)
  })
  
  # data<-beta.df
  # location<-location.df
  # clust.locations<-window.df
  # nwin<-length(unique(bvals$windows))
  # datsub<-data[data$window == 1,]
  # cols<-ncol(datsub)
  # n<-nrow(datsub)
  # x1<-c(1:n)
  # yall<-datsub[, 2:cols]
  
  #### Log-Transforming Methylation data ####
  
  yall = as.matrix(logit2(yall))
  
  #### Separating Normal and Cancer groups data ####
  
  ynorm<-yall[,grep(grplab1,names(yall))]
  ycanc<-yall[,grep(grplab2,names(yall))]
  
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
  
  #### Summarizing output and projection particles ####
  
  results.indep<-data.frame(Window=1,
                            M1_bayes=M1_bayes,
                            M2_bayes=M2_bayes,
                            logbayesfactor=logbayesfactor,
                            time.min = time.taken)
  results.pass<-list(lambda_overall=lambda_overall,lambda_norm=lambda_norm,lambda_canc=lambda_canc)
  
  #### Dependent Regions modeling ####
  
  results.dep<-NULL
  for(w in 2:length(unique(bvals$windows))){
    print(paste0("window ", w, " is running"))
    
    #### Extracting Inputs ####
    
    ptm.dep.start<- proc.time()
    lambda_overall=results.pass$lambda_overall
    lambda_norm=results.pass$lambda_norm
    lambda_canc=results.pass$lambda_canc
    
    #### Organizing input data ####
    
    datsub=data[data$window==w,]
    cols=ncol(datsub)
    nt=nrow(datsub)
    xt=c(1:nt)
    yallt=datsub[,2:cols]
    
    #### Log-Transforming Methylation data ####
    
    yallt<- as.matrix(logit2(yallt))
    gm_allt<-rep(0,nt)
    gm_allt<- apply(yallt,1,mean)
    
    #### Separating subsequent Normal and Cancer groups data ####
    
    ynormt=yallt[,grep(control.lab,names(yallt))]
    ycanct=yallt[,grep(treat.lab,names(yallt))]
    
    #### Gathering Means of subsequent Normal and Cancer groups data ####
    
    gm_normt<-rep(0,nt)
    gm_canct<-rep(0,nt)
    gm_normt<-apply(ynormt,1,mean)
    gm_canct<-apply(ycanct,1,mean)
    
    #### Projection part overall ####
    
    Nt.overall1<-nrow(lambda_overall)
    e1<-rnorm(Nt.overall1,0,0.1)
    e2<-rnorm(Nt.overall1,0,0.1)
    pltau<- -log(lambda_overall[,1])+e1  
    plsig2<-log(lambda_overall[,2])+e2  
    particles_overall1<-cbind(exp(-pltau),exp(plsig2))
    
    #### Projection part Normal ####
    
    Nt.norm1<-nrow(lambda_norm)
    e1<-rnorm(Nt.norm1,0,0.1)
    e2<-rnorm(Nt.norm1,0,0.1)
    pltau<- -log(lambda_norm[,1])+e1  
    plsig2<-log(lambda_norm[,2])+e2  
    particles_norm1<-cbind(exp(-pltau),exp(plsig2))
    
    #### Projection part Cancer ####
    
    Nt.canc1<-nrow(lambda_canc)
    e1<-rnorm(Nt.canc1,0,0.1)
    e2<-rnorm(Nt.canc1,0,0.1)
    pltau<- -log(lambda_canc[,1]) + e1  
    plsig2<-log(lambda_canc[,2]) + e2  
    particles_canc1<-cbind(exp(-pltau),exp(plsig2))
    
    #### Getting weights for particles overall ####
    
    weight_overall<-rep(0,Nt.overall1)
    K_overall<-getK(xt)
    dim(gm_allt)<-c(NROW(gm_allt),NCOL(gm_allt))
    gen.weight.overall<-getweight(particles_overall1,gm_allt,K_overall)
    if (sum(is.na(gen.weight.overall))>0){
      ll<-sum(is.na(gen.weight.overall))
      gen.weight.overall[which(is.na(gen.weight.overall))]<-rep(mean(gen.weight.overall,na.rm = T),ll)
    }
    overall.weight1<-weight_overall+gen.weight.overall
    
    #### Getting weights for particles Normal ####
    
    weight_norm <- rep(0,Nt.norm1)
    K_norm <- getK(xt)
    dim(gm_normt) <- c(NROW(gm_normt),NCOL(gm_normt))
    gen.weight.norm <- getweight(particles_norm1,gm_normt,K_norm)
    if (sum(is.na(gen.weight.norm))>0){
      ll<-sum(is.na(gen.weight.norm))
      gen.weight.norm[which(is.na(gen.weight.norm))]<-rep(mean(gen.weight.norm,na.rm = T),ll)
    }
    norm.weight1<-weight_norm+gen.weight.norm
    
    #### Getting weights for particles Cancer ####
    
    weight_canc <- rep(0,Nt.canc1)
    K_canc <- getK(xt)
    dim(gm_canct) <- c(NROW(gm_canct),NCOL(gm_canct))
    gen.weight.canc <- getweight(particles_canc1,gm_canct,K_canc)
    if (sum(is.na(gen.weight.canc))>0){
      ll<-sum(is.na(gen.weight.canc))
      gen.weight.canc[which(is.na(gen.weight.canc))]<-rep(mean(gen.weight.canc,na.rm = T),ll)
    }
    canc.weight1<-weight_canc+gen.weight.canc
    
    #### Population and weight control step ####
    
    pop_overall<-particle_cntl(particles_overall1,overall.weight1)
    pop_norm<-particle_cntl(particles_norm1,norm.weight1)
    pop_canc<-particle_cntl(particles_canc1,canc.weight1)
    
    #### Getting new particles and weights for overall ####
    
    particles_overall2<-pop_overall$lamb[,1:2]
    overall.weight2<-pop_overall$logw
    Nt.overall2<-length(overall.weight2)
    nt.overall2<-length(xt)
    tau_overall2<-particles_overall2[,1]
    sig2_overall2<-particles_overall2[,2]
    sig_overall2<-sqrt(sig2_overall2)
    
    #### Getting new particles and weights for Normal ####
    
    particles_norm2<-pop_norm$lamb[,1:2]
    norm.weight2<-pop_norm$logw
    Nt.norm2<-length(norm.weight2)
    nt.norm2<-length(xt)
    tau_norm2<-particles_norm2[,1]
    sig2_norm2<-particles_norm2[,2]
    sig_norm2<-sqrt(sig2_norm2)
    
    #### Getting new particles and weights for Cancer ####
    
    particles_canc2<-pop_canc$lamb[,1:2]
    canc.weight2<-pop_canc$logw
    Nt.canc2<-length(canc.weight2)
    nt.canc2<-length(xt)
    tau_canc2<-particles_canc2[,1]
    sig2_canc2<-particles_canc2[,2]
    sig_canc2<-sqrt(sig2_canc2)
    
    #### Getting fitted values using projected particles Overall ####
    
    overall.project.mean.fit<-matrix(rep(0,nt.overall2*Nt.overall2),ncol=nt.overall2)
    alp  <- tau_overall2*sig2_overall2
    I <- diag(rep(1,nt.overall2))
    for (k in 1:Nt.overall2){
      Aalp <- solve(I+c(alp[k])*K_overall)
      Halp <- try(chol(Aalp),silent = T)
      if(class(Halp)=="try-error"){
        iAalp <- I + c(alp[k])*K_overall
        s <- svd(iAalp)
        D <- solve(diag(s$d))
        Aalp <- s$v%*%D%*%t(s$u)
        Halp <- s$v%*%sqrt(D)
      }
      Z <- rnorm(nt, 0, 1)
      overall.project.mean.fit[k,]<-(Aalp%*%gm_allt)+c(sig_overall2[k])*(Halp%*%Z)
    }
    sig2.particle.overall2 <- particles_overall2[,2]
    tau.particle.overall2<- particles_overall2[,1]
    
    #### Getting fitted values using projected particles Normal ####
    
    norm.project.mean.fit<- matrix(rep(0,nt.norm2*Nt.norm2),ncol=nt.norm2)
    alp  <- tau_norm2*sig2_norm2
    I <- diag(rep(1,nt.norm2))
    for (k in 1:Nt.norm2){
      Aalp <- solve(I+c(alp[k])*K_norm)
      Halp <- try(chol(Aalp),silent = T)
      if(class(Halp)=="try-error"){
        iAalp <- I + c(alp[k])*K_norm
        s<- svd(iAalp)
        D<- solve(diag(s$d))
        Aalp <- s$v%*%D%*%t(s$u)
        Halp <- s$v%*%sqrt(D)
      }
      Z <- rnorm(nt, 0, 1)
      norm.project.mean.fit[k,]<-(Aalp%*%gm_normt)+c(sig_norm2[k])*(Halp%*%Z)
    }
    sig2.particle.norm2 <- particles_norm2[,2]
    tau.particle.norm2<- particles_norm2[,1]
    
    #### Getting fitted values using projected particles Cancer ####
    
    canc.project.mean.fit<- matrix(rep(0,nt.canc2*Nt.canc2),ncol=nt.canc2)
    alp  <- tau_canc2*sig2_canc2
    I <- diag(rep(1,nt.canc2))
    for (k in 1:Nt.canc2){
      Aalp <- solve(I+c(alp[k])*K_canc)
      Halp <- try(chol(Aalp),silent = T)
      if(class(Halp)=="try-error"){
        iAalp <- I + c(alp[k])*K_canc
        s<- svd(iAalp)
        D<- solve(diag(s$d))
        Aalp<- s$v%*%D%*%t(s$u)
        Halp<- s$v%*%sqrt(D)
      }
      Z <- rnorm(nt, 0, 1)
      canc.project.mean.fit[k,]<-(Aalp%*%gm_canct)+c(sig_canc2[k])*(Halp%*%Z)
    }
    sig2.particle.canc2 <- particles_canc2[,2]
    tau.particle.canc2<- particles_canc2[,1]
    
    #### Sigma matrix creation Overall ####
    
    sig_over_mat.proj<- matrix(0,Nt.overall2,nt.overall2)
    for(i in 1:Nt.overall2){
      sig_over_mat.proj[i,]=rep(sig2.particle.overall2[i],each=nt.overall2)
    }
    
    #### Sigma matrix creation Normal ####
    
    sig_norm_mat.proj<- matrix(0,Nt.norm2,nt.norm2)
    for(i in 1:Nt.norm2){
      sig_norm_mat.proj[i,]=rep(sig2.particle.norm2[i],each=nt.norm2)
    } 
    
    #### Sigma matrix creation Cancer ####
    
    sig_canc_mat.proj<- matrix(0,Nt.canc2,nt.canc2)
    for(i in 1:Nt.canc2){
      sig_canc_mat.proj[i,]=rep(sig2.particle.canc2[i],each=nt.canc2)
    }
    
    #### Calculating Marginal Likelihood Overall #### 
    
    lik_overall.proj<-matrix(0,Nt.overall2,nt.overall2)
    sum_lik_overall.proj<-rep(0,Nt.overall2)
    for(l in 1:Nt.overall2){
      for(k in 1:nt.overall2){
        m<- overall.project.mean.fit[l,k]
        s<- sqrt(sig_over_mat.proj[l,k])
        lik_overall.proj[l,k]<- log(dnorm(gm_allt[k],m,s))
      }
      sum_lik_overall.proj[l]<-sum(lik_overall.proj[l,])
    }
    marginal_lik_overall.proj<- mean(sum_lik_overall.proj)
    
    #### Calculating Marginal Likelihood Normal ####
    
    lik_norm.proj<-matrix(0,Nt.norm2,nt.norm2)
    sum_lik_norm.proj<-rep(0,Nt.norm2)
    for(l in 1:Nt.norm2){
      for(k in 1:nt.norm2){
        m<- norm.project.mean.fit[l,k]
        s<- sqrt(sig_norm_mat.proj[l,k])
        lik_norm.proj[l,k]<- log(dnorm(gm_normt[k],m,s))
      }
      sum_lik_norm.proj[l]<-sum(lik_norm.proj[l,])
    }
    marginal_lik_norm.proj<- mean(sum_lik_norm.proj)
    
    #### Calculating Marginal Likelihood Cancer ####
    
    lik_canc.proj<-matrix(0,Nt.canc2,nt.canc2)
    sum_lik_canc.proj<-rep(0,Nt.canc2)
    for(l in 1:Nt.canc2){
      for(k in 1:nt.canc2){
        m<- canc.project.mean.fit[l,k]
        s<- sqrt(sig_canc_mat.proj[l,k])
        lik_canc.proj[l,k]<- log(dnorm(gm_canct[k],m,s))
      }
      sum_lik_canc.proj[l]<-sum(lik_canc.proj[l,])
    }
    marginal_lik_canc.proj<- mean(sum_lik_canc.proj)
    
    #### Calculating Log Bayes Factors ####
    
    M1_bayes.proj<-marginal_lik_overall.proj
    M2_bayes.proj<-marginal_lik_canc.proj+ marginal_lik_norm.proj
    logbayesfactor.proj<-M1_bayes.proj-M2_bayes.proj
    
    ptm.dep.end<- proc.time()
    time.taken.dep<- (ptm.end - ptm.start)[3]/60
    
    #### Summarizing outputs ####
    
    tmp<-data.frame(Window=w,
                    M1_bayes=M1_bayes.proj,
                    M2_bayes=M2_bayes.proj,
                    logbayesfactor=logbayesfactor.proj,
                    time.min= time.taken.dep)
    results.dep<-data.frame(rbind(results.dep,tmp))
  }
  
  #### Final Output ####
  
  results<-data.frame(rbind(results.indep,results.dep))
  
  #### Calling DMRs and Non-DMRs ####
  
  if(is.null(clust.locations$Truth)){
    bayes.cutoff<- -5
    Pred=ifelse(results$logbayesfactor<bayes.cutoff,'Dmr','Non-Dmr')
    len.agg=data.frame(ID=paste('win',1:nwin,sep=''),
                       Length=dim(datsub)[1])
    output=data.frame(Method="BFDAM2",
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
    output=data.frame(Method="BFDAM2",
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