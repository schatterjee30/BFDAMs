
#### Description ####

# beta.df = Data frame of methylation beta values
# chr.pos = Chromosome and position of CpG sites
# nCpGs = Number of CpGs desired to be in one region
# iter = Number of MCMC iterations to be run
# burn = Number of burn-ins
# seed = A random number for reproducibilty of results
# control.group = Label for samples in group 1 in methylation data frame(beta.df)
# case.group = Label for samples in group 2 in methylation data frame(beta.df)

#### Required Libraries ####

suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(minfi)))

#### The BFDAM1 Function ####

fit.BFDAM1 = function(beta.df, 
                      chr.pos, 
                      nCpGs = 100, 
                      iter = 500, 
                      burn = 50, 
                      seed = 1234, 
                      control.group = 'Normal', 
                      case.group = 'Tumor'){
  
  #### Organizing Input Data ####
  
  set.seed(seed)
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
    x$windows<-rep(seq_len(nrow(x)%/%nCpGs+1L),each=nCpGs,len=nrow(x))
    return(x)
  })
  
  #### Running in Serial by Chromosome ####
  
  result.list<-list()
  for(c in 1:length(chrpos.win)){
    run.chrpos<-chrpos.win[[c]]
    
    #### Fitting each Window by Chromosome ####
    
    results<-NULL
    for(w in 1:length(unique(run.chrpos$windows))){
      win.chrpos<-run.chrpos[run.chrpos$windows == w,]
      win.beta<-beta.df[rownames(beta.df) %in% rownames(win.chrpos),]
      print(paste("Running Genomic Window", w, "of Chromosome", unique(win.chrpos$chr), sep = ' '))
      
      #### Log-Transforming Methylation data ####
      
      n<-nrow(win.beta)
      wincpgs<-rownames(win.chrpos)
      wincpgs<-paste("'",as.character(wincpgs),"'",collapse=", ",sep="")
      x1<-c(1:n)
      yall<-data.frame(logit2(win.beta))
      
      #### Separating Normal and Cancer groups data ####
      
      ynorm<-yall[,grep(control.group,names(yall))]
      ycanc<-yall[,grep(case.group,names(yall))]
      
      #### Combined group mean fitting ####
      
      gm_all<-rep(0,n)
      gm_all<- apply(yall,1,meanfn)
      fitg <- ncs(x1,gm_all,burn,iter)
      overall_mean_fit<-fitg$fhat
      overall_sigma_fit<-fitg$sig2
      overall_tau_fit<-fitg$tau
      lambda_overall<-cbind(overall_tau_fit, overall_sigma_fit)
      
      #### Normal group mean fitting ####
      
      gm_norm<-rep(0,n)
      gm_norm<- apply(ynorm,1,meanfn)
      fitgnorm <- ncs(x1,gm_norm,burn,iter)
      norm_mean_fit<-fitgnorm$fhat
      norm_sigma_fit<-fitgnorm$sig2
      norm_tau_fit<-fitg$tau
      lambda_norm<-cbind(norm_tau_fit, norm_sigma_fit)
      
      #### Cancer group mean fitting ####
      
      gm_canc<-rep(0,n)
      gm_canc<- apply(ycanc,1,meanfn)
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
      logbayesfactor<-M1_bayes/M2_bayes
      
      #### Summarizing output ####
      
      tmp<-data.frame(chr=unique(win.chrpos$chr),
                      window=w,
                      pos.start=min(win.chrpos$pos),
                      pos.end=max(win.chrpos$pos),
                      M1_bayes=M1_bayes,
                      M2_bayes=M2_bayes,
                      logbayesfactor=logbayesfactor,
                      CpGs=wincpgs)
      results<-data.frame(rbind(results,tmp))
    }

      #### Calling DMRs ####
      
      bayes.cutoff<-median(results$logbayesfactor)
      status<-ifelse(results$logbayesfactor<bayes.cutoff,'DMR','Non-DMR')
      results.shr<-data.frame(chr=results$chr,
                              window=results$window,
                              pos.start=results$pos.start,
                              pos.end=results$pos.end,
                              status=status,
                              CpGs=results$CpGs)
      result.list[[c]]<-results.shr
  }
  return(result.list)
}

####################################################################################
################################# Utility Functions ################################
####################################################################################

#### Mean Calculator ####

meanfn<-function(x)
{
  x[is.infinite(x)]<-1  
  return(sum(x)/length(x))
}

#### Getting K Matrix function for BFDAM1 and BFDAM2 methods ####

getK <- function(x){
  n <- length(x)
  dx <- diff(x)
  dx2 <- diff(x,lag=2)
  Qt <- matrix(rep(0,(n-2)*n),n-2,n)
  for (k in 1:(n-2)){
    Qt[k,k]   <- 1/dx[k]
    Qt[k,k+1] <- -1/dx[k+1]-1/dx[k]
    Qt[k,k+2] <- 1/dx[k+1]
  }
  R  <- matrix(rep(0,(n-2)*(n-2)),n-2,n-2)
  R[1,1] <- dx[1]
  R[1,2] <- dx[2]
  for (k in 2:(n-3)){
    R[k,k-1] <- dx[k]
    R[k,k]   <- 2*dx2[k]
    R[k,k+1] <- dx[k+1]
  }
  R[n-2,n-3] <- dx[n-2]
  R[n-2,n-2] <- dx[n-1]
  R <- R/6
  invR <- solve(R)
  K <- t(Qt)%*%invR%*%Qt
  return(K)
}

#### Natural Cubic Smoothing Spline Function for BFDAM1 and BFDAM2 methods ####

ncs <- function(x,y,burn,iter)
{
  n <- length(y)
  if (nargs() < 3){
    burn <- 100  
    nrep <- 1000
  }else{
    burn<-burn
    nrep<-iter
  }
  At <- 1
  Bt <- 10
  As <- 1
  Bs <- 1
  fhat <- rep(0,n)
  alp  <- .01
  sig2 <- 1
  sig  <- sqrt(sig2)
  tau  <- alp/sig2
  K <- getK(x)
  I <- diag(rep(1,n))
  niter <- burn + nrep
  isamp <- 0
  fhatall <- matrix(rep(0,n*nrep),nrow=nrep)
  tauall  <- rep(0,nrep)
  sig2all <- rep(0,nrep)
  for (iter in 1:niter){
    RSS <- t(y-fhat)%*%(y-fhat)
    sig2 <- (RSS/2+Bs)/rgamma(1, n/2+As, 1)
    sig <- sqrt(sig2)
    iAalp <- I + c(alp)*K
    s <- svd(iAalp)
    D <- solve(diag(s$d))
    Aalp <- s$v%*%D%*%t(s$u)
    Halp <- s$v%*%sqrt(D)
    Z <- rnorm(n, 0, 1)
    fhat <- (Aalp%*%y) + c(sig)*(Halp%*%Z)
    rough <- t(fhat)%*%K%*%fhat
    tau <- rgamma(1, (n-2)/2+At, 1)/(rough/2+Bt)
    alp <- tau*sig2
    if (iter>burn){
      isamp <- isamp +1
      fhatall[isamp,] <- fhat
      sig2all[isamp] <- sig2
      tauall[isamp] <- tau
    }
    if ((iter%%2000)==0){
      cat(iter,"/",niter," is done.\n")
    }
  }
  ret <- list("fhat"=fhatall, "sig2"=sig2all, "tau"=tauall)
  return(ret)
}
