
###################################
#### List of Helper Functions #####
###################################

quiet = function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

spikeIn.dmrs = function(x){
  xsub1 = x[, grep('Normal', names(x))]
  xsub2 = x[, grep('Tumor', names(x))]
  xsub1 = xsub1 + effectSize
  dat = data.frame(windows = x$windows, xsub2, xsub1)
  return(dat)
}

GenerateWinInfo = function(x){
  infodf = data.frame(windows = unique(x$windows),
                      chr = unique(x$chr),
                      start.position = min(x$pos),
                      end.position = max(x$pos),
                      lenreg = nrow(x),
                      Truth = unique(x$Truth))
  return(infodf)
}

#### Processing outputs ####

processOutputs = function(reslist){
  dmrtab = reslist$DMRs
  if(is.null(dmrtab)){
    dmrtab = NULL
  }
  if(is.null(dmrtab)){
    dmrtab = NULL
  }else{
    dmrtab = dmrtab[dmrtab$ncpgs > min.len & dmrtab$dmr.pval < 0.05, ]
  }
  if(is.null(dmrtab)){
    sumreslist = list(DMRs = NULL,
                      clustLoc = reslist$clust.locations,
                      Method = reslist$method,
                      Time = reslist$time.min)
  }else{
    sumreslist = list(DMRs = dmrtab,
                      clustLoc = reslist$clust.locations,
                      Method = reslist$method,
                      Time = reslist$time.min)
  }
  return(sumreslist)
}

#### Summarizing outputs ####

summarizeOutputs = function(sumreslist){
  
  if(sumreslist$Method == "BFDAM1" | sumreslist$Method == "BFDAM2"){
    truthtab = sumreslist
    truthtab$status = ifelse(truthtab$Truth %in% "DMR" & truthtab$predicted %in% "DMR", 'TP',
                             ifelse(truthtab$Truth %in% "DMR" & truthtab$predicted %in% "Non-DMR", 'FN',
                                    ifelse(truthtab$Truth %in% "Non-DMR" & truthtab$predicted %in% "DMR", 'FP',
                                           ifelse(truthtab$Truth %in% "Non-DMR" & truthtab$predicted == "Non-DMR", 'TN', NA))))
  }else{
    if(!is.null(sumreslist$DMRs)){
      groupA = data.table(windows = sumreslist$clustLoc$windows,
                          chr = sumreslist$clustLoc$chr,
                          start = sumreslist$clustLoc$start.position,
                          end = sumreslist$clustLoc$end.position)
      groupB = data.table(chr = sumreslist$DMRs$dmr.chr,
                          start = sumreslist$DMRs$dmr.start,
                          end = sumreslist$DMRs$dmr.end)
      setkey(groupA, chr, start, end)
      setkey(groupB, chr, start, end)
      over = foverlaps(groupA, groupB, nomatch = 0)
      over2 = data.table(windows = over$windows,
                         chr = over$chr,
                         start = over[, ifelse(start > i.start, start, i.start)],
                         end = over[, ifelse(end < i.end, end, i.end)])
      over2 = over2[!duplicated(over2$windows), ]
      truthtab = data.frame(windows = sumreslist$clustLoc$windows,
                            Truth = sumreslist$clustLoc$Truth)
      truthtab$predicted = ifelse(truthtab$windows %in% over2$windows, 
                                  "DMR", "Non-DMR")
    }else{
      truthtab = data.frame(windows = sumreslist$clustLoc$windows,
                            Truth = sumreslist$clustLoc$Truth)
      truthtab$predicted = rep('Non-DMR', nrow(truthtab))
    }
    truthtab$status = ifelse(truthtab$Truth %in% "DMR" & truthtab$predicted %in% "DMR", 'TP',
                             ifelse(truthtab$Truth %in% "DMR" & truthtab$predicted %in% "Non-DMR", 'FN',
                                    ifelse(truthtab$Truth %in% "Non-DMR" & truthtab$predicted %in% "DMR", 'FP',
                                           ifelse(truthtab$Truth %in% "Non-DMR" & truthtab$predicted == "Non-DMR", 'TN', NA))))
    truthtab$Method = rep(sumreslist$Method, nrow(truthtab))
    truthtab$Time = rep(sumreslist$Time, nrow(truthtab))
  }
  return(truthtab)
}

#### Calculating benchmark metrics ####

benchMethods = function(truthtab, method){
  if (nrow(truthtab) == 0){
    Measures = data.frame(Method = method,
                          Sensitivity = NA,
                          FPR = NA,
                          FDR = NA,
                          FScore = NA,
                          Time = NA)
  }else{
    status = truthtab$status
    FPs = sum(status %in% "FP")
    TPs = sum(status %in% "TP")
    TNs = sum(status %in% "TN")
    FNs = sum(status %in% "FN")
    Method = method
    Time = unique(truthtab$Time)
    Sensitivity = TPs/(TPs + FNs)
    FPR = 1 - (TNs/(TNs + FPs))
    FDR = FPs/(FPs + TPs)
    FScore = 2 * TPs/(2 * TPs + FNs + FPs) 
    Measures = data.frame(Method = Method,
                          Sensitivity = Sensitivity,
                          FPR = FPR,
                          FDR = FDR,
                          FScore = FScore,
                          Time = Time)
  }
  return(Measures)
}

#### Combining benchmark metrics ####

combineMetrics = function(measures){
  Method = unique(measures$Method)
  Power = mean(measures$Sensitivity, na.rm = TRUE)
  FPR = mean(measures$FPR, na.rm = TRUE)
  Replicates = length(rownames(measures))
  lenreg = as.numeric(simuparams[7])
  numreg = as.numeric(simuparams[5])
  effectSize = as.numeric(simuparams[11])
  nsamples = as.numeric(simuparams[3])
  if(simuparams[10] == 'Rho'){
    rho = as.numeric(simuparams[11])
    combination = paste(Method, nsamples, rho, lenreg, numreg, 
                        effectSize, Replicates, sep = '_')
    results = data.frame(Method = Method,
                         Power = Power,
                         FPR = FPR,
                         EffectSize = effectSize,
                         Rho = rho,
                         Len.region = lenreg,
                         Num.region = numreg,
                         Nsamples = nsamples,
                         Replicates = Replicates,
                         Combination = combination)
  }else{
    combination = paste(Method, nsamples, lenreg, numreg, 
                        effectSize, Replicates, sep = '_')
    results = data.frame(Method = Method,
                         Power = Power,
                         FPR = FPR,
                         FDR = FDR,
                         Fscore = FScore,
                         Time.mins = Time.mins,
                         EffectSize = effectSize,
                         Len.region = lenreg,
                         Num.region = numreg,
                         Nsamples = nsamples,
                         Replicates = Replicates,
                         Combination = combination)
  }
  return(results)
}

#### Mean Calculator ####

mean<-function(x)
{
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

ncs <- function(x,y,burn,nrep)
{
  n <- length(y)
  if (nargs() < 3){
    burn <- 1000  
    nrep <- 1000
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
    if ((iter%%1000)==0){
      cat(iter,"/",niter," is done.\n")
    }
  }
  ret <- list("fhat"=fhatall, "sig2"=sig2all, "tau"=tauall)
  return(ret)
}

#### Getting Particle Weights for BFDAM2 method ####

getweight <- function(lamb,y,K){
  N <- nrow(lamb)    
  n.j<-dim(y)[2]
  logw <- rep(0,N)
  logwm<- NULL
  for(j in 1:n.j)
  {
    print(paste("getwt: sample",j,"is running"))
    tau  <- lamb[,(j)]
    sig2 <- lamb[,(2*j)]
    alp<- tau*sig2
    n <- nrow(y)  # number of data points
    I <- diag(rep(1,n))
    logw1 <- rep(0,N)
    for (k in 1:N){
      Aalp <- solve(I+alp[k]*K)
      Halp <- try(chol(Aalp),silent = T)
      if(class(Halp)=="try-error"){
        iAalp <- I + c(alp[k])*K
        s <- svd(iAalp)
        D <- solve(diag(s$d))
        Aalp <- s$v%*%D%*%t(s$u)
        Halp <- s$v%*%sqrt(D)
      }
      logw1[k] <- (n-2)/2*log(tau[k]) + 0.5*sum(log(diag(Halp))) - 1/(2*sig2[k]) * t(y[,j])%*%(I-Aalp)%*%y[,j]
    }
    logwm<- cbind(logwm,logw1)
  }
  logw<- apply(logwm,1,sum)
  return(logw)
}

#### Population control function for BFDAM2 method ####

particle_cntl <- function(lamb,logw){
  N <- nrow(lamb)
  n.j<-dim(lamb)[2]/2
  Wlow <- -5   
  Wup  <- 5
  Nmin <- 10000
  Nmax <- 30000
  Nlow <- 15000
  Nup  <- 25000
  lambda <- 2
  kappa  <- 10
  zi <- 0
  c1 <- 1
  c2 <- 1
  c3 <- 1
  if (N<Nlow){
    Wlow <- Wlow - log(lambda)
    Wup  <- Wup  - log(lambda)
  }else if(N>Nup){
    Wlow <- Wlow + log(lambda)
    Wup  <- Wup  + log(lambda)
  }
  cont_loop <- 1
  while(cont_loop){
    prune_flag  <- rep(0,N)
    enrich_flag <- rep(1,N)
    W <- rep(0,N)
    for(samp in 1:N){
      W[samp] <- logw[samp]
      if(W[samp] < Wlow){          
        if(runif(1) < (1-exp(W[samp]-Wlow))){
          prune_flag[samp] <- 1
        }else{
          W[samp] <- Wlow
        }
      }else if(W[samp] > Wup){    
        d <- ceiling(exp(W[samp]-Wup)+1)
        W[samp] <- W[samp] - log(d)
        enrich_flag[samp] <- d
      }
    }
    Nnew <- sum(enrich_flag[prune_flag==0])
    if(Nnew > Nmax){
      Wlow <- Wlow + log(lambda)
      Wup <- Wup + log(lambda)
      print(paste(Nnew,'is too large.'))
    }else if(Nnew < Nmin){
      Wlow <- Wlow-log(lambda)
      Wup <- Wup -log(lambda)
    }else{
      cont_loop <- 0
    }
  }
  lambnew<-NULL
  for(j in 1:n.j)
  {
    lambnew <- cbind(lambnew,rep(0,Nnew),rep(0,Nnew))
  }
  logwnew <- cbind(rep(0,Nnew))
  sampnew <- 1
  for(samp in 1:N){
    if(prune_flag[samp]==0){
      for(rep in 1:enrich_flag[samp]){
        lambnew[sampnew,] <- lamb[samp,]
        logwnew[sampnew] <- W[samp]
        sampnew <- sampnew+1
      }
    }
  }
  ret <- list("lamb"=lambnew,"logw"=logwnew)
  return(ret)
}