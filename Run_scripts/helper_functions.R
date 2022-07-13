
#################################################
#### Functions for Generating Synthetic Data ####
#################################################

#### Assigning Region numbers ####

spikeIn.dmrs = function(x){
  grp1 =  x[, grep('Normal', names(x))]
  grp2 = grp1 + effectSize
  names(grp1) = paste('Tumor', 1:ncol(grp1), sep = '_')
  syn = data.frame(cbind(grp1, grp2))
  return(syn)
}


#### Assigning Region numbers ####

createRegions = function(betavals, cpglocs){
  chrlist = split(cpglocs, cpglocs$chr.num)
  windf = do.call('rbind', pblapply(chrlist, function(x){
    x$windows = paste(x$chr, rep(seq_len(nrow(x)%/% lenreg + 1L), 
                                 each = lenreg, len = nrow(x)), sep = '_')
    return(x)
  }))
  rownames(windf) = gsub("^.*\\.","", rownames(windf))
  return(windf)
}

#### Generate Correlated Beta Random Variates ####

rmbeta = function(n, a, b, Sigma) {
  require(MASS)
  Z = mvrnorm(n, rep(0, nrow(Sigma)), Sigma)
  t(qbeta(t(pnorm(Z, 0, sqrt(diag(Sigma)))), a, b))
}

#### Estimate Shape parameters of Beta Distribution ####

estBetaParams = function(mu, var) {
  alpha = (((1 - mu) / var) - (1 / mu)) * (mu ^ 2)
  beta = alpha * ((1 / mu) - 1)
  vals = data.frame(alpha = alpha, beta = beta)
  return(vals)
}

#### Create DMRs with Rho parameter ####

GenerateDmrs_wRho = function(x){
  xsub = x[, !names(x) %in% 'windows']
  k = nrow(xsub)
  Sigma = rho * matrix(1, k, k) + (1 - rho) * diag(k)
  params = estBetaParams(mu = as.numeric(rowMeans(xsub)),
                         var = as.numeric(rowVars(as.matrix(xsub))))
  betas = data.frame(t(rmbeta(n = nsamples, 
                              a = as.numeric(params$alpha),
                              b = as.numeric(params$beta),
                              Sigma = Sigma)))
  names(betas) = names(xsub)
  spiked.betas = apply(betas, 1, function(x){
    xcanc = x[1:9]
    xnorm = x[10:18]
    mu1 = mean(xcanc, na.rm = T)
    mu2 = mean(xnorm, na.rm = T)
    if(mu1 > mu2){
      xnorm = xnorm + effectSize
    }else{
      xcanc = xcanc + effectSize
    }
    beta.vec = c(xcanc, xnorm)
    return(beta.vec)
  })
  spiked.betas = data.frame(t(spiked.betas))
  spiked.betas[spiked.betas > 0.99] = 0.999
  rownames(spiked.betas) = rownames(xsub)
  return(spiked.betas)
}

#### Create Non-DMRs with Rho parameter ####

GenerateNonDmrs_wRho = function(x){
  xsub = x[, !names(x) %in% 'windows']
  k = nrow(xsub)
  Sigma = rho * matrix(1, k, k) + (1 - rho) * diag(k)
  params = estBetaParams(mu = as.numeric(rowMeans(xsub)),
                         var = as.numeric(rowVars(as.matrix(xsub))))
  betas = data.frame(t(rmbeta(n = nsamples, 
                              a = as.numeric(params$alpha),
                              b = as.numeric(params$beta),
                              Sigma = Sigma)))
  names(betas) = names(xsub)
  betas[betas > 0.99] = 0.999
  rownames(betas) = rownames(xsub)
  return(betas)
}

#### Create DMRs with no Rho parameter ####

GenerateDmrs_woRho = function(x){
  xcanc = x[1:9]
  xnorm = x[10:18]
  xcanc.m = mean(xcanc, na.rm = T)
  xnorm.m = mean(xnorm, na.rm = T)
  if(xcanc.m > xnorm.m){
    xnorm = xnorm + effectSize
  }else{
    xcanc = xcanc + effectSize
  }
  b = c(xcanc, xnorm)
  return(b)
}

#### Replace NAN beta values with median beta values ####

nan.rm = function(x){
  x = as.numeric(x)
  if(length(which(is.nan(x))) >= 1){
    xsub = x[!is.nan(x)]
    x[is.nan(x)] = median(xsub, na.rm = TRUE)
  }
  return(x)
}

#### Generating Windows Information ####

GenerateWinInfo = function(x){
  infodf = data.frame(windows = unique(x$windows),
                      chr = unique(x$chr),
                      start.position = min(x$pos),
                      end.position = max(x$pos),
                      Truth = unique(x$Truth))
  return(infodf)
}

################################################
### Functions for Comparing Existing Methods ###
################################################

#### Function to suppress any messages from other methods ####

quiet = function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

#### Processing outputs from competing methods ####

processOutputs = function(reslist){
  dmrtab = reslist$DMRs
  
  #### Checking NULL DMRs Output ####
  
  if(is.null(dmrtab)){
    dmrtab = NULL
  }
  
  #### Checking DMRs with length > 5 CpGs and significant at 5% ####
  
  if(is.null(dmrtab)){
    dmrtab = NULL
    }else{
      dmrtab = dmrtab[dmrtab$ncpgs > min.len & dmrtab$dmr.pval < 0.05, ]
    }
  
  #### Extracting cluster locations and annotations for DMRs ####
  
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

#### Summarizing outputs from competing methods ####

summarizeOutputs = function(sumreslist){
  
  if(!is.null(sumreslist$DMRs)){
    
    #### Creating Ranges for True and Predicted DMRs ####
    
    groupA = data.table(windows = sumreslist$clustLoc$windows,
                        chr = sumreslist$clustLoc$chr,
                        start = sumreslist$clustLoc$start.position,
                        end = sumreslist$clustLoc$end.position)
    groupB = data.table(chr = sumreslist$DMRs$dmr.chr,
                        start = sumreslist$DMRs$dmr.start,
                        end = sumreslist$DMRs$dmr.end)
    setkey(groupA, chr, start, end)
    setkey(groupB, chr, start, end)
    
    #### Finding Overlaps with True DMRs ####
    
    over = foverlaps(groupA, groupB, nomatch = 0)
    over2 = data.table(windows = over$windows,
                       chr = over$chr,
                       start = over[, ifelse(start > i.start, start, i.start)],
                       end = over[, ifelse(end < i.end, end, i.end)])
    over2 = over2[!duplicated(over2$windows), ]
    
    #### Creating Table with True and Predicted DMRs ####
    
    truthtab = data.frame(windows = sumreslist$clustLoc$windows,
                          Truth = sumreslist$clustLoc$Truth)
    truthtab$predicted = ifelse(truthtab$windows %in% over2$windows, 
                                "DMR", "Non-DMR")
  }else{
    truthtab = data.frame(windows = sumreslist$clustLoc$windows,
                          Truth = sumreslist$clustLoc$Truth)
    truthtab$predicted = rep('Non-DMR', nrow(truthtab))
  }
  
  #### Creating Status for Benchmarking purposes ####
  
  truthtab$status = ifelse(truthtab$Truth %in% "DMR" & truthtab$predicted %in% "DMR", 'TP',
                           ifelse(truthtab$Truth %in% "DMR" & truthtab$predicted %in% "Non-DMR", 'FN',
                                  ifelse(truthtab$Truth %in% "Non-DMR" & truthtab$predicted %in% "DMR", 'FP',
                                         ifelse(truthtab$Truth %in% "Non-DMR" & truthtab$predicted == "Non-DMR", 'TN', NA))))
  truthtab$Method = rep(sumreslist$Method, nrow(truthtab))
  truthtab$Time = rep(sumreslist$Time, nrow(truthtab))
  return(truthtab)
}

#### Calculating benchmark metrics for competing methods ####

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
    
    #### Getting Basic Measures ####
    
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

#### Combining benchmark metrics for competing methods ####

combineMetrics = function(measures){
  Method = unique(measures$Method)
  Power = mean(measures$Sensitivity, na.rm = TRUE)
  FPR = mean(measures$FPR, na.rm = TRUE)
  FDR = mean(measures$FDR, na.rm = TRUE)
  FScore = mean(measures$FScore, na.rm = TRUE)
  Time.mins = mean(measures$Time, na.rm = TRUE)
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
                         FDR = FDR,
                         Fscore = FScore,
                         Time.mins = Time.mins,
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
