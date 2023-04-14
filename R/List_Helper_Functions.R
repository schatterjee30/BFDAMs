
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