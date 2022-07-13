
#######################################################################
#### Defining Simulation Function for Dependent Methylated regions ####
#######################################################################

simuData.mval.shuffle = function(betavals,
                                 cpglocs,
                                 nsamples = 18,
                                 numreg = 50,
                                 lenreg = 100,
                                 effectSize = 0){
  
  #### Extracting Normal Samples ####
  
  betavalsub = betavals[, grep("Normal", names(betavals))]
  betavalsub = betavalsub[, 1:nsamples]
  
  #### Selecting CpGs for simulating purposes ####
  
  ncpgs = numreg * lenreg
  betavalsub = betavalsub[1:ncpgs, ]
  betavalsub = betavalsub[, sample(names(betavalsub), nsamples)]
  cpglocsub = cpglocs[rownames(cpglocs) %in% rownames(betavalsub), ]
  names(betavalsub) = c(paste("Tumor_", 1:(nsamples/2), sep = ''),
                        paste("Normal_", 1:(nsamples/2), sep = ''))
  
  
  #### Creating Regions for simulation CpGs ####
  
  windows = rep(seq_len(nrow(cpglocsub)%/% lenreg + 1L), each = lenreg, len = nrow(cpglocsub))
  cpglocsub = data.frame(windows, cpglocsub)
  cpglocsub$Truth = "Non-DMR"
  betavalsub = data.frame(windows = windows, betavalsub)
  betavalsub = betavalsub[order(betavalsub$windows), ]
  
  #### Gathering Window information ####
  
  wininfolist = split(cpglocsub, betavalsub$windows)
  wininfo = do.call('rbind', pblapply(wininfolist, 
                                      function(x) GenerateWinInfo(x)))
  
  #### Saving Synthetic Data ####
  
  simdata = list(betavals = betavalsub,
                 cpglocationdf = cpglocsub,
                 windowinfo = wininfo)
  return(simdata)
}