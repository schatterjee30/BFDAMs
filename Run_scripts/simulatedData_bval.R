
#######################################################################
#### Defining Simulation Function for Dependent Methylated regions ####
#######################################################################

simuData.bval = function(betavals,
                         cpglocs,
                         nsamples = 18,
                         numreg = 50,
                         lenreg = 100,
                         pdmr = 0.2,
                         effectSize = 0.4){
  
  #### Extracting Normal Samples ####
  
  betavalsub = betavals[, grep("Normal", names(betavals))]
  betavalsub = betavalsub[, sample(names(betavalsub), nsamples)]
  
  #### Selecting CpGs for simulating purposes ####
  
  ncpgs = numreg * lenreg
  betavalsub = betavalsub[1:ncpgs, ]
  cpglocsub = cpglocs[rownames(cpglocs) %in% rownames(betavalsub), ]
  names(betavalsub) = c(paste("Tumor_", 1:(nsamples/2), sep = ''),
                        paste("Normal_", 1:(nsamples/2), sep = ''))
  
  #### Creating Regions for simulation CpGs ####
  
  windows = rep(seq_len(nrow(cpglocsub)%/% lenreg + 1L), each = lenreg, len = nrow(cpglocsub))
  cpglocsub = data.frame(windows, cpglocsub)
  
  #### Selecting Regions with their Locations to make DMRs ####
  
  ndmr = pdmr * numreg
  dmrids = sample(unique(cpglocsub$windows), ndmr)
  cpglocdmr = cpglocsub[cpglocsub$windows %in% dmrids, ]
  cpglocdmr$Truth = rep("DMR", nrow(cpglocdmr))
  betvaldmr = betavalsub[rownames(betavalsub) %in% rownames(cpglocdmr), ]
  
  #### Selecting Regions with their Locations to make Non-DMRs ####
  
  nondmrids = unique(cpglocsub$windows)[!unique(cpglocsub$windows) %in% dmrids]
  cpglocnondmr = cpglocsub[cpglocsub$windows %in% nondmrids, ]
  cpglocnondmr$Truth = rep("Non-DMR", nrow(cpglocnondmr))
  betvalnondmr = betavalsub[rownames(betavalsub) %in% rownames(cpglocnondmr), ]
  rm(betavals, cpglocs, betavalsub, cpglocsub)
  
  #### Creating DMRs & Non-DMRs ####
  
  bvaldmr = data.frame(windows = cpglocdmr$windows, betvaldmr)
  bvaldmrlist = split(bvaldmr, bvaldmr$windows)
  bvaldmr.spiked = do.call('rbind', pblapply(bvaldmrlist, function(x) spikeIn.dmrs(x)))
  rownames(bvaldmr.spiked) = rownames(bvaldmr)
  syndmrs = bvaldmr.spiked
  syndmrs[syndmrs > 0.99] = 0.99
  synnondmrs = betvalnondmr
  
  #### Organizing DMRs & Non-DMRs ####
  
  syndmrs = data.frame(windows = cpglocdmr$windows, syndmrs)
  synnondmrs = data.frame(windows = cpglocnondmr$windows, synnondmrs)
  
  #### Organizing output files ####
  
  synthetic.betas = data.frame(rbind(synnondmrs, syndmrs))
  synthetic.betas = synthetic.betas[order(synthetic.betas$windows), ]
  synthetic.locs = data.frame(rbind(cpglocdmr, cpglocnondmr))
  synthetic.locs = synthetic.locs[order(synthetic.locs$windows), ]
  
  #### Gathering Window information ####
  
  wininfolist = split(synthetic.locs, synthetic.locs$windows)
  wininfo = do.call('rbind', pblapply(wininfolist, 
                                      function(x) GenerateWinInfo(x)))
  
  #### Saving Synthetic Data ####
  
  simdata = list(betavals = synthetic.betas,
                 cpglocationdf = synthetic.locs,
                 windowinfo = wininfo)
  return(simdata)
}

