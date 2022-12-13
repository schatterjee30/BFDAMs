
###########################################################################################
#### 450k Methylation Data Simulator Function for Power Analysis (Simulation Study 2) #####
###########################################################################################

# betavals = Methylation 450K array data
# cpglocs = Spatial location of Methylation sites
# grplab = Group Label from which the data will be subsetted for downstream simulation
# nsamples = Number of Samples in the array Data
# numreg = Number of Regions to be simulated
# lenreg = Length of regions to be simulated
# pdmr = proportion of differentially methylated regions
# effectSize = Two Group mean difference in mvalue scale

simuData.mval = function(betavals,
                         cpglocs,
                         grplab = 'Normal',
                         nsamples = 18,
                         numreg = 50,
                         lenreg = 100,
                         pdmr = 0.3,
                         effectSize = 2){
  
  #### Extracting Normal Samples ####
  
  betavalsub = betavals[, grep(grplab, names(betavals))]
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
  
  mvaldmr = beta2m(betvaldmr)
  mvaldmr = data.frame(windows = cpglocdmr$windows, mvaldmr)
  mvaldmrlist = split(mvaldmr, mvaldmr$windows)
  mvaldmr.spiked = do.call('rbind', pblapply(mvaldmrlist, function(x) spikeIn.dmrs(x)))
  rownames(mvaldmr.spiked) = rownames(mvaldmr)
  mvaldmr.spiked = mvaldmr.spiked[,-1]
  syndmrs = m2beta(mvaldmr.spiked)
  syndmrs[syndmrs > 0.99] = 0.99
  
  #### Organizing DMRs & Non-DMRs ####
  
  syndmrs = data.frame(windows = cpglocdmr$windows, syndmrs)
  synnondmrs = data.frame(windows = cpglocnondmr$windows, betvalnondmr)
  
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