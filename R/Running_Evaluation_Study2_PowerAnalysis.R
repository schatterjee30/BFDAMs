
###########################################
#### Running Simulation Data analysis #####
###########################################

suppressMessages(suppressWarnings(library(bumphunter)))
suppressMessages(suppressWarnings(library(DMRcate)))
suppressMessages(suppressWarnings(library(ENmix)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(pbapply)))
suppressMessages(suppressWarnings(library(Repitools)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(lumi)))

#### Setting Path ####

dir = "D:/research/Biometrics_revision/simulation/Chatterjee_2022BFDA"  

#### Loading Necessary Functions ####

source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "helper_functions.R", sep = '/'))
source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "run_bumphunter.R", sep = '/'))
source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "run_dmrcate.R", sep = '/'))
source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "run_combp.R", sep = '/'))
source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "run_BFDAM1.R", sep = '/'))
source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "run_BFDAM2.R", sep = '/'))

#### Loading Required Files ####

files = list.files(paste(dir, "simulatedData", sep = '/'), "*.RData")

#### Sequentially Performing Benchmarking ####

all.results = NULL
for(i in 1:length(files)){
  message("Processing Synthetic Data ", i)
  load(paste(dir, "simulatedData", files[i], sep = '/'))
  simuparams = str_split(files[i], pattern = "_")[[1]]
  lenreg = as.numeric(simuparams[7])
  min.len = lenreg * 0.20
  
  #### Running All Methods ####
  
  message("Running Bumphunter")
  bump.res = pblapply(simlist, function(x) fit.bumphunter(x))
  
  message("Running DMRCate")
  dmrcat.res = pblapply(simlist, function(x) fit.dmrcate(x))
  
  message("Running Combp")
  combp.res = pblapply(simlist, function(x) fit.combp(x))
  
  message("Running bfdam1")
  bfdam1.res = pblapply(simlist, function(x) fit.BFDAM1(x, iter = 20000, burn = 1000))
  
  message("Running bfdam2")
  bfdam2.res = pblapply(simlist, function(x) fit.BFDAM2(x, iter = 20000, burn = 1000))
  
  #### Processing Results from Existing Methods ####
  
  message("Processing Bumphunter Results")
  bump.proccess = pblapply(bump.res, function(x) processOutputs(x))
  
  message("Processing DMRcate Results")
  dmrcat.proccess = pblapply(dmrcat.res, function(x) processOutputs(x))
  
  message("Processing Combp Results")
  combp.proccess = pblapply(combp.res, function(x) processOutputs(x))
  
  #### Summarizing Results ####
  
  message("Summarizing Bumphunter Results")
  bump.summary = pblapply(bump.proccess, function(x) summarizeOutputs(x))
  
  message("Summarizing DMRcate Results")
  dmrcat.summary = pblapply(dmrcat.proccess, function(x) summarizeOutputs(x))
  
  message("Summarizing Combp Results")
  combp.summary = pblapply(combp.proccess, function(x) summarizeOutputs(x))
  
  message("Summarizing BFDAM1 Results")
  bfdam1.summary = pblapply(bfdam2.res, function(x) summarizeOutputs(x))
  
  message("Summarizing BFDAM2 Results")
  bfdam2.summary = pblapply(bfdam2.res, function(x) summarizeOutputs(x))
  
  #### Benchmarking Metrics ####
  
  message("Benchmarking Bumphunter Results")
  bump.bench = do.call('rbind', pblapply(bump.summary, function(x) benchMethods(x, method = 'Bumphunter')))
  bump.bench = combineMetrics(bump.bench)
  
  message("Benchmarking DMRcate Results")
  dmrcat.bench = do.call('rbind', pblapply(dmrcat.summary, function(x) benchMethods(x, method = 'DMRcate')))
  dmrcat.bench = combineMetrics(dmrcat.bench)
  
  message("Benchmarking Combp Results")
  combp.bench = do.call('rbind', pblapply(combp.summary, function(x) benchMethods(x, method = 'Combp')))
  combp.bench = combineMetrics(combp.bench)
  
  message("Benchmarking BFDAM1 Results")
  bfdam1.bench = do.call('rbind', pblapply(bfdam1.summary, function(x) benchMethods(x, method = 'BFDAM1')))
  bfdam1.bench = combineMetrics(bfdam1.bench)
  
  message("Benchmarking BFDAM2 Results")
  bfdam2.bench = do.call('rbind', pblapply(bfdam2.summary, function(x) benchMethods(x, method = 'BFDAM2')))
  bfdam2.bench = combineMetrics(bfdam2.bench)
  
  #### Appending Results ####
  
  results = data.frame(rbind(bump.bench, dmrcat.bench, combp.bench, 
                             bfdam1.bench, bfdam2.bench))
  all.results = data.frame(rbind(all.results, results))
}

#### Saving All Results ####

save(all.results, file = paste(dir, 'outputs', 
                               paste("AllResults_CompetingMethods", '.RData', sep = ''), sep = '/'))

