
##########################################################################
#### Pipeline for Generating Synthetic Data for Permutation Analysis #####
##########################################################################

suppressMessages(suppressWarnings(library(doParallel)))
suppressMessages(suppressWarnings(library(foreach)))
suppressMessages(suppressWarnings(library(matrixStats)))
suppressMessages(suppressWarnings(library(pbapply)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(lumi)))

#### Setting Path ####

dir = "D:/research/Biometrics_revision/simulation/Chatterjee_2022BFDA"

#### Loading Helper Functions ####

source(paste(dir, "library", "helper_functions.R", sep = '/'))
source(paste(dir, "library", "simulatedData_mval_shuffle.R", sep = '/'))

#### Loading LUAD Real Data ####

load(paste(dir, "Data", "Luad_2015_Fleischer.RData", sep = '/'))
betavals = Luad_data$betavals
cpglocs = Luad_data$annotation
rm(Luad_data)

#### Extracting common CpGs ####

bcpgs = data.frame(cpg = rownames(betavals), df = 'betadf')
lcpgs = data.frame(cpg = rownames(cpglocs), df = 'locdf')
comcpgs = merge(bcpgs, lcpgs, by = 'cpg')
betavals = betavals[rownames(betavals) %in% comcpgs$cpg, ]
cpglocs = cpglocs[rownames(cpglocs) %in% comcpgs$cpg, ]
rm(bcpgs, comcpgs, lcpgs)

#### Specify Simulation Parameters ####

nsamples = 18
lenreg = c(100)
numreg = c(50)
effectSize = c(0) 
seed = seq(100, 10000, by = 100)
istrings = apply(expand.grid(nsamples, lenreg, numreg, effectSize), 
                 1, paste0, collapse = "_")
istringlabs = c("nsamples", "lenreg", "numreg", "effectSize")

#### Generating Specific Simulation Settings ####

for(i in 1:length(istrings)){
  message("Scenario ", i, " Running")
  istring = str_split(istrings[i], pattern = "_")[[1]]
  names(istring) = istringlabs
  numtreat = 9
  pdmr = 0
  lenreg = as.numeric(istring["lenreg"])
  numreg = as.numeric(istring["numreg"])
  effectSize = as.numeric(istring["effectSize"])
  simlist = list()
  for(j in 1:length(seed)){
    message("Replicate ", j, " Running")
    set.seed(seed[j])
    
    #### Simulating ####
    
    simlist = simuData.mval.shuffle(betavals = betavals,
                                    cpglocs = cpglocs,
                                    nsamples = nsamples,
                                    numreg = numreg,
                                    lenreg = lenreg,
                                    effectSize = effectSize)
    ndmrs = numreg * pdmr
    lstrings = paste(istrings[i], seed, sep = "_")
    outstring = paste(paste("LUAD", "nSamples", nsamples,
                            "numRegions", numreg,
                            "lenRegions", lenreg,
                            "nDMRs", ndmrs, 
                            "effectSize", effectSize,
                            "seed", seed[j],
                            sep = '_'),
                      ".RData",  sep = '') 
    
    message('Saving Synthetic Data ', istrings[i])
    save(simlist, file = paste(dir, "simulatedData", "shuffleData", outstring, sep = '/'))
  }
}