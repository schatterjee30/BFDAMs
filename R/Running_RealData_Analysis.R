
#####################################
#### Running Real Data analysis #####
#####################################

suppressMessages(suppressWarnings(library(bumphunter)))
suppressMessages(suppressWarnings(library(DMRcate)))
suppressMessages(suppressWarnings(library(limma)))
suppressMessages(suppressWarnings(library(ENmix)))
suppressMessages(suppressWarnings(library(doParallel)))

#### Set Path ####

dir = "D:/research/Biometrics_revision/Real_data_analysis"

#### Loading Necessary Functions ####

source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "helper_functions.R", sep = '/'))
source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "run_BFDAM1.R", sep = '/'))
source(paste(dir, "simulation/Chatterjee_2022BFDA/library", "run_BFDAM2.R", sep = '/'))

#### Load Real Data ####

load(paste(dir, "Data", "Luad_methylation.RData", sep = '/'))

#### Configure Real Data ####

location = luad_data$locations
treat.lab = 'Cancer'
control.lab = 'Normal'
data.raw = as.matrix(luad_data$betaVals)
data.raw = logit2(data.raw)
ncontrols = ncol(data.raw[, grep(control.lab, colnames(data.raw))])
ntreats = ncol(data.raw[, grep(treat.lab, colnames(data.raw))])
labels = factor(c(rep("Normal", ncontrols), rep("Cancer", ntreats)))
design = model.matrix(~ labels)
registerDoParallel(cores = 6)

#### Configuring Locations ####

pos = location$Coordinate_37
chr = as.numeric(location$CHR)

#### Fitting Bumphunter ####

fit = bumphunter(as.matrix(data.raw),
                 design = design,
                 chr = chr,
                 pos = pos,
                 maxGap = 1000,
                 pickCutoff = TRUE,
                 pickCutoffQ = 0.95,
                 nullMethod = "permutation",
                 B = 10,
                 type = "M",
                 verbose = TRUE)

#### Organizing Fit Results ####

fit_df = fit$table
fit_df$dmr.pval = fit_df$p.value
fit_df = fit_df[fit_df$dmr.pval < 0.05,]

#### save results ####

save(fit_df, file = paste(dir, "Results", "Bumphunter_Luad.RData", sep = '/'))

#### Fitting DMRcate ####

myannotation = cpg.annotate("array", 
                            data.raw, 
                            what = "M",
                            arraytype = "450K",
                            analysis.type = "differential",
                            design = design,
                            coef = 2, 
                            fdr = 0.05)
fit= dmrcate(myannotation, 
             lambda = 1000, 
             C = 2)

extractCoords <- function(x){ 
  xx = x@coord
  coords = sapply(xx, strsplit, '[:-]') 
  coords = as.data.frame(do.call(rbind, coords), stringsAsFactors=F) 
  colnames(coords) = c('chrom', 'chromStart', 'chromEnd') 
  rownames(coords) = NULL
  coords$Pval = x@Stouffer
  return(coords)
}
fit_df = extractCoords(fit)
fit_df = fit_df[fit_df$Pval < 0.05,]

#### save results ####

save(fit_df,
     file="D:/research/Biometrics_revision/simulation/Chatterjee_2022BFDA/DMRcate_RealData_Result.RData")

#### Fitting comb-p ####

fit = lmFit(data.raw, design)
fit = eBayes(fit)
Pvals = data.frame(fit$p.value)
fit_limma_df = data.frame(chr = location$CHR,
                          start = location$Coordinate_37,
                          end = (location$Coordinate_37) + 1,
                          p = Pvals[, 2],
                          probe = rownames(Pvals))
fit_limma_df = fit_limma_df[complete.cases(fit_limma_df), ]

fit = combp(fit_limma_df,
            dist.cutoff = 1000,
            bin.size = 310,
            seed = 0.01,
            region_plot = FALSE,
            mht_plot = FALSE,
            nCores = 1,
            verbose = TRUE)
fit_df = read.csv("D:/research/Biometrics_revision/simulation/Chatterjee_2022BFDA/resu_combp.csv")
fit_df = fit_df[fit_df$sidak < 0.05,]

#### save results ####

save(fit_df,
     file="D:/research/Biometrics_revision/simulation/Chatterjee_2022BFDA/combp_RealData_Result.RData")

#### Fitting BFDAM1 ####

data = luad_data$betaVals
fit_df = fit.BFDAM1(x = data)
save(fit_df,
     file="D:/research/Biometrics_revision/simulation/Chatterjee_2022BFDA/BFDAM1_RealData_Result.RData")

#### Fitting BFDAM2 ####

data = luad_data$betaVals
fit_df = fit.BFDAM2(x = data)
save(fit_df,
     file="D:/research/Biometrics_revision/simulation/Chatterjee_2022BFDA/BFDAM2_RealData_Result.RData")

