
## This code differential analysis between 2 groups to detect Differentially Methylated Regions (DMRs) using the Comb-P Method ##
## The funcion fit.bumphunter performs the DMR analysis and returns a list object which contains a data frame of DMRs among others ##
## The only input to this function is a Methylation data where rows are CpG sites and columns are normal or tumor samples ##

#### Loading required Packages ####

if(! require("ENmix")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("ENmix")
}
if(! require("limma")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("limma")
}
if(! require("DMRcate")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("DMRcate")
}

#### Fitting Comb-p Function ####

fit.combp = function(x){
  
  #### Input Data ####
  
  dir = "/Path/to where/you want/the results to be stored"  
  setwd(paste(dir, "outputs", sep = '/'))
  treat.lab = 'Tumor'
  control.lab = 'Normal'
  ptm = proc.time()
  bvals = x$betavals
  location = x$cpglocationdf
  clust.locations = x$windowinfo
  data.raw = as.matrix(bvals[, !names(bvals) %in% c('windows')])
  
  #### Configuring group labels ####
  
  ncontrols = ncol(data.raw[, grep(control.lab, colnames(data.raw))])
  ntreats = ncol(data.raw[, grep(treat.lab, colnames(data.raw))])
  labels = factor(c(rep("Tumor", ncontrols), rep("Normal", ntreats)))
  design = model.matrix(~ labels)
  
  #### Fitting Limma to get site wise P-values ####
  
  fit = lmFit(data.raw, design)
  fit = eBayes(fit)
  Pvals = data.frame(fit$p.value)
  fit_limma_df = data.frame(chr = location$chr,
                            start = location$pos,
                            end = (location$pos) + 1,
                            p = Pvals[, 2],
                            probe = rownames(Pvals))
  fit_limma_df = fit_limma_df[complete.cases(fit_limma_df), ]
  
  #### Fitting combp ####
  
  quiet(combp(fit_limma_df,
              dist.cutoff = 1000,
              bin.size = 310,
              seed = 0.01,
              region_plot = FALSE,
              mht_plot = FALSE,
              nCores = 1,
              verbose = FALSE))
  elapsedtime = proc.time() - ptm
  
  #### Reading in Output ####
  
  tmp = read.csv(paste("/Path/to where/you stored the results", 'resu_combp.csv', sep = '/'))
  if(!is.null(fit)){
    tmp$dmr.pval = tmp$sidak
    tmp$dmr.chr = tmp$chr
    tmp$dmr.start = tmp$start
    tmp$dmr.end = tmp$end
    tmp$ncpgs = tmp$nprobe
    output = list(DMRs = tmp,
                  method = "combp",
                  location = location,
                  clust.locations = clust.locations,
                  time.min = as.numeric(elapsedtime[3])/60)
  }else{
    output = NULL
  }
  return(output)
}



