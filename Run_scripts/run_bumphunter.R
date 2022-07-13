
#### Loading required Packages ####

if(! require("bumphunter")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("bumphunter")
}
if(! require("DMRcate")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("DMRcate")
}
suppressMessages(suppressWarnings(library(bumphunter)))
suppressMessages(suppressWarnings(library(DMRcate)))

#### Fitting Bumphunter Function ####

fit.bumphunter = function(x){
  
  #### Input Data ####
  
  treat.lab = 'Tumor'
  control.lab = 'Normal'
  ptm = proc.time()
  bvals = x$betavals
  location = x$cpglocationdf
  clust.locations = x$windowinfo
  data.raw = bvals[, !names(bvals) %in% c('windows')]
  #data.raw = as.matrix(logit2(data.raw))
  
  #### Configuring group labels ####
  
  ncontrols = ncol(data.raw[, grep(control.lab, colnames(data.raw))])
  ntreats = ncol(data.raw[, grep(treat.lab, colnames(data.raw))])
  labels = factor(c(rep("Tumor", ncontrols), rep("Normal", ntreats)))
  design = model.matrix(~ labels)
  
  #### Configuring Locations ####
  
  pos = location$pos
  chr = location$chr.num
  
  #### Fitting Bumphunter ####
  
  fit = tryCatch(
    suppressMessages(
      suppressWarnings(
        bumphunter(as.matrix(data.raw),
                   design = design,
                   chr = chr,
                   pos = pos,
                   maxGap = 1000,
                   pickCutoff = TRUE,
                   pickCutoffQ = 0.95,
                   nullMethod = "permutation",
                   B = 10,
                   type = "M",
                   verbose = FALSE))),
    error = function(e1){NULL})
  elapsedtime = proc.time() - ptm
  
  #### Organizing Fit Results ####
  
  if(!is.null(fit)){
    fit_df = fit$table
    fit_df$dmr.pval = fit_df$p.valueArea
    fit_df$dmr.chr = paste0("chr", as.character(fit_df$chr))
    fit_df$dmr.start = fit_df$start
    fit_df$dmr.end = fit_df$end
    fit_df$ncpgs = fit_df$clusterL
    output = list(DMRs = fit_df,
                  method = "Bumphunter",
                  location = location,
                  clust.locations = clust.locations,
                  time.min = as.numeric(elapsedtime[3])/60)
  } else {
    output = NULL
  }
  return(output)
}
