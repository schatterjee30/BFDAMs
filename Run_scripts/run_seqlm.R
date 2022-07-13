

#### Loading required Packages ####

if(! require("seqlm")) {
  if(! require("devtools")) {
    install.packages("devtools")
  }
  install_github("raivokolde/seqlm")
}
if(! require("Repitools")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("Repitools")
}

#### Fitting seqlm Function ####

fit.seqlm = function(x, genome_information){
  
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
  
  #### Fitting seqlm ####
  
  fit = tryCatch(
    suppressMessages(
      suppressWarnings(
        quiet(seqlm(values = as.matrix(data.raw), 
                    genome_information = genome_information, 
                    annotation = labels)))),
    error = function(e1){NULL})
  elapsedtime = proc.time() - ptm
  
  #### Organizing Fit Results ####
  
  if(!is.null(fit)){
    fit_df = annoGR2DF(fit)
    fit_df$dmr.pval = fit_df$fdr
    fit_df$dmr.chr = paste0("chr", as.character(fit_df$chr))
    fit_df$dmr.start = fit_df$start
    fit_df$dmr.end = fit_df$end
    ncpgs = str_split(fit_df$probes, pattern = ";")
    ncpgs = do.call('c', lapply(ncpgs, length))
    fit_df$ncpgs = ncpgs
    output = list(DMRs = fit_df,
                  method = "seqlm",
                  location = location,
                  clust.locations = clust.locations,
                  time.min = as.numeric(elapsedtime[3])/60)
  } else {
    output = NULL
  }
  return(output)
}
