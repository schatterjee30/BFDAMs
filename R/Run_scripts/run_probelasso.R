
#### Loading required Packages ####

if(! require("ChAMP")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("ChAMP")
}
if(! require("DMRcate")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("DMRcate")
}

#### Fitting ProbeLasso Function ####

fit.probelasso = function(x){
  
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
  labels = factor(c(rep("Normal", ncontrols), rep("Cancer", ntreats)))
  
  #### Fitting ProbeLasso ####
  
  fit = tryCatch(
    suppressMessages(
      suppressWarnings(
        quiet(champ.DMR(beta = as.matrix(data.raw),
                        pheno = labels,
                        method = "ProbeLasso",
                        minProbes = 5,
                        adjPvalDmr = 0.05,
                        meanLassoRadius = 375,
                        minDmrSep = 1000,
                        adjPvalProbe = 0.05,
                        PDFplot = FALSE,
                        Rplot = FALSE)))),
    error = function(e1){NULL})
  elapsedtime = proc.time() - ptm
  
  #### Organizing Fit Results ####
  
  if(!is.null(fit)){
    fit_df = fit$ProbeLassoDMR
    fit_df$dmr.pval = fit_df$dmrP
    fit_df$dmr.chr = fit_df$seqnames
    fit_df$dmr.start = fit_df$start
    fit_df$dmr.end = fit_df$end
    fit_df$ncpgs = rep(20, nrow(fit_df))
    output = list(DMRs = fit_df,
                  method = "ProbeLasso",
                  location = location,
                  clust.locations = clust.locations,
                  time.min = as.numeric(elapsedtime[3])/60)
  } else {
    output = NULL
  }
  return(output)
}
