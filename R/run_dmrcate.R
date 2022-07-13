
#### Loading required Packages ####

if(! require("DMRcate")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("DMRcate")
}
if(! require("IlluminaHumanMethylation450kanno.ilmn12.hg19")) {
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
}

#### Fitting DMRcate Function ####

fit.dmrcate = function(x){
  
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
  
  #### Annotating Probes ####
  
  myannotation = suppressMessages(cpg.annotate("array", 
                                               as.matrix(data.raw), 
                                               what = "M",
                                               arraytype = "450K",
                                               analysis.type = "differential",
                                               design = design,
                                               coef = 2, 
                                               fdr = 0.05))
  
  #### Fitting DMRcate ####
  
  fit = tryCatch(
    suppressMessages(
      suppressWarnings(
        quiet(dmrcate(myannotation, 
                      lambda = 1000, 
                      C = 3)))),
    error = function(e1){NULL})
  elapsedtime = proc.time() - ptm
  
  #### Organizing Fit Results ####
  
  if(!is.null(fit)){
    fit_df = suppressMessages(suppressWarnings(data.frame(extractRanges(fit, genome = "hg19"))))
    fit_df$dmr.pval = fit_df$Stouffer
    fit_df$dmr.chr = fit_df$seqnames
    fit_df$dmr.start = fit_df$start
    fit_df$dmr.end = fit_df$end
    fit_df$ncpgs = fit_df$no.cpgs
    output = list(DMRs = fit_df,
                  method = "DMRcate",
                  location = location,
                  clust.locations = clust.locations,
                  time.min = as.numeric(elapsedtime[3])/60)
  } else {
    output = NULL
  }
  return(output)
}

