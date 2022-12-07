
## This code differential analysis between 2 groups to detect Differentially Methylated Regions (DMRs) using the DMRcate Method 
## The funcion fit.bumphunter performs the DMR analysis and returns a list object which contains a data frame of DMRs among others 
## The only input to this function is a list which has 3 data frames (DF). 
## The 1st DF (betavals) should contain the Methylation data where rows are CpG sites and the columns are the window information of the respective CpGs along with normal or tumor samples. 
## The 2nd DF (cpglocationdf) should contain the annotations of the CpGs obtained using the Illumina manifest file. 
## The 3rd DF (windowinfo) should contain the information about the span (region start - end) of each genomic window created by the user.

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


