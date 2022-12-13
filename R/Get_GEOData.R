
#########################################################################
#### Downloading & Pre-processing TCGA Lung Adenocarcinoma Real Data ####
#########################################################################

Get_GEOData = function(workDirectory, geoID = GSE66836, outfile_name = Luad_2015_Fleischer){
  library(GEOquery)
  library(minfi)
  library(stringr)
  library(data.table)
  
  #### Set Path ####
  
  dir = workDirectory
  setwd(paste(dir, "Data", sep = '/'))
  
  #### Download Methylation Files ####
  
  geoID = geoID
  if(dir.exists(paste(dir, "Data", geoID, sep = '/'))){
    stop("Data with Specified GEO ID exists")
  }else{
    getGEOSuppFiles(geoID)
    untar(paste(geoID, paste(geoID, "RAW.tar", sep = '_'), sep = '/'), 
          exdir = paste(geoID, "idat", sep = '/'))
  }
  
  #### Decompressing Methylation Files ####
  
  idatFiles = list.files(paste(geoID, "idat", sep = '/'), 
                         pattern = "idat.gz$", full = TRUE)
  sapply(idatFiles, gunzip, overwrite = TRUE)
  
  #### Reading idat files ####
  
  rgSet = read.metharray.exp(paste(geoID, "idat", sep = '/'))
  
  #### Normalizing Data ####
  
  grSet = preprocessFunnorm(rgSet)
  
  #### Getting CpG Location Data ####
  
  annotation = data.frame(getAnnotation(grSet))
  
  #### Remove CpGs Near SNPs ####
  
  grSet = addSnpInfo(grSet)
  grSet = dropLociWithSnps(grSet, snps=c("SBE","CpG"), maf = 0)
  
  #### Extracting Beta values ####
  
  betavals = data.frame(getBeta(grSet))
  
  #### Extracting all Autosomes ####
  
  annotation = annotation[!annotation$chr %in% c('chrX', 'chrY'), ]
  annotation$chr.num = as.numeric(str_extract(annotation$chr, "[[:digit:]]+"))
  
  #### Extracting Pheno Data ####
  
  phenofile = data.frame(fread(paste(dir, "Data", "450k_Metadata.csv", 
                                     sep = '/'))) 
  phenofile$Gender = gsub(".*:","", phenofile$Gender)
  phenofile$Samples_IDs = gsub(".*_","", phenofile$Samples_IDs)
  phenofile$colIDs = paste(phenofile$Geo_IDs, phenofile$Samples_IDs, sep = '_')
  
  #### Organizing Beta Values Table ####
  
  beta_cols = data.frame(Geo_IDs = gsub("_.*","", names(betavals)))
  merged_pheno = merge(beta_cols, phenofile, by = 'Geo_IDs')
  names(betavals) = merged_pheno$colIDs
  betavals = betavals[rownames(betavals) %in% rownames(annotation), ]
  
  #### Saving Data ####
  
  Luad_data = list(betavals = betavals,
                   annotation = annotation)
  outfile_string = paste(outfile_name, ."RData", sep = '')
  save(Luad_data, file = paste(dir, "Data", 
                               outfile_string, sep = '/'))
}