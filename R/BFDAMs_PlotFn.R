
#### Description ####

# beta.df = Data frame of methylation beta values
# chr.pos = Chromosome and position of CpG sites
# nCpGs = Number of CpGs desired to be in one region
# iter = Number of MCMC iterations to be run
# burn = Number of burn-ins
# seed = A random number for reproducibilty of results
# control.group = Label for samples in group 1 in methylation data frame(beta.df)
# case.group = Label for samples in group 2 in methylation data frame(beta.df)
# chr = Chromosome number where the desired window is located
# window.plot = Desired window/genomic region to plot

#### Required Libraries ####

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(minfi)))

#### Plotting Fitted Curves ####

plot.fittedCurves<-function(beta.df, 
                            chr.pos, 
                            nCpGs = 100, 
                            iter = 500, 
                            burn = 50, 
                            seed = 1234, 
                            control.group = 'Normal', 
                            case.group = 'Tumor',
                            chr=12, 
                            window.plot=1){
  
  #### Organizing Input Data ####
  
  set.seed(seed)
  beta.df<-beta.df[!duplicated(rownames(beta.df)),]
  chr.pos<-chr.pos[!duplicated(rownames(chr.pos)),]
  beta.df<-beta.df[rownames(beta.df) %in% rownames(chr.pos),]
  if(is.character(chr.pos$chr)){
    chr.pos$chr<-readr::parse_number(chr.pos$chr)
  }
  if(is.character(chr.pos$pos)){
    chr.pos$pos<-as.numeric(chr.pos$pos)
  }
  
  #### Assigning Equally Spaced and Sized Windows ####
  
  split.chrpos<-split(chr.pos, chr.pos$chr)
  chrpos.win<-lapply(split.chrpos, function(x){
    x$windows<-rep(seq_len(nrow(x)%/%nCpGs+1L),each=nCpGs,len=nrow(x))
    return(x)
  })
  
  #### Fitting Desired Window by Desired Chromosome ####
  
  chr<-as.character(chr)
  id<-intersect(names(chrpos.win),chr)
  if(identical(id, character(0))){
    stop("Chromosome number mentioned does not exist in data")
  }else{
    run.chrpos<-chrpos.win[[id]]
  }
  win.chrpos<-run.chrpos[run.chrpos$windows == window.plot,]
  win.beta<-beta.df[rownames(beta.df) %in% rownames(win.chrpos),]
  
  #### Log-Transforming Methylation data ####
  
  n<-nrow(win.beta)
  x1<-c(1:n)
  yall<-data.frame(logit2(win.beta))
  
  #### Separating Normal and Cancer groups data ####
  
  ynorm<-yall[,grep(control.group,names(yall))]
  ycanc<-yall[,grep(case.group,names(yall))]
  
  #### Combined group mean fitting ####
  
  gm_all<-rep(0,n)
  gm_all<- apply(yall,1,meanfn)
  fitg <- ncs(x1,gm_all,burn,iter)
  overall_mean_fit<-fitg$fhat
  overall_sigma_fit<-fitg$sig2
  overall_tau_fit<-fitg$tau
  lambda_overall<-cbind(overall_tau_fit, overall_sigma_fit)
  
  #### Normal group mean fitting ####
  
  gm_norm<-rep(0,n)
  gm_norm<- apply(ynorm,1,meanfn)
  fitgnorm <- ncs(x1,gm_norm,burn,iter)
  norm_mean_fit<-fitgnorm$fhat
  norm_sigma_fit<-fitgnorm$sig2
  norm_tau_fit<-fitg$tau
  lambda_norm<-cbind(norm_tau_fit, norm_sigma_fit)
  
  #### Cancer group mean fitting ####
  
  gm_canc<-rep(0,n)
  gm_canc<- apply(ycanc,1,meanfn)
  fitgcanc <- ncs(x1,gm_canc,burn,iter)
  canc_mean_fit<-fitgcanc$fhat
  canc_sigma_fit<-fitgcanc$sig2
  canc_tau_fit<-fitg$tau
  lambda_canc<-cbind(canc_tau_fit, canc_sigma_fit)
  
  #### Plotting Curves ####
  
  normfit<-apply(norm_mean_fit,2,mean)
  cancfit<-apply(canc_mean_fit,2,mean)
  CpG.Site<-c(1:length(normfit))
  pdata<-data.frame(CpG.Site=CpG.Site,Mean.Methylation=gm_all, fit.normal=normfit,fit.cancer=cancfit)
  colors<-c("fit.normal" = "blue", "fit.cancer" = "red")
  p<-suppressMessages(ggplot(pdata, aes(x=CpG.Site)) + 
    geom_point(aes(y = Mean.Methylation), color = "green") +
    geom_smooth(aes(y = fit.normal), color = "blue", se=FALSE) + 
    geom_smooth(aes(y = fit.cancer), color="red", se=FALSE) + 
    theme_bw() + 
    labs(color = "Legend") +
    scale_color_manual(values = colors) +
    scale_x_continuous("CpG.Site") + 
    scale_y_continuous("Mean.Methylation(M-values)"))
  return(p)
}

