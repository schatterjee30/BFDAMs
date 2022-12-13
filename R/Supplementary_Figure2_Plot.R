
#### Supplementary Figure 2 ####

##### For each chromosome fit NSC and plot the fitted values #####

Supp.figure2 = function(bvals){
  burn=1000
  iter= 20000
  treat.lab = 'Tumor'
  control.lab = 'Normal'
  nwin<- length(unique(bvals$windows))
  
  #### Applying Independent method for all windows ####
  
  set.seed(1234)
  for(i in 1:22){
    message('chromosome', i, ' Running')
    for(w in 1:length(unique(bvals$windows))){
      data<- bvals
      datsub=data[data$windows==w,]
      cols=ncol(datsub)
      n=nrow(datsub)
      x1=c(1:n)
      yall=datsub[,2:cols]
      
      #### Log-Transforming Methylation data ####
      
      yall<- as.matrix(logit2(yall))
      
      #### Separating Normal and Cancer groups data ####
      
      ynorm=yall[,grep(control.lab,names(yall))]
      ycanc=yall[,grep(treat.lab,names(yall))]
      
      #### Combined group mean fitting ####
      
      gm_all<-rep(0,n)
      gm_all<- apply(yall,1,mean)
      fitg <- ncs(x1,gm_all,burn,iter)
      overall_mean_fit<-apply(fitg$fhat, 2, mean)
      
      #### Normal group mean fitting ####
      
      gm_norm<-rep(0,n)
      gm_norm<- apply(ynorm,1,mean)
      fitgnorm <- ncs(x1,gm_norm,burn,iter)
      norm_mean_fit<-apply(fitgnorm$fhat, 2, mean)
      
      #### Cancer group mean fitting ####
      
      gm_canc<-rep(0,n)
      gm_canc<- apply(ycanc,1,mean)
      fitgcanc <- ncs(x1,gm_canc,burn,iter)
      canc_mean_fit<-apply(fitgcanc$fhat, 2, mean)
    }
    
    #### Plotting fit results by chromosome ####
    
    df.overall<- cbind.data.frame("windows" = data$windows, "fitted" = overall_mean_fit)
    df.norm<- cbind.data.frame("windows" = data$windows, "fitted" = norm_mean_fit)
    df.canc<- cbind.data.frame("windows" = data$windows, "fitted" = canc_mean_fit)
    
    plot(df.overall$windows, df.overall$fitted, type="b")
    lines(df.norm$windows, df.norm$fitted)
    lines(df.canc$windows, df.canc$fitted)
  }
}