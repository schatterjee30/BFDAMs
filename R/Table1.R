
#### Table 1 ####

########## load the results table of log-bases factor ############

tab<- read.csv("0_corr_indep.csv")
tab<- read.csv("0.3_corr_indep.csv")
tab<- read.csv("0.5_corr_indep.csv")
tab<- read.csv("0.7_corr_indep.csv")

#### for each tab run the following ######

tab<- tab[,-1]
truth<- c("Non-DMR",rep("DMR",6),"Non-DMR", rep("DMR",2))

pred<- matrix(NA,100,10)
for(i in 1:100){
  for(j in 1:10){
    if(tab[i,j] < (-5)){
      pred[i,j]<- "DMR"
    }else{
      pred[i,j]<- "Non-DMR"
    }
  }
}

for(i in 1:100){
  win1.miss<- sum(pred[,1]=="DMR")/100
  win2.miss<- sum(pred[,2]=="Non-DMR")/100
  win3.miss<- sum(pred[,3]=="Non-DMR")/100
  win4.miss<- sum(pred[,4]=="Non-DMR")/100
  win5.miss<- sum(pred[,5]=="Non-DMR")/100
  win6.miss<- sum(pred[,6]=="Non-DMR")/100
  win7.miss<- sum(pred[,7]=="Non-DMR")/100
  win8.miss<- sum(pred[,8]=="DMR")/100
  win9.miss<- sum(pred[,9]=="Non-DMR")/100
  win10.miss<- sum(pred[,10]=="Non-DMR")/100
}

