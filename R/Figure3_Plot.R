
#### Figure 3 ####

load("realData_windowInfo.Rdata")

#### Loading Bumphunter Results ####

load("Bumphunter_Luad.RData")
c$pred_bump<- rep(NA, length(c$windows))
for(i in 1:length(c$windows)){
  print(i)
  l=0
  for(j in 1:length(fit_df$chr)){
    if(fit_df$start[j] > c$start[i] & fit_df$end[j] < c$end[i]){
        l=l+1
    }
  }
  if(l > 20){
    c$pred_bump[i]<-  "DMR"
  }else{
    c$pred_bump[i]<- "Non-DMR"
  }
}

#### Loading Comb-p Results ####

load("Comb-p_Luad.RData")
c$pred_combp<- rep(NA, length(c$windows))
for(i in 1:length(c$windows)){
  print(i)
  l=0
  for(j in 1:length(fit_df$chr)){
    if(fit_df$start[j] > c$start[i] & fit_df$end[j] < c$end[i]){
      l=l+1
    }
  }
  if(l > 20){
    c$pred_combp[i]<-  "DMR"
  }else{
    c$pred_combp[i]<- "Non-DMR"
  }
}

#### Loading DMRcate Results ####

load("DMRcate_Luad.RData")
c$pred_dmrcate<- rep(NA, length(c$windows))
for(i in 1:length(c$windows)){
  print(i)
  l=0
  for(j in 1:length(fit_df$chr)){
    if(fit_df$start[j] > c$start[i] & fit_df$end[j] < c$end[i]){
      l=l+1
    }
  }
  if(l > 20){
    c$pred_dmrcate[i]<-  "DMR"
  }else{
    c$pred_dmrcate[i]<- "Non-DMR"
  }
}

######### BFDAM2 result ########

b1<-read.csv("bayes1.csv")
b2<-read.csv("bayes2.csv")
b3<-read.csv("bayes3.csv")
b4<-read.csv("bayes4.csv")
b5<-read.csv("bayes5.csv")
b6<-read.csv("bayes6.csv")
b7<-read.csv("bayes7.csv")
b8<-read.csv("bayes8.csv")
b9<-read.csv("bayes9.csv")
b10<-read.csv("bayes10.csv")
b11<-read.csv("bayes11.csv")
b12<-read.csv("bayes12.csv")
b13<-read.csv("bayes13.csv")
b14<-read.csv("bayes14.csv")
b15<-read.csv("bayes15.csv")
b16<-read.csv("bayes16.csv")
b17<-read.csv("bayes17.csv")
b18<-read.csv("bayes18.csv")
b19<-read.csv("bayes19.csv")
b20<-read.csv("bayes20.csv")
b21<-read.csv("bayes21.csv")
b22<-read.csv("bayes22.csv")
b23<-read.csv("bayes23.csv")

bayes<- rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23)
bayes$window_new<- paste0("Window",1:dim(bayes)[1])
cutoff<- -5
bayes$pred_bfdm<- ifelse(bayes$bayesfactor < cutoff, "DMR", "Non-DMR")
c$window_new<- paste0("Window", 1:length(c$window))

########## upsetR plot #############

library(UpSetR)
genelist<- list("BFDAM2"=bayes$window_new[bayes$pred_bfdm=="DMR"],"Bumphunter"=c$window_new[c$pred_bump=="DMR"], "DMRcate"=c$window_new[c$pred_dmrcate=="DMR"], "Comb-p"=c$window_new[c$pred_combp=="DMR"])
upset(fromList(genelist), order.by = "freq", main.bar.color = "black", sets.bar.color = "black", text.scale = 2)
