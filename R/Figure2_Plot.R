
#### Figure 2 ####

library(ggplot2)
library(ggpubr)

d0.ind<-read.csv("plots/0_corr_indep.csv",header=T)
d0.dep<-read.csv("plots/0_corr_dep.csv",header=T)

d0.ind<-read.csv("plots/0.3_corr_indep.csv",header=T)
d0.dep<-read.csv("plots/0.3_corr_dep.csv",header=T)

d0.ind<-read.csv("plots/0.5_corr_indep.csv",header=T)
d0.dep<-read.csv("plots/0.5_corr_dep.csv",header=T)

d0.ind<-read.csv("plots/0.7_corr_indep.csv",header=T)
d0.dep<-read.csv("plots/0.7_corr_dep.csv",header=T)

df0<-cbind.data.frame("bf"=c(d0.ind$win1,d0.ind$win2,d0.ind$win3,d0.ind$win4,d0.ind$win5,d0.ind$win6,d0.ind$win7,d0.ind$win8,d0.ind$win9,d0.ind$win10, d0.dep$win1,d0.dep$win2,d0.dep$win3,d0.dep$win4,d0.dep$win5,d0.dep$win6,d0.dep$win7,d0.dep$win8,d0.dep$win9,d0.dep$win10),"windows"=c(rep(1:10,each=100), rep(1:10,each = 100)), "method" = c(rep("independent",1000), rep("dependent", 1000)))

g0<-ggplot(data = df0, aes(x=windows, y=bf, fill=method)) + 
  geom_boxplot(aes(fill=windows)) + labs(x = "Windows", y = "log(Bayes_Factor)") + geom_abline(slope=0,intercept=-5,col="blue") + theme_bw() + 
  theme(axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=0.5,vjust=0,face="plain"),
        legend.text=element_text(colour="black",size=20), legend.title =element_text(colour="black",size=20)) + ggtitle("Correlation: 0.7") +   theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none")

##### similarly g1 (for corr 0.3), g2 (for corr 0.5) and g3 (for corr 0.7) ########

g<-ggarrange(g0, g1, g2, g3, ncol = 2, nrow = 2, common.legend = TRUE) + theme_bw() + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position = "top")
annotate_figure(g,top = text_grob("Distribution of log(Bayes_Factor)", color = "black", face = "bold", size = 20))
