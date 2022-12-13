
#### Supplementary Figure 1 ####

bump<-read.csv("plots/bumps_0corr_simul.csv",header=T)  
bump<-read.csv("plots/bumps_0.3corr_simul.csv",header=T)  
bump<-read.csv("plots/bumps_0.5corr_simul.csv",header=T)  
bump<-read.csv("plots/bumps_0.7corr_simul.csv",header=T) 

df0<-cbind.data.frame("bf"=c(bump$win1,bump$win2,bump$win3,bump$win4,bump$win5,bump$win6,bump$win7,bump$win8,bump$win9,bump$win10),"windows"=rep(1:10,each=100))
df0$windows<-as.factor(df0$windows)

#### Plotting ####

g0<-ggplot(data = df0, aes(x=windows, y=bf, fill=windows)) + 
  geom_boxplot(aes(fill=windows)) + labs(x = "Windows", y = "Number of bumps") + geom_abline(slope=0,intercept=6,col="blue") + theme_bw() + 
  theme(axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=0.5,vjust=0,face="plain"),
        legend.text=element_text(colour="black",size=20), legend.title =element_text(colour="black",size=20)) + ggtitle("Correlation: 0.7") +   theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none")

##### similarly g1 (for corr 0.3), g2 (for corr 0.5) and g3 (for corr 0.7) ########

g<-ggarrange(g0, g1, g2, g3, ncol = 2, nrow = 2, common.legend = TRUE) + theme_bw() + theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(legend.position = "none")
annotate_figure(g,top = text_grob("Distribution of number of bumps", color = "black", face = "bold", size = 20))
