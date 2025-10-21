
#Panel A

wormsizer <- read.csv("wormsizer.csv", header = T)
wormsizer<-subset(wormsizer,pass=="TRUE")
wormsizer$Strain <- factor(wormsizer$Strain,levels = c("WT", "daf-2", "daf-16", "daf-16; daf-2"))
wormsizer$day <- as.character(wormsizer$day)
wormsizer <- subset(wormsizer, day == "1" | day == "6")

#Plot worm lenth
Wormsizer <- ggplot(wormsizer,aes(x=day,y=length))+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=1,binwidth=1,stackdir = "center", dotsize = 5)+
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 15, color = "black"),  # Adjust size of x-axis text
    axis.text.y = element_text(size = 12, color = "black"),  # Adjust size of y-axis text
    axis.title = element_text(size = 14, color = "black"),   # Adjust size of axis titles
    plot.title = element_text(size = 16, color = "black", face = "bold"),  # Adjust size and style of plot title
    strip.text = element_text(size = 12, color = "black"),   # Adjust size of facet strip text
    aspect.ratio = 1
  ) +
  labs(x="Strain",y="Length (um)") +
  facet_grid(.~Strain)
Wormsizer


#Statistics
library(nlme)
library(emmeans)

#Day 1 only
d1<-subset(wormsizer,day == "1")
d1<-subset(d1, Strain=="daf-16" | Strain=="daf-16; daf-2")
interaction<-lme(length~Strain,random=~1|replicate, data=d1) #the * notation is what we use for the interaction
summary(interaction)


#Interaction statistics
d1d6<-subset(wormsizer,Strain=="daf-16" | Strain=="daf-16; daf-2")
interaction<-lme(length~Strain*day,random=~1|replicate, data=d1d6) #the * notation is what we use for the interaction
summary(interaction)
#Change depending on comparison of interest



#Panel B

#Plot brood size
broodsize <- read.csv("Broodsize_EtOH.csv", header = T)
broodsize$Days_Starved<- as.character(broodsize$Days_Starved)
broodsize_plot<-ggplot(broodsize,aes(x=Days_Starved,y=Total))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=1,binwidth=2,stackdir = "center", dotsize = 1)+
  facet_grid(.~Strain) +
  theme_bw()+
  theme(legend.position="none")+
  stat_summary(
    fun = mean,
    geom = 'line',
    aes(group= Strain),
    position = position_dodge(width = 0.9),
    color = "black") +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=1, color="black")
broodsize_plot

#Day 1 only
d1<-subset(broodsize,Days_Starved == "1")
d1<-subset(d1, Strain=="daf-16; daf-2" | Strain=="daf-16")
interaction<-lme(Total~Strain,random=~1|Rep, data=d1) #the * notation is what we use for the interaction
summary(interaction)


#Interaction statistics
d1d6<-subset(broodsize,Strain=="WT" | Strain=="daf-2")
interaction<-lme(Total~Strain*Days_Starved,random=~1|Rep, data=d1d6) #the * notation is what we use for the interaction
summary(interaction)
#Change depending on comparison of interest
