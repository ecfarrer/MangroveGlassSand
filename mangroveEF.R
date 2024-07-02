#Emily's mangrove code

library(tidyverse)
library(ggplot2)
library(nlme)
library(plotrix)
library(multcomp)
library(cowplot)

options(contrasts=c("contr.helmert","contr.poly"))
options("contrasts") #see what contrasts are set


##### Read in data #####

mangrove <- read.csv("~/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Data/MangroveSheetCorrectedEF.csv",stringsAsFactors = T)
mangrove$Substrate<-factor(mangrove$Substrate,levels=c("Glass","Mix","Dredge"))
mangrove$SubstrateInoculum<-factor(paste(mangrove$Substrate,mangrove$Inoculum,sep=""))

#View(mangrove)
head(mangrove)

#to use rhizovision explorer
#load image, push the blackwhite inverse button (very important) to make it white
#set the region of interest with the dotted line cutty button at the top
#in image preprocessing: Broken root, 600dpi, 200-215 threshold, filter non root up to 150mm (no fill holes or edge smoothing)
#root diameter categories: 0-.5, .5-1, 1-2, >2 
#preview segmented image
#then if it looks ok click the go arrow button at the top to process

  
#Removing measurements that were suspect
mangrove$RootVolumeML[which(mangrove$RootVolumeConfidence=="suspect")]<-NA
mangrove[which(mangrove$RootScanConfidence=="suspect"),35:62]<-NA

#Convert all units into grams rather than mg
mangrove$BelowGroundBiomassg<-mangrove$BelowGroundBiomassDry/1000
mangrove$AboveGroundBiomassg<-mangrove$AboveGroundBiomass/1000
mangrove$TotalBiomassg<-mangrove$TotalBiomass/1000

#Calculate root:shoot ratio
mangrove$RootShootRatio<-mangrove$BelowGroundBiomassg/mangrove$AboveGroundBiomassg

#See Kramer-Walter et al 2016 for calculations, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.12562#:~:text=SRL%20was%20calculated%20as%20root,divided%20by%20fresh%20root%20volume.
mangrove$FineRootLength<-mangrove$Root.Length.Diameter.Range.1.mm+mangrove$Root.Length.Diameter.Range.2.mm+mangrove$Root.Length.Diameter.Range.3.mm

#SRL (specific root length) was calculated as root length divided by root dry mass, m/g
mangrove$SRL<-mangrove$Total.Root.Length.mm/1000/mangrove$BelowGroundBiomassg

#Root tissue density (RTD), RTD was calculated as root dry mass divided by fresh root volume, g/ml
mangrove$RTDlab<-mangrove$BelowGroundBiomassg/mangrove$RootVolumeML
mangrove$RTDwinrhizo<-mangrove$BelowGroundBiomassg/mangrove$Volume.mm3

#Root branching intensity (RBI) was calculated as the number of root tips divided by root length, tips/cm
mangrove$RBI<-mangrove$Number.of.Root.Tips/(mangrove$Total.Root.Length.mm/10)




##### Growth results #####
#BelowGroundBiomassDry, AboveGroundBiomass, TotalBiomass, root:shoot ratio, Survival

###### Total Biomass ######  
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(TotalBiomassg),
    TotalBiomassg = mean(TotalBiomassg,na.rm=T))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/totalbiomass.pdf",width=2.2,height=2.2)
ptot<-ggplot(mangrove, aes(x=Substrate, y=TotalBiomassg,color=Inoculum)) +
  labs(x = "Substrate",y="Total biomass (g)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = TotalBiomassg-se, ymax = TotalBiomassg+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
dev.off()

ggplot(mangrove, aes(x=Fungus, y=TotalBiomassg)) +
  geom_boxplot()

m0<-gls(TotalBiomassg~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(TotalBiomassg~Substrate*Inoculum, random=~1|Fungus,data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
anova(m0,m1)


###### Aboveground biomass ######  
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(AboveGroundBiomassg),
    AboveGroundBiomassg = mean(AboveGroundBiomassg,na.rm=T))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/totalbiomass.pdf",width=2.2,height=2.2)
ggplot(mangrove, aes(x=Substrate, y=AboveGroundBiomassg,color=Inoculum)) +
  labs(x = "Substrate",y="Aboveground biomass (g)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = AboveGroundBiomassg-se, ymax = AboveGroundBiomassg+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
dev.off()

ggplot(mangrove, aes(x=Fungus, y=AboveGroundBiomassg)) +
  geom_boxplot()

m0<-gls(AboveGroundBiomassg~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(AboveGroundBiomassg~Substrate*Inoculum, random=~1|Fungus,data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
anova(m0,m1)

###### Belowground biomass ######  
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(BelowGroundBiomassg),
    BelowGroundBiomassg = mean(BelowGroundBiomassg,na.rm=T))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/totalbiomass.pdf",width=2.2,height=2.2)
ggplot(mangrove, aes(x=Substrate, y=BelowGroundBiomassg,color=Inoculum)) +
  labs(x = "Substrate",y="Aboveground biomass (g)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = BelowGroundBiomassg-se, ymax = BelowGroundBiomassg+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
dev.off()

ggplot(mangrove, aes(x=Fungus, y=BelowGroundBiomassg)) +
  geom_boxplot()

m0<-gls(BelowGroundBiomassg~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(BelowGroundBiomassg~Substrate*Inoculum, random=~1|Fungus,data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
anova(m0,m1)


###### Root:Shoot ######  
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(RootShootRatio),
    RootShootRatio = mean(RootShootRatio,na.rm=T))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/rootshootratio.pdf",width=2.2,height=2.2)
prootshoot<-ggplot(mangrove, aes(x=Substrate, y=RootShootRatio,color=Inoculum)) +
  labs(x = "Substrate",y="Root:shoot ratio") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = RootShootRatio-se, ymax = RootShootRatio+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
dev.off()
  
ggplot(mangrove, aes(x=Fungus, y=RootShootRatio)) +
  geom_boxplot()

m0<-gls(RootShootRatio~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(RootShootRatio~Substrate*Inoculum, random=~1|Fungus,data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
anova(m0,m1)


###### Survival ######
#standard deviation of a binary variable is sqrt(p*(1-p)), where p is the mean
# stderror is sqrt(p*(1-p)/n), see wikipedia, https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

dfsummary<-mangrove%>%
  group_by(Substrate,Inoculum)%>%
  summarize(seold=std.error(Survival, na.rm=T),Survival=mean(Survival,na.rm=T),n=sum(Survival))
dfsummary$se<-sqrt(dfsummary$Survival*(1-dfsummary$Survival)/10)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/survival.pdf",width=2.2,height=2.2)
psurv<-ggplot(mangrove, aes(x=Substrate, y=Survival,color=Inoculum)) +
  scale_y_continuous(name="Survival (proportion)",breaks=c(0,.25,.5,.75,1))+
  labs(x = "Substrate",y="Survival (proportion)") +
  #ylim(c(0,1.1))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = Survival-se, ymax = Survival+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
dev.off()

ms<-mangrove%>%
  group_by(Fungus)%>%
  summarize(mean=mean(Survival,na.rm=T),se=std.error(Survival, na.rm=T))

ggplot(ms, aes(x=Fungus, y=mean,col=Fungus)) +
  geom_point()+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se))


m1<-glm(Survival~Substrate*Inoculum,family=binomial(link="logit"),data=mangrove)
drop1(m1,test="Chisq",.~.)

hist(resid(m1))

library(lme4)
#There are convergence issues b/c the random effects variance is so close to 0, because there is no affect of fungus on survival. so I can't even do an anova on it.
m2<-glmer(Survival~Substrate*Inoculum+(1|Fungus),family=binomial,data=mangrove)
library(car)
Anova(m2,type="III")
anova(m2,m1)

mc<-glm(Survival~SubstrateInoculum,family=binomial(link="logit"),data=mangrove)
summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")))
mc2<-glm(Survival~Substrate,family=binomial(link="logit"),data=mangrove)
summary(glht(mc2, linfct = mcp(Substrate = "Tukey")))
mc3<-glm(Survival~Inoculum,family=binomial(link="logit"),data=mangrove)
summary(glht(mc3, linfct = mcp(Inoculum = "Tukey")))


###### Final plot of total biomass, root:shoot, survival ######
plegend<-ggplot(mangrove, aes(x=Substrate, y=Survival,color=Inoculum)) +
  labs(x = "Substrate",y="Survival (proportion)") +
  ylim(c(0,1.1))+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10))+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = Survival-se, ymax = Survival+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
legend <- cowplot::get_legend(plegend)

#for legend on the right
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/biomassratiosurvival.pdf",width=7.2,height=2.2)
plot_grid(ptot,prootshoot,psurv,legend, nrow = 1,labels=c("a) Total biomass","b) Root:shoot ratio","c) Survival",""),label_size=10,hjust=-.17,vjust=1.2,scale=1,label_fontface = "plain",rel_widths=c(1,1,1,.5))
dev.off()

#for legend on the bottom
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/biomassratiosurvival2.pdf",width=6.9,height=2.2)
plot_grid(ptot,prootshoot,psurv, nrow = 1,labels=c("a) Total biomass","b) Root:shoot ratio","c) Survival"),label_size=10,hjust=-.17,vjust=1.2,scale=1,label_fontface = "plain",rel_widths=c(1,1,1))
dev.off()

#note the panel a is wider b/c the y axis labels are smaller, i could try to fix this if I was feeling anal





##### Root scan results #####

#Total root length
ggplot(mangrove, aes(x=Substrate, y=Total.Root.Length.mm,fill=Inoculum)) +
  geom_boxplot()
#ggplot(subset(mangrove2,mangrove2$SRL<400), aes(x=Substrate, y=SRL,fill=Inoculum)) +geom_boxplot()
#ggplot(mangrove2, aes(x=Substrate, y=Total.Root.Length.mm,fill=Fungus)) +
#  geom_boxplot()

m1<-gls(Total.Root.Length.mm~Substrate*Inoculum+Fungus, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
m1<-lme(Total.Root.Length.mm~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
#summary(m1)

#need to add posthoc tests



#Branching frequency
ggplot(mangrove, aes(x=Substrate, y=Branching.frequency.per.mm,fill=Inoculum)) +
  geom_boxplot()

m1<-gls(Branching.frequency.per.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Avg diameter
ggplot(mangrove, aes(x=Substrate, y=Average.Diameter.mm,fill=Inoculum)) +
  geom_boxplot()

m1<-gls(Average.Diameter.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")


#Winrhizo's calculation of volume
ggplot(mangrove, aes(x=Substrate, y=Volume.mm3,fill=Inoculum)) +
  geom_boxplot()

m1<-gls(Volume.mm3~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Kathryn's measurement in lab
ggplot(mangrove, aes(x=Substrate, y=RootVolumeML,fill=Inoculum)) +
  geom_boxplot()

m1<-gls(RootVolumeML~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")


#Surface.Area.mm2
ggplot(mangrove, aes(x=Substrate, y=Surface.Area.mm2,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(Surface.Area.mm2~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Root length range 1 (0-.5mm)
ggplot(mangrove, aes(x=Substrate, y=Root.Length.Diameter.Range.1.mm,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(Root.Length.Diameter.Range.1.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Root length range 4 (>2)
ggplot(mangrove, aes(x=Substrate, y=Root.Length.Diameter.Range.4.mm,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(Root.Length.Diameter.Range.4.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Fine roots (<2)
ggplot(mangrove, aes(x=Substrate, y=FineRootLength,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(FineRootLength~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Specific root length
ggplot(mangrove, aes(x=Substrate, y=SRL,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(SRL~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Root tissue density
ggplot(mangrove, aes(x=Substrate, y=RTDwinrhizo,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(RTDwinrhizo~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Root branching intensity
ggplot(mangrove, aes(x=Substrate, y=RBI,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(RBI~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")












