#Emily's mangrove code

library(tidyverse)
library(ggplot2)
library(nlme)
library(plotrix)
library(multcomp)
library(cowplot)
library(piecewiseSEM)
library(performance)

options(contrasts=c("contr.helmert","contr.poly"))
options("contrasts") #see what contrasts are set


##### Read in data #####

mangrove <- read.csv("~/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Data/MangroveSheetCorrectedEF.csv",stringsAsFactors = T)
mangrove$Substrate<-factor(mangrove$Substrate,levels=c("Glass","Mix","Dredge"))
mangrove$SubstrateInoculum<-factor(paste(mangrove$Substrate,mangrove$Inoculum,sep=""))
mangrove$Fungus<-case_match(mangrove$Fungus,FALSE~"Uninfected",TRUE~"Infected")
mangrove$Fungus<-factor(mangrove$Fungus,levels=c("Uninfected","Infected"))

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
mangrove$FineRootLengthmm<-mangrove$Root.Length.Diameter.Range.1.mm+mangrove$Root.Length.Diameter.Range.2.mm+mangrove$Root.Length.Diameter.Range.3.mm

#SRL (specific root length) was calculated as root length divided by root dry mass, m/g
mangrove$SRL<-mangrove$Total.Root.Length.mm/1000/mangrove$BelowGroundBiomassg

#Root tissue density (RTD), RTD was calculated as root dry mass divided by fresh root volume, g/ml
mangrove$RTDlab<-mangrove$BelowGroundBiomassg/mangrove$RootVolumeML
mangrove$RTDwinrhizo<-mangrove$BelowGroundBiomassg/mangrove$Volume.mm3*1000

#Root branching intensity (RBI) was calculated as the number of root tips divided by root length, tips/cm
mangrove$RBI<-mangrove$Number.of.Root.Tips/(mangrove$Total.Root.Length.mm/10)


#Convert some root lengths to m rather than mm, and volume to ML
mangrove$Root.Length.Diameter.Range.4.m<-mangrove$Root.Length.Diameter.Range.4.m/1000
mangrove$FineRootLengthm<-mangrove$FineRootLengthmm/1000
mangrove$Volume.ml<-mangrove$Volume.mm3/1000




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
hist(resid(m1))
random.effects(m1)
anova(m0,m1)
anova(m1,type="marginal")
r2_nakagawa(m1)
rsquared(m1,method=NULL)

#Fungal infection figure
dfsummary <- mangrove %>%
  group_by(Fungus) %>%
  summarise(
    se = std.error(TotalBiomassg),
    TotalBiomassg = mean(TotalBiomassg,na.rm=T))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/totalbiomass.pdf",width=2.2,height=2.2)
ptotfungus<-
  ggplot(mangrove, aes(x=Fungus, y=TotalBiomassg)) +
  labs(x = "Fungal infection",y="Total biomass (g)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8)+
  geom_point(size=.7,alpha=rep(.4,60),position=position_jitter(width=0.2))+
  geom_errorbar(aes(ymin = TotalBiomassg-se, ymax = TotalBiomassg+se), data = dfsummary, width = 0.2)
#dev.off()



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
hist(resid(m1))
random.effects(m1)
anova(m0,m1)
anova(m1,type="marginal")
r2_nakagawa(m1)
rsquared(m1,method=NULL)


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
hist(resid(m1))
random.effects(m1)
anova(m0,m1)
anova(m1,type="marginal")
r2_nakagawa(m1)
rsquared(m1,method=NULL)



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
hist(resid(m1))
random.effects(m1)
anova(m0,m1)
anova(m1,type="marginal")
r2_nakagawa(m1)
rsquared(m1,method=NULL)



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

mangrove%>%
  group_by(Inoculum)%>%
  summarize(mean=mean(Survival,na.rm=T))
mangrove%>%
  group_by(Substrate)%>%
  summarize(mean=mean(Survival,na.rm=T))


m1<-glm(Survival~Substrate*Inoculum,family=binomial(link="logit"),data=mangrove)
drop1(m1,test="Chisq",.~.)
hist(resid(m1))
rsquared(m1,method=NULL)


library(lme4)
#There are convergence issues b/c the random effects variance is so close to 0, because there is no affect of fungus on survival. so I can't even do an anova on it.
m2<-glmer(Survival~Substrate*Inoculum+(1|Fungus),family=binomial,data=mangrove)
library(car)
Anova(m2,type="III")
anova(m2,m1)



mc<-glm(Survival~SubstrateInoculum,family=binomial(link="logit"),data=mangrove)
summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")))
summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")), adjusted('none'))
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
#To use:
# Fine root length - good
# Coarse root length - good
# Average.Diameter.mm - might be good not correlated with fine root length or total root length
# Volume.mm3 or RootVolumeML - good
# SRL - (specific root length) root length divided by root dry mass, m/g
# RTDlab or RTDwinrhizo - Root tissue density (RTD), root dry mass divided by fresh root volume, g/ml
# RBI - Root branching intensity (RBI), number of root tips divided by root length, tips/cm

#Not using
# Total root length - good, but very tightly correlated with fine root length
# Number.of.Root.Tips - correlated with total root length
# Number.of.Branch.Points - correlated with total root length
# Branching frequency - results are odd, weird interaction
# Network.Area.mm2 - correlated with total root length
# Median.Diameter.mm - kind of corrlated with avg diameter
# Maximum.Diameter.mm - meh
# Perimeter.mm - very tightly correlated with total root length 
# Surface.Area.mm2 - highly correlated with total root length
# projected area by size class - not sure what this is
# surface area by size class
# volume by size class


###### Fine root length (<2) ######

dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(FineRootLengthm),
    FineRootLengthm = mean(FineRootLengthm,na.rm=T))

pfine<-ggplot(mangrove, aes(x=Substrate, y=FineRootLengthm,color=Inoculum)) +
  labs(x = "Substrate",y="Fine root length (m)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = FineRootLengthm-se, ymax = FineRootLengthm+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) +
  annotate("text",x = c(1,2,3), y = 35, label = c("a","ab","b"),size=3)
             
mangrove %>%
  group_by(Substrate) %>%
  summarise(mean = mean(FineRootLengthm,na.rm=T))

ggplot(mangrove, aes(x=Substrate, y=FineRootLengthm,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Substrate, y=FineRootLengthm,fill=Fungus)) +
  geom_boxplot()

m0<-gls(FineRootLengthm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(FineRootLengthm~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)
# r2_nakagawa(m1)
# random.effects(m1)
# check_singularity(m1)

#with fungus as fixed, doesn't change anything
m0<-gls(FineRootLengthm~Substrate*Inoculum+Fungus, data=mangrove,na.action = na.omit)
anova(m0,type="marginal")
rsquared(m0)
summary(glht(m0, linfct = mcp(Substrate = "Tukey")))


mc<-lme(FineRootLengthm~Substrate*Inoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
summary(glht(mc, linfct = mcp(Substrate = "Tukey")))
mc<-lme(FineRootLengthm~SubstrateInoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")))

#Fungus figure
dfsummary <- mangrove %>%
  group_by(Fungus) %>%
  summarise(
    se = std.error(FineRootLengthm),
    FineRootLengthm = mean(FineRootLengthm,na.rm=T))

pfinefungus<-
  ggplot(mangrove, aes(x=Fungus, y=FineRootLengthm)) +
  labs(x = "Fungal infection",y="Fine root length (m)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8)+
  geom_point(size=.7,alpha=rep(.4,60),position=position_jitter(width=0.2))+
  geom_errorbar(aes(ymin = FineRootLengthm-se, ymax = FineRootLengthm+se), data = dfsummary, width = 0.2)


###### Coarse root length, range 4 (>2mm) ######
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(Root.Length.Diameter.Range.4.m),
    Root.Length.Diameter.Range.4.m = mean(Root.Length.Diameter.Range.4.m,na.rm=T))

pcoarse<-ggplot(mangrove, aes(x=Substrate, y=Root.Length.Diameter.Range.4.m,color=Inoculum)) +
  labs(x = "Substrate",y="Coarse root length (m)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = Root.Length.Diameter.Range.4.m-se, ymax = Root.Length.Diameter.Range.4.m+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
#  annotate("text",x = c(1,2,3), y = 35, label = c("a","ab","b"),size=3)

ggplot(mangrove, aes(x=Substrate, y=Root.Length.Diameter.Range.4.m,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Substrate, y=Root.Length.Diameter.Range.4.m,fill=Fungus)) +
  geom_boxplot()

m0<-gls(Root.Length.Diameter.Range.4.m~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(Root.Length.Diameter.Range.4.m~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)
r2_nakagawa(m1)

ggplot(mangrove, aes(x=Substrate, y=Root.Length.Diameter.Range.4.mm,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(Root.Length.Diameter.Range.4.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Fungus figure
dfsummary <- mangrove %>%
  group_by(Fungus) %>%
  summarise(
    se = std.error(Root.Length.Diameter.Range.4.m),
    Root.Length.Diameter.Range.4.m = mean(Root.Length.Diameter.Range.4.m,na.rm=T))

pcoarsefungus<-
  ggplot(mangrove, aes(x=Fungus, y=Root.Length.Diameter.Range.4.m)) +
  labs(x = "Fungal infection",y="Coarse root length (m)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8)+
  geom_point(size=.7,alpha=rep(.4,60),position=position_jitter(width=0.2))+
  geom_errorbar(aes(ymin = Root.Length.Diameter.Range.4.m-se, ymax = Root.Length.Diameter.Range.4.m+se), data = dfsummary, width = 0.2)



###### Average diameter ######

dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(Average.Diameter.mm),
    Average.Diameter.mm = mean(Average.Diameter.mm,na.rm=T))

pdia<-
  ggplot(mangrove, aes(x=Substrate, y=Average.Diameter.mm,color=Inoculum)) +
  labs(x = "Substrate",y="Average root dimater (mm)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = Average.Diameter.mm-se, ymax = Average.Diameter.mm+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) +
  annotate("text",x = c(1,2,3), y = .7, label = c("a","ab","b"),size=3)


ggplot(mangrove, aes(x=Substrate, y=Average.Diameter.mm,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Fungus, y=Average.Diameter.mm,fill=Fungus)) +
  geom_boxplot()

m0<-gls(Average.Diameter.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(Average.Diameter.mm~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)
# r2_nakagawa(m1)
# random.effects(m1)

mc<-lme(Average.Diameter.mm~Substrate*Inoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
summary(glht(mc, linfct = mcp(Substrate = "Tukey")))
# mc<-lme(FineRootLengthm~SubstrateInoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
# summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")))


###### Winrhizo's calculation of volume ######
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(Volume.ml),
    Volume.ml = mean(Volume.ml,na.rm=T))

pvolw<-
  ggplot(mangrove, aes(x=Substrate, y=Volume.ml,color=Inoculum)) +
  labs(x = "Substrate",y="Root volume (ml)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = Volume.ml-se, ymax = Volume.ml+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 

ggplot(mangrove, aes(x=Substrate, y=Volume.ml,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Fungus, y=Volume.ml,fill=Fungus)) +
  geom_boxplot()

m0<-gls(Volume.ml~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(Volume.ml~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)
# r2_nakagawa(m1)
# random.effects(m1)


###### Kathryn's volume measurement in lab ######
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(RootVolumeML),
    RootVolumeML = mean(RootVolumeML,na.rm=T))

pvolk<-
  ggplot(mangrove, aes(x=Substrate, y=RootVolumeML,color=Inoculum)) +
  labs(x = "Substrate",y="Root volume (ml)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = RootVolumeML-se, ymax = RootVolumeML+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 

ggplot(mangrove, aes(x=Substrate, y=RootVolumeML,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Fungus, y=RootVolumeML,fill=Fungus)) +
  geom_boxplot()

m0<-gls(RootVolumeML~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(RootVolumeML~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)

#Fungus figure
dfsummary <- mangrove %>%
  group_by(Fungus) %>%
  summarise(
    se = std.error(RootVolumeML),
    RootVolumeML = mean(RootVolumeML,na.rm=T))

pvolkfungus<-
  ggplot(mangrove, aes(x=Fungus, y=RootVolumeML)) +
  labs(x = "Fungal infection",y="Root volume (ml)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8)+
  geom_point(size=.7,alpha=rep(.4,60),position=position_jitter(width=0.2))+
  geom_errorbar(aes(ymin = RootVolumeML-se, ymax = RootVolumeML+se), data = dfsummary, width = 0.2)




###### Specific root length ######
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(SRL),
    SRL = mean(SRL,na.rm=T))

psrl<-
  ggplot(mangrove, aes(x=Substrate, y=SRL,color=Inoculum)) +
  labs(x = "Substrate",y="Specific root length (m/g)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = SRL-se, ymax = SRL+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) +
  annotate("text",x = c(1,2,3), y = 95, label = c("a","b","b"),size=3)

ggplot(mangrove, aes(x=Substrate, y=SRL,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Fungus, y=SRL,fill=Fungus)) +
  geom_boxplot()

m0<-gls(SRL~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(SRL~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)
# r2_nakagawa(m1)

mc<-lme(SRL~Substrate*Inoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
summary(glht(mc, linfct = mcp(Substrate = "Tukey")))
# mc<-lme(SRL~SubstrateInoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
# summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")))


#Root tissue density - winrhizo
ggplot(mangrove, aes(x=Substrate, y=RTDwinrhizo,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(RTDwinrhizo~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")


###### Root tissue density - lab/kathryn ######
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(RTDlab),
    RTDlab = mean(RTDlab,na.rm=T))

prtdk<-
  ggplot(mangrove, aes(x=Substrate, y=RTDlab,color=Inoculum)) +
  labs(x = "Substrate",y="Root tissue density (g/ml)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = RTDlab-se, ymax = RTDlab+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
#  annotate("text",x = c(1,2,3), y = 95, label = c("a","b","b"),size=3)

ggplot(mangrove, aes(x=Substrate, y=RTDlab,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Fungus, y=RTDlab,fill=Fungus)) +
  geom_boxplot()

m0<-gls(RTDlab~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(RTDlab~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)
# r2_nakagawa(m1)

mc<-lme(RTDlab~Substrate*Inoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
summary(glht(mc, linfct = mcp(Substrate = "Tukey")))
# mc<-lme(SRL~SubstrateInoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
# summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")))


###### Root branching intensity ######
dfsummary <- mangrove %>%
  group_by(Substrate, Inoculum) %>%
  summarise(
    se = std.error(RBI),
    RBI = mean(RBI,na.rm=T))

mangrove %>%
  group_by(Substrate, Inoculum,Fungus) %>%
  summarise(n = sum(SRL>0,na.rm=T)) #low survival and the scan was bad for mix autoclaved (n=4)
    
prbi<-
  ggplot(mangrove, aes(x=Substrate, y=RBI,color=Inoculum)) +
  labs(x = "Substrate",y="Root branching intensity (tips/cm)") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),legend.position = "none")+
  geom_point(data=dfsummary,size=1.8,aes(color=Inoculum,group=Inoculum),position = position_dodge(0.8))+
  geom_point(size=.7,aes(fill=Inoculum,group=Inoculum),position = position_jitterdodge(jitter.width=.2),alpha=rep(.4,60))+
  geom_errorbar(aes(ymin = RBI-se, ymax = RBI+se), data = dfsummary, width = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values = c("#6c66be", "#8ca54f")) 
#  annotate("text",x = c(1,2,3), y = 95, label = c("a","b","b"),size=3)

ggplot(mangrove, aes(x=Substrate, y=RBI,fill=Inoculum)) +
  geom_boxplot()
ggplot(mangrove, aes(x=Fungus, y=RBI,fill=Fungus)) +
  geom_boxplot()

m0<-gls(RBI~Substrate*Inoculum, data=mangrove,na.action = na.omit)
m1<-lme(RBI~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m0,m1)
anova(m1,type="marginal")
rsquared(m1)
# r2_nakagawa(m1)

#mc<-lme(RBI~Substrate*Inoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
#summary(glht(mc, linfct = mcp(Substrate = "Tukey")))
mc<-lme(RBI~SubstrateInoculum,random=~1|Fungus,data=mangrove,na.action=na.omit)
summary(glht(mc, linfct = mcp(SubstrateInoculum = "Tukey")))
summary(glht(mc, linfct = mcp(SubstrateInoculum=c("GlassAutoclaved-GlassLive=0","MixAutoclaved-MixLive=0","DredgeAutoclaved-DredgeLive=0"))))
summary(glht(mc, linfct = mcp(SubstrateInoculum=c("GlassAutoclaved-DredgeAutoclaved=0","GlassAutoclaved-MixAutoclaved=0","DredgeAutoclaved-MixAutoclaved=0","GlassLive-DredgeLive=0","GlassLive-MixLive=0","DredgeLive-MixLive=0"))))




###### Final plot of root properties ######
dfsummary<-mangrove%>%
  group_by(Substrate,Inoculum)%>%
  summarize(seold=std.error(Survival, na.rm=T),Survival=mean(Survival,na.rm=T),n=sum(Survival))
dfsummary$se<-sqrt(dfsummary$Survival*(1-dfsummary$Survival)/10)

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
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/rootproperties.pdf",width=7.2,height=7.2)
plot_grid(pfine,pcoarse,pdia,pvolk,psrl,prtdk,prbi,legend, nrow = 3,labels=c("a) Fine root length","b) Coarse root length","c) Average root diameter","d) Root volume","e) Specific root length","f) Root tissue density","g) Root branching intensity",""),label_size=10,hjust=-.17,vjust=1.2,scale=1,label_fontface = "plain")
dev.off()

#for legend on the bottom
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/biomassratiosurvival2.pdf",width=6.9,height=2.2)
plot_grid(ptot,prootshoot,psurv, nrow = 1,labels=c("a) Total biomass","b) Root:shoot ratio","c) Survival"),label_size=10,hjust=-.17,vjust=1.2,scale=1,label_fontface = "plain",rel_widths=c(1,1,1))
dev.off()

#note the panel a is wider b/c the y axis labels are smaller, i could try to fix this if I was feeling anal



# Fungus figure
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Manuscript/Figs/fungus.pdf",width=7.2,height=1.8)
plot_grid(ptotfungus,pcoarsefungus,pfinefungus,pvolkfungus, nrow = 1,labels=c("a) Total biomass","b) Coarse root length","c) Fine root length","d) Root volume"),label_size=10,hjust=-.17,vjust=1.2,scale=1,label_fontface = "plain")
dev.off()




###### Other root properties ######

#Total root length
ggplot(mangrove, aes(x=Substrate, y=Total.Root.Length.mm,fill=Inoculum)) +
  geom_boxplot()
#ggplot(subset(mangrove2,mangrove2$SRL<400), aes(x=Substrate, y=SRL,fill=Inoculum)) +geom_boxplot()
#ggplot(mangrove2, aes(x=Substrate, y=Total.Root.Length.mm,fill=Fungus)) +
#  geom_boxplot()
plot(mangrove$Total.Root.Length.mm,mangrove$FineRootLength)

m1<-gls(Total.Root.Length.mm~Substrate*Inoculum+Fungus, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
m1<-lme(Total.Root.Length.mm~Substrate*Inoculum, random=~1|Fungus, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")
#summary(m1)


# Number of root tips
ggplot(mangrove, aes(x=Substrate, y=Number.of.Root.Tips,fill=Inoculum)) +
  geom_boxplot()
plot(mangrove$Number.of.Root.Tips,mangrove$Total.Root.Length.mm)

# Number of branch points
ggplot(mangrove, aes(x=Substrate, y=Number.of.Branch.Points,fill=Inoculum)) +
  geom_boxplot()
plot(mangrove$Number.of.Branch.Points,mangrove$Total.Root.Length.mm)

# Network area
ggplot(mangrove, aes(x=Substrate, y=Network.Area.mm2,fill=Inoculum)) +
  geom_boxplot()
plot(mangrove$Network.Area.mm2,mangrove$Total.Root.Length.mm)

# Perimeter.mm 
ggplot(mangrove, aes(x=Substrate, y=Perimeter.mm,fill=Inoculum)) +
  geom_boxplot()
plot(mangrove$Perimeter.mm,mangrove$Total.Root.Length.mm)

#Branching frequency
ggplot(mangrove, aes(x=Substrate, y=Branching.frequency.per.mm,fill=Inoculum)) +
  geom_boxplot()

m1<-gls(Branching.frequency.per.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

#Surface.Area.mm2
ggplot(mangrove, aes(x=Substrate, y=Surface.Area.mm2,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(Surface.Area.mm2~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

plot(mangrove$Surface.Area.mm2,mangrove$Total.Root.Length.mm)

#Root length range 1 (0-.5mm)
ggplot(mangrove, aes(x=Substrate, y=Root.Length.Diameter.Range.1.mm,fill=Inoculum)) +
  geom_boxplot()
m1<-gls(Root.Length.Diameter.Range.1.mm~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")













