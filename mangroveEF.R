#Emily's mangrove code

library(tidyverse)
library(ggplot2)
library(nlme)
library(plotrix)

options(contrasts=c("contr.helmert","contr.poly"))
options("contrasts") #see what contrasts are set


##### Read in data #####

mangrove <- read.csv("~/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Mangroves/Data/MangroveSheetCorrectedEF.csv")
mangrove$Substrate<-factor(mangrove$Substrate,levels=c("Glass","Mix","Dredge"))

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


#See Kramer-Walter et al 2016 for calculations, https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.12562#:~:text=SRL%20was%20calculated%20as%20root,divided%20by%20fresh%20root%20volume.
mangrove$FineRootLength<-mangrove$Root.Length.Diameter.Range.1.mm+mangrove$Root.Length.Diameter.Range.2.mm+mangrove$Root.Length.Diameter.Range.3.mm

#SRL (specific root length) was calculated as root length divided by root dry mass
mangrove$SRL<-mangrove$Total.Root.Length.mm/mangrove$BelowGroundBiomassDry

#Root tissue density (RTD), RTD was calculated as root dry mass divided by fresh root volume
mangrove$RTDlab<-mangrove$BelowGroundBiomassDry/mangrove$RootVolumeML
mangrove$RTDwinrhizo<-mangrove$BelowGroundBiomassDry/mangrove$Volume.mm3

#Root branching intensity (RBI) was calculated as the number of root tips divided by root length
mangrove$RBI<-mangrove$Number.of.Root.Tips/mangrove$Total.Root.Length.mm




##### Growth results #####
#BelowGroundBiomassDry, AboveGroundBiomass, TotalBiomass, SurvivalEF

ggplot(mangrove, aes(x=Substrate, y=TotalBiomass,fill=Inoculum)) +
  geom_boxplot()


m1<-gls(TotalBiomass~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

m1<-gls(BelowGroundBiomassDry~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")

m1<-gls(AboveGroundBiomass~Substrate*Inoculum, data=mangrove,na.action = na.omit)
anova(m1,type="marginal")


#Survival
ms<-mangrove%>%
  group_by(Substrate,Inoculum)%>%
  summarize(mean=mean(Survival,na.rm=T),se=std.error(Survival, na.rm=T))

ggplot(ms, aes(x=Substrate, y=mean,col=Inoculum)) +
  geom_point()+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se))

m1<-glm(Survival~Substrate*Inoculum,family=binomial(link="logit"),data=mangrove)
drop1(m1,test="Chisq",.~.)



##### Root results #####

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












