### Figure 1 - Environmental Data ###
### This script will show you have to make plots that show changes in environmental variables (e.g., temperature, oxygen, etc.) as a function of depth. ###
### Last Updated: July 23, 2023 ###

# Load libraries
library(patchwork)
library(ggplot2)
library(reshape2)

# Load CTD data
ctd <- read.csv("CTD_Falkor.csv",header=TRUE)

# Chlorophyll
chl <- ggplot(ctd, aes(x=PRS, y=CHL, group=Eddy, color=Eddy)) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y = expression("Chlorophyll"~italic(a)~"("*mu*"g L"^{"-1"}*")"), x="Depth (m)")+scale_color_manual(breaks=c("Cyclonic","Anticyclonic"),labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))+geom_line()+coord_flip()+scale_y_continuous(position = "right")

# Temperature
temp <- ggplot(ctd, aes(x=PRS, y=TMP, group=Eddy, color=Eddy)) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y = "Temperature (°C)", x="Depth (m)")+scale_color_manual(breaks=c("Cyclonic","Anticyclonic"),labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))+geom_line()+coord_flip()+scale_y_continuous(position = "right",limits=c(0,25))

# Oxygen
oxy <- ggplot(ctd, aes(x=PRS, y=OXY, group=Eddy, color=Eddy)) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y = expression("Oxygen ("*mu*"mol L"^{"-1"}*")"), x="Depth (m)")+scale_color_manual(breaks=c("Cyclonic","Anticyclonic"),labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))+geom_line()+coord_flip()+scale_y_continuous(position = "right",limits=c(190,230))

# Density
den <- read.csv("Density.csv",header=TRUE)
density <- ggplot(den, aes(x=Pressure, y=Sigma, group=Eddy, color=Eddy)) +scale_x_reverse(position = "bottom")+theme_classic()+labs(y = expression("Potential Density (kg m"^{-3}*")"), x="Depth (m)")+scale_color_manual(breaks=c("Cyclonic","Anticyclonic"),labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))+geom_line()+coord_flip()+scale_y_continuous(position = "right")

# Nitrate
no3 <- read.csv("InorganicNutrients_Final.csv",header=TRUE)
no3$Site.. <- NULL
no3$Density <- NULL
no3m <- melt(no3,id.vars=c("EddyRep","Depth"))
n <- ggplot(no3m, aes(x=Depth, y=value,color=EddyRep))+scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+scale_color_manual(name="Eddy",breaks=c("Cyclonic","Anticyclonic"),labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))+geom_line(aes(linetype=variable))+coord_flip()+scale_y_continuous(position = "right")+geom_point()+labs(x="Depth (m)",y=expression("Nutrient ("*mu*"M)"))+scale_linetype_manual(values=c("dotted","solid"),labels=c(expression("PO"["4"]^{"3-"}),expression("NO"["3"]^{"-"}*"+ NO"["2"]^{"-"})),name="Nutrient")

# TOC 
toc <- read.csv("Falkor_Other_Data.csv",header=TRUE)
toc <- subset(toc,Parameter=="TOC")
carb <- ggplot(toc, aes(x=as.numeric(Depth), y=Measurement, group=as.factor(Site..), color=as.factor(Site..)))+scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+geom_line()+coord_flip()+scale_y_continuous(position = "right",limits=c(50,75))+geom_point()+labs(x="Depth (m)",y=expression("TOC ("*mu*"M)"))+scale_color_manual(name="Eddy",breaks=c(59,74),labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))
carb

# Combine all plots together, add a common legend, and add panel labels (b-g)
density+chl+temp+oxy+n+carb+plot_layout(guides = "collect")+plot_annotation(tag_levels=list(c("b","c","d","e","f","g")))
# ggsave("Figure1b_1g.pdf",width=13,height=8)
