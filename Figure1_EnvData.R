### Figure 1 - Environmental Data ###
### This script will show you have to make plots that show changes in environmental variables (e.g., temperature, oxygen, etc.) as a function of depth. ###
### Last Updated: November 5, 2023 ###

# Load Data and Add Empty Points
ctd <- read.csv("Falkor_CTD_Oct2023.csv",header=TRUE)
ctd$Point <- ifelse(ctd$Depth==25 & ctd$Site==85,1,0)
ctd$Point <- ifelse(ctd$Depth==25 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==150 & ctd$Site==85,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==150 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==250 & ctd$Site==85,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==250 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==112 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==122 & ctd$Site==85,1,ctd$Point)
ctd$Point <- ifelse(ctd$Point!=0,ctd$Density,NA)

# Density
Den <- ggplot(ctd, aes(x=Depth, y=Density, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y=expression("Density (kg m"^{-3}*")"), x="Depth (m)")+geom_line()+coord_flip(xlim=c(250,0), ylim=c(23.5,26), clip="off")+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+scale_color_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=ctd,aes(x=Depth,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_y_continuous(position = "right")
Den

# Chlorophyll
ctd$CHL2 <- (ctd$CHL-0.10896)/1.7389
ctd$Point <- ifelse(!is.na(ctd$Point),ctd$CHL2,NA)
Chl <- ggplot(ctd, aes(x=Depth, y=CHL2, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y = expression("Chlorophyll Fluorescence (RFU)"), x="Depth (m)")+geom_line()+coord_flip(xlim=c(250,0), ylim=c(-0.05,1), clip="off")+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+scale_color_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=ctd,aes(x=Depth,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_y_continuous(position = "right")
Chl

# PAR
# new <- subset(ctd,Site!=73)
# new$Point <- ifelse(!is.na(new$Point),new$PAR,NA)
# Par <- ggplot(new,aes(x=Depth, y=PAR, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y = expression("PAR ("*mu*"mol photons m"^{"-2"}*"s"^{"-1"}*")"), x="Depth (m)")+geom_line()+coord_flip(xlim=c(250,0), ylim=c(0,628.91), clip="off")+scale_y_continuous(position = "right")+scale_color_manual(name="Site",labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))+geom_point(data=new,aes(x=Depth,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_y_continuous(position = "right")+scale_fill_manual(name="Site",labels=c("Cyclonic","Anticyclonic"), values=c("blue","red"))
# Par

# Temp
ctd$Point <- ifelse(!is.na(ctd$Point),ctd$TMP,NA)
Tmp <- ggplot(ctd, aes(x=Depth, y=TMP, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y = expression("Temperature (°C)"), x="Depth (m)")+geom_line()+coord_flip(xlim=c(250,0), ylim=c(0,25), clip="off")+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+scale_color_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=ctd,aes(x=Depth,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_y_continuous(position = "right")
Tmp

# Oxygen
ctd$Point <- ifelse(!is.na(ctd$Point),ctd$OXY,NA)
Oxy <- ggplot(ctd, aes(x=Depth, y=OXY, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+labs(y = expression("Oxygen ("*mu*"mol L"^{-1}*")"), x="Depth (m)")+geom_line()+coord_flip(xlim=c(250,0), ylim=c(190,230), clip="off")+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+scale_color_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=ctd,aes(x=Depth,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_y_continuous(position = "right")
Oxy

# Nutrients
no3 <- read.csv("Falkor_Nutrients_Oct2023.csv",header=TRUE)
no3m <- melt(no3,id.vars=c("Site","Depth","Density"))
tmp <- c(74,122,NA,"P",NA)
tmp2 <- c(59,112,NA,"P",NA)
tmp3 <- c(74,122,NA,"N",NA)
tmp4 <- c(59,112,NA,"N",NA)
no3m <- rbind(no3m,tmp,tmp2,tmp3,tmp4)
no3m$Depth <- as.numeric(no3m$Depth)
no3m$value <- as.numeric(no3m$value)
no3m$Site <- as.factor(no3m$Site)
no3m$Density <- as.numeric(no3m$Density)

no3mNC <- subset(no3m,variable=="N" & Site=="59")
Intersect <- approxfun(no3mNC$Depth,no3mNC$value)
no3mNC <- Intersect(112)

no3mPC <- subset(no3m,variable=="P" & Site=="59")
Intersect <- approxfun(no3mPC$Depth,no3mPC$value)
no3mPC <- Intersect(112)

no3mNA <- subset(no3m,variable=="N" & Site=="74")
Intersect <- approxfun(no3mNA$Depth,no3mNA$value)
no3mNA <- Intersect(122)

no3mPA <- subset(no3m,variable=="P" & Site=="74")
Intersect <- approxfun(no3mPA$Depth,no3mPA$value)
no3mPA <- Intersect(122)

no3m$Point <- ifelse(no3m$Depth==25 & no3m$Site==74,1,0)
no3m$Point <- ifelse(no3m$Depth==25 & no3m$Site==59,1,no3m$Point)
no3m$Point <- ifelse(no3m$Depth==150 & no3m$Site==74,1,no3m$Point)
no3m$Point <- ifelse(no3m$Depth==150 & no3m$Site==59,1,no3m$Point)
no3m$Point <- ifelse(no3m$Depth==250 & no3m$Site==74,1,no3m$Point)
no3m$Point <- ifelse(no3m$Depth==250 & no3m$Site==59,1,no3m$Point)
no3m$Point <- ifelse(no3m$Depth==112 & no3m$Site==59,1,no3m$Point)
no3m$Point <- ifelse(no3m$Depth==122 & no3m$Site==74,1,no3m$Point)
no3m$value <- ifelse(no3m$Site==59 & no3m$Depth==112 & no3m$variable=="N",no3mNC,no3m$value)
no3m$value <- ifelse(no3m$Site==59 & no3m$Depth==112 & no3m$variable=="P",no3mPC,no3m$value)
no3m$value <- ifelse(no3m$Site==74 & no3m$Depth==122 & no3m$variable=="N",no3mNA,no3m$value)
no3m$value <- ifelse(no3m$Site==74 & no3m$Depth==122 & no3m$variable=="P",no3mPA,no3m$value)

no3m$Point <- ifelse(no3m$Point!=0,no3m$value,NA)
no3m$Point <- as.numeric(no3m$Point)

n <- ggplot(data =  subset(no3m, !is.na(value)), aes(x=Depth, y=value,color=as.factor(Site),linetype=variable),group=variable)+scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+geom_line()+coord_flip(xlim=c(250,0), ylim=c(0,3), clip="off")+scale_y_continuous(position = "right")+geom_point(shape=15)+labs(x="Depth (m)",y=expression("Nutrient ("*mu*"M)"))+scale_linetype_manual(values=c("dotted","solid"),labels=c(expression("PO"["4"]^{"3-"}),expression("NO"["3"]^{"-"}*"+ NO"["2"]^{"-"})),name="Nutrient")+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+scale_color_manual(name="Site",labels=c("Cyclonic","SLA ≈0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=no3m,aes(x=Depth,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)
n

# TOC
toc <- read.csv("Falkor_TOC_Oct2023.csv",header=TRUE)
toc <- left_join(toc,no3)
tmp <- c(74,122,NA,"P","N",NA)
tmp2 <- c(59,112,NA,"P","N",NA)
toc <- rbind(toc,tmp,tmp2)
toc$Depth <- as.numeric(toc$Depth)
toc$TOC <- as.numeric(toc$TOC)
toc$Site <- as.factor(toc$Site)
toc$Density <- as.numeric(toc$Density)

tocNC <- subset(toc,Site=="59")
Intersect <- approxfun(tocNC$Depth,tocNC$TOC)
tocNC <- Intersect(112)

tocNA <- subset(toc,Site=="74")
Intersect <- approxfun(tocNA$Depth,tocNA$TOC)
tocNA <- Intersect(122)

toc$Point <- ifelse(toc$Depth==25,1,NA)
toc$Point <- ifelse(toc$Depth==150,1,toc$Point)
toc$Point <- ifelse(toc$Depth==250,1,toc$Point)
toc$Point <- ifelse(toc$Site==59 & toc$Depth==112,tocNC,toc$Point)
toc$Point <- ifelse(toc$Site==74 & toc$Depth==122,tocNA,toc$Point)

toc$Point <- ifelse(toc$Point==1,toc$TOC,toc$Point)
toc$Point <- ifelse(toc$Site==73,NA,toc$Point)
toc$Point <- as.numeric(toc$Point)

tocP <- ggplot(data =  subset(toc, !is.na(TOC)), aes(x=Depth, y=TOC,color=as.factor(Site)))+scale_x_reverse(position = "bottom",limits=c(250,0))+theme_classic()+geom_line()+coord_flip(xlim=c(250,0), ylim=c(50,80), clip="off")+scale_y_continuous(position = "right")+geom_point(shape=15)+labs(x="Depth (m)",y=expression("TOC ("*mu*"M)"))+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))+scale_color_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=toc,aes(x=Depth,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)
tocP

library(patchwork)
Den+Chl+Tmp+Oxy+n+tocP+plot_layout(guides = "collect",nrow=2) #+plot_annotation(tag_levels=list(c("b","c","d","e","f","g","h")))
ggsave("Depth_v_CTD.pdf",width=12,height=6)
