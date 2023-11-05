### Figure S2 - Environmental Data vs. Density ###
### This script will show you have to make plots that show changes in environmental variables (e.g., temperature, oxygen, etc.) as a function of density. ###
### Last Updated: November 5, 2023 ###

# Load data and add empty points
ctd <- read.csv("Falkor_CTD_Oct2023.csv",header=TRUE)
ctd$Point <- ifelse(ctd$Depth==25 & ctd$Site==85,1,0)
ctd$Point <- ifelse(ctd$Depth==25 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==150 & ctd$Site==85,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==150 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==250 & ctd$Site==85,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==250 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==112 & ctd$Site==69,1,ctd$Point)
ctd$Point <- ifelse(ctd$Depth==122 & ctd$Site==85,1,ctd$Point)

# Chlorophyll
ctd$CHL2 <- (ctd$CHL-0.10896)/1.7389
ctd$Point <- ifelse(ctd$Point!=0,ctd$CHL2,NA)
denChl <- ggplot(ctd, aes(x=Density, y=CHL2, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom")+theme_classic()+labs(y = expression("Chlorophyll fluorescence (RFU)"), x=expression("Density (kg m"^{-3}*")"))+geom_line()+coord_flip(ylim=c(-0.05,1), clip="off")+scale_y_continuous(position = "right")+scale_color_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=ctd,aes(x=Density,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))
denChl


# Temperature
ctd$Point <- ifelse(!is.na(ctd$Point),ctd$TMP,NA)
denTmp <- ggplot(ctd, aes(x=Density, y=TMP, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom")+theme_classic()+labs(y = expression("Temperature (°C)"),x=expression("Density (kg m"^{-3}*")"))+geom_line()+coord_flip(ylim=c(10,25), clip="off")+scale_y_continuous(position = "right")+scale_color_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=ctd,aes(x=Density,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))
denTmp

# Oxygen
ctd$Point <- ifelse(!is.na(ctd$Point),ctd$OXY,NA)
denOxy <- ggplot(ctd, aes(x=Density, y=OXY, group=as.factor(Site), color=as.factor(Site))) +scale_x_reverse(position = "bottom")+theme_classic()+labs(y = expression("Oxygen ("*mu*"mol L"^{-1}*")"), x=expression("Density (kg m"^{-3}*")"))+geom_line()+coord_flip(ylim=c(190,230), clip="off")+scale_y_continuous(position = "right")+scale_color_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))+geom_point(data=ctd,aes(x=Density,y=Point,fill=as.factor(Site)),shape=21,color="black",size=3)+scale_fill_manual(name="Site",labels=c("Cyclonic","SLA Approx. 0","Anticyclonic"), values=c("blue","black","red"))
denOxy


denChl+denTmp+denOxy+plot_layout(guides = "collect")
ggsave("Den_v_CTD.pdf",width=10,height=4)
