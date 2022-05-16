library(tidyverse)
library(readxl)
library(RColorBrewer)
library(ggthemes)
library(multcompView)
library(car)
library(emmeans)
library(ggpubr)
library(ggtext)
Height_Biomass<- read_excel("Plantheight_Biomass.xlsx",
                            na="NA")
Height_Biomass$treatment=as.factor(Height_Biomass$treatment)
Height_Biomass$AMF_C=as.factor(Height_Biomass$AMF_C)
Height_Biomass$Hg_C=as.factor(Height_Biomass$Hg_C)
#Height_Biomass <- Height_Biomass %>%
  #filter(Hg_C=="0"|Hg_C=="25")
table(Height_Biomass$AMF_C,Height_Biomass$Hg_C) #to see balanced or unbalanced
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
colnames(Height_Biomass)
#root/shoot
root_shoot_ratio <- Height_Biomass %>%
  select(c(treatment,AMF_C,Hg_C,Root_dry,Shoot_dry)) %>%
  mutate(root_shoot_ratio=Root_dry/Shoot_dry)
leveneTest(root_shoot_ratio ~ AMF_C*Hg_C,data = root_shoot_ratio)
aov_root_shoot_ratio <- aov(root_shoot_ratio ~ AMF_C*Hg_C,data = root_shoot_ratio)
Anova(aov_root_shoot_ratio,type="III")
shapiro.test(residuals(aov_root_shoot_ratio))
tukey.root_shoot_ratio <- TukeyHSD(aov_root_shoot_ratio)
multcompLetters4(aov_root_shoot_ratio,tukey.root_shoot_ratio)

aov_ratio <- aov(root_shoot_ratio ~ treatment,data = root_shoot_ratio)
summary(aov_ratio)
tukey_ratio <- TukeyHSD(aov_ratio)
cld_ratio <- multcompLetters4(aov_ratio,tukey_ratio)
cld_ratio


root_shoot_ratio_25 <- root_shoot_ratio %>%
  filter(Hg_C=="25")
aov_root_shoot_ratio_25 <- aov(root_shoot_ratio ~ AMF_C,data = root_shoot_ratio_25)
Anova(aov_root_shoot_ratio_25,type="III")
root_shoot_ratio_AMF <- root_shoot_ratio %>%
  filter(AMF_C=="AMF")
aov_root_shoot_ratio_AMF <- aov(root_shoot_ratio ~ Hg_C,data = root_shoot_ratio_AMF)
Anova(aov_root_shoot_ratio_AMF,type="III")
tukey_root_shoot_ratio <- TukeyHSD(aov_root_shoot_ratio_AMF)
cld_root_shoot_ratio <- multcompLetters4(aov_root_shoot_ratio_AMF, tukey_root_shoot_ratio)
cld_root_shoot_ratio
root_shoot_ratio_data <- data_summary(root_shoot_ratio,varname = "root_shoot_ratio",
                           groupnames = c("AMF_C","Hg_C"))
root_shoot_ratio_data$AMF_C <- factor(root_shoot_ratio_data$AMF_C,levels = c("Control","AMF"))
root_shoot_ratio_data$Hg_C <- as.factor(root_shoot_ratio_data$Hg_C)
root_shoot_ratio_data <- root_shoot_ratio_data %>% mutate(cld=c("a","ab","b","ab","ab","ab"))

P_root_shoot_ratio <- ggplot(root_shoot_ratio_data,aes(x=Hg_C,y=root_shoot_ratio,fill=AMF_C))+
  geom_errorbar(aes(ymin=root_shoot_ratio-sd,ymax=root_shoot_ratio+sd),width=0.2,
                position = position_dodge(0.6))+
  geom_bar(position = position_dodge(0.6),stat = "identity",
            width = 0.5,color="black")+
  geom_richtext(x=2.0,y=0.59,label=c("*P*(AMF)=0.53<br>*P*(Hg)=0.01<br>*P*(AMFXHg)=0.14"),fill=NA,
                hjust=0,vjust=0,
                label.color=NA,size=3,family = "Arial",fontface="plain")+
  geom_text(aes(x=Hg_C,y=root_shoot_ratio+sd,label=cld),position=position_dodge(0.6),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,
    legend.position = "none",axis.title.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10, family = "Arial"),plot.margin = margin(0,0,0,0,"pt") 
  )+
  labs(y="Root-shoot ratio")+
  geom_text(x=3.4, y=0.72,label="(D)",size=3)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.75),position = "right")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values = c("white","grey80"))
P_root_shoot_ratio
ggsave("P_root_shoot_ratio_190322.tiff",width = 5,heigh=4)
#shoot/root
shoot_root_ratio <- Height_Biomass %>%
  select(c(AMF_C,Hg_C,Root_dry,Shoot_dry)) %>%
  mutate(shoot_root_ratio=Shoot_dry/Root_dry)
leveneTest(shoot_root_ratio ~ AMF_C*Hg_C,data = shoot_root_ratio)
aov_shoot_root_ratio <- aov(shoot_root_ratio ~ AMF_C+Hg_C,data = shoot_root_ratio)
Anova(aov_shoot_root_ratio,type="III")
shapiro.test(residuals(aov_shoot_root_ratio))
tukey_shoot_root_ratio <- TukeyHSD(aov_shoot_root_ratio)
tukey_shoot_root_ratio
#hers is to know Hg increase or decrease the value
shoot_root_ratio_Hg <- data_summary(shoot_root_ratio,varname = "shoot_root_ratio",
                                      groupnames = "Hg_C")
shoot_root_ratio_25 <- shoot_root_ratio %>%
  filter(Hg_C=="25")
aov_shoot_root_ratio_25 <- aov(shoot_root_ratio ~ AMF_C,data = shoot_root_ratio_25)
Anova(aov_shoot_root_ratio_25,type="III")
shoot_root_ratio_AMF <- shoot_root_ratio %>%
  filter(AMF_C=="AMF")
aov_shoot_root_ratio_AMF <- aov(shoot_root_ratio ~ Hg_C,data = shoot_root_ratio_AMF)
Anova(aov_shoot_root_ratio_AMF,type="III")
tukey_shoot_root_ratio <- TukeyHSD(aov_shoot_root_ratio_AMF)
cld_shoot_root_ratio <- multcompLetters4(aov_shoot_root_ratio_AMF, tukey_shoot_root_ratio)
cld_shoot_root_ratio
shoot_root_ratio_data <- data_summary(shoot_root_ratio,varname = "shoot_root_ratio",
                                      groupnames = c("AMF_C","Hg_C"))
shoot_root_ratio_data$AMF_C <- factor(shoot_root_ratio_data$AMF_C,levels = c("Control","AMF"))
shoot_root_ratio_data$Hg_C <- as.factor(shoot_root_ratio_data$Hg_C)
shoot_root_ratio_data <- shoot_root_ratio_data %>% mutate(cld=c("a","a","a","a","a","a"))
write.csv(shoot_root_ratio_data,"shoot_root_ratio_data_summary.csv")
P_shoot_root_ratio <- ggplot(shoot_root_ratio_data,aes(x=Hg_C,y=shoot_root_ratio,fill=AMF_C))+
  geom_errorbar(aes(ymin=shoot_root_ratio-sd,ymax=shoot_root_ratio+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  #geom_text(aes(x=Hg_C,y=shoot_root_ratio+sd,label=cld),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
        legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10, family = "Arial"),plot.margin = margin(0,0,0,0,"pt")
  )+
  labs(y="shoot_root_ratio")+
  #geom_text(x=3.2, y=2.4,label="Hg*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.5),position = "right")+
  scale_x_discrete(breaks=c("0","25","50"),labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))
P_shoot_root_ratio
ggsave("P_shoot_root_ratio_bw.tiff",width = 5,heigh=4)
#root and shoot total dry weight
root_shoot_total <- Height_Biomass %>%
  select(c(treatment,AMF_C,Hg_C,Root_dry,Shoot_dry)) %>%
  mutate(root_shoot_total=Root_dry+Shoot_dry)
write.csv(root_shoot_total,file="root_shoot_total.csv")
leveneTest(root_shoot_total ~ AMF_C*Hg_C,data=root_shoot_total)
aov_root_shoot_total <- aov(root_shoot_total ~ AMF_C*Hg_C,data=root_shoot_total)
Anova(aov_root_shoot_total,type="III")
shapiro.test(residuals(aov_root_shoot_total))
tukey.root_shoot_total <- TukeyHSD(aov_root_shoot_total)
multcompLetters4(aov_root_shoot_total,tukey.root_shoot_total)


leveneTest(root_shoot_total~treatment,data = root_shoot_total)
aov_total <- aov(root_shoot_total~treatment,data = root_shoot_total)
summary(aov_total)
shapiro.test(residuals(aov_total))
tukey_total <- TukeyHSD(aov_total)
cld_total <- multcompLetters4(aov_total,tukey_total)
cld_total


root_shoot_total_25 <- root_shoot_total %>%
  filter(Hg_C=="25")
leveneTest(root_shoot_total ~ AMF_C,data = root_shoot_total_25)
aov_root_shoot_total_25 <- aov(root_shoot_total ~ AMF_C,data = root_shoot_total_25)
shapiro.test(residuals(aov_root_shoot_total_25))
Anova(aov_root_shoot_total_25,type= "III")


root_shoot_total_Control <- root_shoot_total %>%
  filter(AMF_C=="Control")
aov_root_shoot_total_Control <- aov(root_shoot_total ~ Hg_C,data = root_shoot_total_Control)
Anova(aov_root_shoot_total_Control,type= "III")
root_shoot_total_data <- data_summary(root_shoot_total,varname = "root_shoot_total",
                                      groupnames = c("AMF_C","Hg_C"))
root_shoot_total_data$AMF_C <- factor(root_shoot_total_data$AMF_C,levels = c("Control","AMF"))
root_shoot_total_data$Hg_C <- as.factor(root_shoot_total_data$Hg_C)
root_shoot_total_data <- root_shoot_total_data %>% mutate(cld=c("ab","ab","b","a","ab","ab"))
write.csv(root_shoot_total_data,"root_shoot_total_data_summary.csv")
P_root_shoot_total <- ggplot(root_shoot_total_data,aes(x=Hg_C,y=root_shoot_total,fill=AMF_C))+
  geom_errorbar(aes(ymin=root_shoot_total-sd,ymax=root_shoot_total+sd),width=0.2,
                position = position_dodge(0.6))+
  geom_bar(position = position_dodge(0.6),stat = "identity",
           width = 0.5,color="black")+
  geom_richtext(x=2.0,y=2.1,label=c("*P*(AMF)<0.01<br>*P*(Hg)=0.88<br>*P*(AMFXHg)=0.54"),fill=NA,
                hjust=0,vjust=0,
                label.color=NA,size=3,family = "Arial",fontface="plain")+
  geom_text(aes(x=Hg_C,y=root_shoot_total+sd,label=cld),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme( axis.title.x = element_blank(),aspect.ratio = 1,
    legend.position ="none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10, family = "Arial"),plot.margin = margin(0,0,0,0,"pt") 
  )+
  labs(y="Total dry biomass (g)")+
  geom_text(x=3.4,y=2.6,label="(C)",size=3)+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.7))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values = c("white","grey80"))
P_root_shoot_total
ggsave("P_root_shoot_totaldryweight_190322.tiff",width = 5,heigh=4)

MBR <- read_csv("root_shoot_total.csv",col_types = "cffddd")
MBR <- MBR %>%
  select(c(AMF_C,Hg_C,MBR))%>%
  filter(AMF_C=="AMF")
leveneTest(MBR ~ Hg_C,data = MBR)
aov_MBR <- aov(MBR ~ Hg_C,data = MBR)
Anova(aov_MBR,type="III")
tukey_MBR <- TukeyHSD(aov_MBR)
cld_MBR <- multcompLetters4(aov_MBR, tukey_MBR)
cld_MBR
MBR_data <- data_summary(MBR,varname = "MBR",
                           groupnames = c("Hg_C"))
MBR_data$Hg_C <- as.factor(MBR_data$Hg_C)
MBR_data <- MBR_data %>%mutate(cld=c("b","a","a"))
print(MBR_data)
p_MBR <- ggplot(MBR_data)+
  geom_errorbar(aes(x=Hg_C,ymin=MBR-sd,ymax=MBR+sd),width=0.2,
                 position = position_dodge(0.7))+
  geom_bar(aes(x=Hg_C,y=MBR),position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black",fill="grey70")+
  geom_text(aes(x=Hg_C,y=MBR-sd-0.3,label=cld),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10, family = "Arial") 
  )+
  labs(x="Hg_Concentration (mg/kg)",y="Mycorrhizal biomass response (%)")+
  geom_text(x=3.2, y=26,label="Hg*")+
  scale_y_continuous(expand = c(0,0),limits = c(-8,0))
p_MBR
ggsave("p_MBR.tiff",width = 5,heigh=4)


##height
leveneTest(PlantHeight ~ Hg_C*AMF_C,Height_Biomass) #function from car package to test data homogeneity,p=0.4
shapiro.test(Height_Biomass$PlantHeight) #test data normality, p=0.65
##meet anov assumption, then doing anova analysis
plant_height_0 <- Height_Biomass %>%
  select(c(1:7)) %>%
  filter(Hg_C=="0")
avo_plantheight_0<- aov(PlantHeight ~ AMF_C,plant_height_0)
Anova(avo_plantheight_0,type="III")
plant_height_AMF <- Height_Biomass %>%
  select(c(1:7)) %>%
  filter(AMF_C=="AMF")
avo_plantheight_AMF <- aov(PlantHeight ~ Hg_C,plant_height_AMF)
Anova(avo_plantheight_AMF,type="III")
avo_plantheight <- aov(PlantHeight ~ AMF_C+Hg_C,Height_Biomass)
Anova(avo_plantheight,type="III") #function from car package, because my daya is unbanlanced
shapiro.test(residuals(avo_plantheight))
tukey_plantheight<- TukeyHSD(avo_plantheight)
#here is to get letters 
cld_plantheight <- multcompLetters4(avo_plantheight, tukey_plantheight)
MPN_height <- data_summary(Height_Biomass,varname = "PlantHeight",
                           groupnames = c("AMF_C","Hg_C"))
MPN_height$AMF_C <- factor(MPN_height$AMF_C,levels = c("Control","AMF"))
MPN_height$Hg_C <- as.factor(MPN_height$Hg_C)
MPN_height <- MPN_height %>%mutate(cld=c("ab","ab","b","ab","ab","a"))
print(MPN_height )
write.csv(MPN_height,"MPN_plantheight_summary.csv")
P_height <- ggplot(MPN_height,aes(x=Hg_C,y=PlantHeight,fill=AMF_C))+
  geom_errorbar(aes(ymin=PlantHeight-sd,ymax=PlantHeight+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=PlantHeight+sd,label=cld),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10, family = "Arial") 
        )+
  labs(x="Hg_Concentration (mg/kg)",y="Plant_Height (cm)")+
  geom_text(x=3.2, y=26,label="AMF*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,27))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_height
ggsave("P_height_new.tiff",width = 5,heigh=4)
##dry shoot biomass
leveneTest(Shoot_dry ~ Hg_C*AMF_C,Height_Biomass) #function from car package to test data homogeneity,p=0.5
#hist(Height_Biomass$Shoot_dry) this is to confirm whether data is normal or not
#above, my data is not normally distrubuted,then i go for the following method 
# shoot_biomass <- glm(Shoot_dry ~ Hg_C+AMF_C,Height_Biomass,family = "quasipoisson" )
# summary(shoot_biomass) #AMF have a significant effect (p=0.002)
leveneTest(Shoot_dry ~ treatment,data=Height_Biomass)
aov_shoot <- aov(Shoot_dry ~ treatment,data=Height_Biomass)
summary(aov_shoot)
shapiro.test(residuals(aov_shoot))
tukey_shoot <- TukeyHSD(aov_shoot)
cld_shoot <- multcompLetters4(aov_shoot,tukey_shoot)
cld_shoot 

shoot_dry_0 <- Height_Biomass %>%
  select(c(1:6,14)) %>%
  filter(Hg_C=="0")
avo_dryshoot_0 <- aov(Shoot_dry~ AMF_C,data=shoot_dry_0)
leveneTest(Shoot_dry~ AMF_C,data=shoot_dry_0)
shapiro.test(residuals(avo_dryshoot_0))
Anova(avo_dryshoot_0,type="III")
#here only under 0 is significant
tukey_dryshoot_0<- TukeyHSD(avo_dryshoot_0)
cld_dryshoot_0 <- multcompLetters4(avo_dryshoot_0, tukey_dryshoot_0)
#
shoot_dry_Control <- Height_Biomass %>%
  select(c(1:6,14)) %>%
  filter(AMF_C=="Control")
avo_dryshoot_Control <- aov(Shoot_dry~ Hg_C,data=shoot_dry_Control)
Anova(avo_dryshoot_Control,type="III")
#here is main effect
avo_dryshoot <- aov(Shoot_dry ~ Hg_C*AMF_C,data=Height_Biomass)
Anova(avo_dryshoot,type="III")# AMF 0.002463 ** 
shapiro.test(residuals(avo_dryshoot))#p-value = 0.382
tukey_dryshoot<- TukeyHSD(avo_dryshoot)
#here is to get letters 
cld_dryshoot <- multcompLetters4(avo_dryshoot, tukey_dryshoot)
MPN_dryshoot <- data_summary(Height_Biomass,varname = "Shoot_dry",
                            groupnames = c("AMF_C","Hg_C"))
MPN_dryshoot$AMF_C=factor(MPN_dryshoot$AMF_C,levels = c("Control","AMF"))
MPN_dryshoot$Hg_C=as.factor(MPN_dryshoot$Hg_C)
MPN_dryshoot <- MPN_dryshoot%>%mutate(cld=c("a","a","a","a","a","a"))
write.csv(MPN_dryshoot,"MPN_dryshoot_summary.csv")
P_abovedrybiomass <- ggplot(MPN_dryshoot,aes(x=Hg_C,y=Shoot_dry,fill=AMF_C))+
  geom_errorbar(aes(ymin=Shoot_dry-sd,ymax=Shoot_dry+sd),width=0.2,
                position = position_dodge(0.6))+
  geom_bar(position = position_dodge(0.6),stat = "identity",
           width = 0.5,color="black")+
  geom_richtext(x=2.0,y=2.1,label=c("*P*(AMF)<0.01<br>*P*(Hg)=0.91<br>*P*(AMFXHg)=0.42"),fill=NA,
                hjust=0,vjust=0,
                label.color=NA,size=3,family = "Arial",fontface="plain")+
  geom_text(aes(x=Hg_C,y=Shoot_dry+sd,label=cld),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position ="none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10,family = "Arial"),plot.margin = margin(0,0,0,0,"pt"),
    axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  labs(y="Dry shoot biomass (g)")+
  geom_text(x=3.4, y=2.6,label="(A)",size=3)+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.7))+
  scale_fill_manual(values = c("white","grey80"))
P_abovedrybiomass
ggsave("Dry_shoot_biomass_190322.tiff",width = 5,height =4)

##dry root biomass
leveneTest(Root_dry ~ Hg_C*AMF_C,Height_Biomass) #function from car package to test data homogeneity,p=0.5

aov_root <- aov(Root_dry ~ treatment,data = Height_Biomass)
summary(aov_root)
tukey_root <- TukeyHSD(aov_root)
cld_root <- multcompLetters4(aov_root,tukey_root)
cld_root
root_dry_50 <- Height_Biomass %>%
  select(c(1:6,13)) %>%
  filter(Hg_C=="50")
aov_root_50 <- aov(Root_dry ~ AMF_C, data=root_dry_50)
Anova(aov_root_50,type="III")
#t.test(Root_dry ~ AMF_C, data=root_dry_25, var.equal = TRUE)

# root_dry_50 <- aov(Root_dry~ AMF_C,data=root_dry_50)
# Anova(root_dry_50,type="III")
root_dry_Control <- Height_Biomass %>%
  select(c(1:6,13)) %>%
  filter(AMF_C=="Control")
root_dry_Control <- aov(Root_dry~ Hg_C,data=root_dry_Control)
Anova(root_dry_Control,type="III")
shapiro.test(residuals(root_dry_AMF))
tukey_dryroot_AMF<- TukeyHSD(root_dry_AMF)
cld_dryroot <- multcompLetters4(root_dry_AMF, tukey_dryroot_AMF)
#here is for main effect
avo_dryroot <- aov(Root_dry ~ Hg_C*AMF_C,Height_Biomass)
summary(avo_dryroot)
Anova(avo_dryroot,type = "III")#Both AMF and Hg have a significant effect on dry root biomass 
shapiro.test(residuals(avo_dryroot))#p=0.982
tukey_dryroot<- TukeyHSD(avo_dryroot)
tukey.cld_dryroot <- emmeans(avo_dryroot,specs = pairwise ~ AMF_C+Hg_C) %>% 
   multcomp::cld(Letters=letters)

cld_Root <- multcompLetters4(avo_dryroot, tukey_dryroot)
cld_Root


MPN_dryroot <- data_summary(Height_Biomass,varname = "Root_dry",
                            groupnames = c("AMF_C","Hg_C"))
MPN_dryroot$AMF_C=factor(MPN_dryroot$AMF_C,levels = c("Control","AMF"))
MPN_dryroot$Hg_C=as.factor(MPN_dryroot$Hg_C)
MPN_dryroot <- MPN_dryroot%>%mutate(cld=c("ab","b","b","a","ab","ab"))
write.csv(MPN_dryroot,"MPN_dryroot_summary.csv")
P_dryroot <- ggplot(MPN_dryroot,aes(x=Hg_C,y=Root_dry,fill=AMF_C))+
  geom_errorbar(aes(ymin=Root_dry-sd,ymax=Root_dry+sd),width=0.2,
                position = position_dodge(0.6))+
  geom_bar(position = position_dodge(0.6),stat = "identity",
           width = 0.5,color="black")+
  geom_richtext(x=2.0,y=2.1,label=c("*P*(AMF)<0.01<br>*P*(Hg)=0.03<br>*P*(AMFXHg)=0.54"),fill=NA,
                hjust=0,vjust=0,
                label.color=NA,size=3,family = "Arial",fontface="plain")+
  geom_text(aes(x=Hg_C,y=Root_dry+sd,label=cld),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(axis.title.x = element_blank(),aspect.ratio = 1,
    legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10,family = "Arial"),plot.margin =margin(0,0,0,0,"pt"),
    axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  labs(y="Dry root biomass (g)")+
  geom_text(x=3.4,y=2.6,label="(B)",size=3)+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.7),position = "right")+
  scale_fill_manual(values = c("white","grey80"))
P_dryroot
ggsave("Dry_root_biomass_190322.tiff",width = 5,heigh=4)

library(patchwork)
# col1 <- P_abovedrybiomass/P_root_shoot_total
# col2 <- P_dryroot/P_root_shoot_ratio
# col1|col2


col1 <- (P_abovedrybiomass|P_dryroot)
col2 <- (P_root_shoot_total|P_root_shoot_ratio)
col1/col2
ggsave("shoot_root_total_ratio.tiff",width = 6,height = 6)













##stem dry
leveneTest(Stem_dry ~ AMF_C*Hg_C,data = Height_Biomass ) #p=0.5


stem_dry_50 <- Height_Biomass %>%
  select(c(1:6,Stem_dry)) %>%
  filter(Hg_C=="50")

t.test(Stem_dry ~ AMF_C, data=stem_dry_50, var.equal = TRUE)


aov_drystem <- aov(Stem_dry ~ Hg_C+AMF_C,Height_Biomass)
Anova(aov_drystem,type = "III") #AMF (p=0.003019 **)
shapiro.test(residuals(aov_drystem))#p-value = 0.00794

stem_dry_25 <- Height_Biomass %>%
  select(c(1:6,11)) %>%
  filter(Hg_C=="25")
stem_dry_25 <- aov(Stem_dry~ AMF_C,data=stem_dry_25)
Anova(stem_dry_25,type="III")
stem_dry_Control <- Height_Biomass %>%
  select(c(1:6,11)) %>%
  filter(AMF_C=="Control")
stem_dry_Control <- aov(Stem_dry~ Hg_C,data=stem_dry_Control)
Anova(stem_dry_Control,type="III")
shapiro.test(residuals(root_dry_AMF))
tukey_dryroot_AMF<- TukeyHSD(root_dry_AMF)
cld_dryroot <- multcompLetters4(root_dry_AMF, tukey_dryroot_AMF)
tukey.cld_drystem <- emmeans(aov_drystem,specs = pairwise ~ AMF_C+Hg_C) %>% 
  multcomp::cld(Letters=letters)
MPN_drystem <- data_summary(Height_Biomass,varname = "Stem_dry",
                     groupnames = c("AMF_C","Hg_C"))
MPN_drystem <- MPN_drystem %>% mutate(cld=c("a","a","a","b","b","b"))
write.csv(MPN_drystem,"stem_dry_summary.csv")
p_stem_dry <- ggplot(MPN_drystem,aes(x=Hg_C,y=Stem_dry,fill=AMF_C))+
  geom_errorbar(aes(ymin=Stem_dry-sd,ymax=Stem_dry+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width=0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Stem_dry+sd,label=cld),position = position_dodge(0.7),vjust=-0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10,family = "Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Dry Stem Biomass (g)")+
  geom_text(x=3.2,y=0.95,label="AMF**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  scale_fill_manual(values = c('#999999','#E69F00'),labels=c("Control","AMF"))
p_stem_dry
ggsave("Dry_stem_biomass.tiff",width = 5,heigh=4)
##leaf dry
leveneTest(Leaf_dry ~ AMF_C*Hg_C,data = Height_Biomass) #p=0.2
shapiro.test(Height_Biomass$Leaf_dry)#p-value = 0.02234
leaf_dry_biomass <- aov(Leaf_dry ~ Hg_C+AMF_C,Height_Biomass)
Anova(leaf_dry_biomass,type = "III")#AMF p=0.01216 * 
shapiro.test(residuals(leaf_dry_biomass))#p=0.7182

leaf_dry_25 <- Height_Biomass %>%
  select(c(1:6,Leaf_dry)) %>%
  filter(Hg_C=="25")

t.test(Leaf_dry ~ AMF_C, data=leaf_dry_25, var.equal = TRUE)











leaf_dry <- data_summary(Height_Biomass,varname = "Leaf_dry",
                         groupnames = c("AMF_C","Hg_C"))
leaf_dry <- leaf_dry %>%mutate(cld=c("a","a","a","b","b","b"))
write.csv(leaf_dry,"leaf_dry_summary.csv")
p_leaf_dry <- ggplot(leaf_dry,aes(x=Hg_C,y=Leaf_dry,fill=AMF_C))+
  geom_errorbar(aes(ymin=Leaf_dry-sd,ymax=Leaf_dry+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width=0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Leaf_dry+sd,label=cld),position = position_dodge(0.7),vjust=-0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size = 10,family = "Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Dry Leaf Biomass (g)")+
  geom_text(x=3.2,y=0.95,label="AMF*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  scale_fill_manual(values = c('#999999','#E69F00'),labels=c("Control","AMF"))
p_leaf_dry
ggsave("Dry_leaf_biomass_Hg.tiff",width = 5,heigh=4)
#create stack plot
Height_Biomass<- read_excel("Plantheight_Biomass.xlsx",
                            na="NA")
colnames(Height_Biomass)
Height_Biomass$AMF_C=factor(Height_Biomass$AMF_C,levels = c("Control","AMF"))
Height_Biomass$Hg_C=as.factor(Height_Biomass$Hg_C)
data_0_leaf_stem <- Height_Biomass %>%
  select(c("AMF_C","Hg_C","Leaf_dry","Stem_dry"))%>%
  filter(Hg_C=="0") %>%
  pivot_longer(cols = c("Leaf_dry","Stem_dry"),names_to="Plant_part",
               values_to="Weight")
data_0_leaf_stem$Plant_part <- as.factor(data_0_leaf_stem$Plant_part)
data_0_leaf_stem$AMF_C <- factor(data_0_leaf_stem$AMF_C,levels = c("Control","AMF"))
p_leaf_stem <- ggbarplot(data_0_leaf_stem,x="AMF_C",y="Weight",add="mean_sd",
               fill ="Plant_part",width = 0.5,
               palette=c("forestgreen","gold4"))
p_leaf_stem+geom_text(x=1,y=2,label="a")+geom_text(x=2,y=2,label="b")

ggsave("p_leaf_stem.tiff",width = 5,height = 4)
#p_leaf_stem+rremove("x.axis")+rremove("xlab")+rremove("x.text")+
  #rremove("x.ticks")+rremove("legend")



#test
mm_above <- aov(Weight ~ AMF_C,data=data_0_leaf_stem)
car::Anova(mm_above,type="III")


data_0_root <- Height_Biomass %>%
  select(c("AMF_C","Hg_C","Root_dry"))%>%
  filter(Hg_C=="0")
p_root <- ggbarplot(data_0_root,x="AMF_C",y="Root_dry",add="mean_sd",
                           width = 0.5,
                           fill = "gold4")
p_root
ggsave("p_root.tiff",width = 5,height = 4)

ggplot(aes(x=AMF_C))+
  geom_bar(data=data_0_leaf_stem,aes(y=Weight,fill=Plant_part))+
  geom_bar(data=data_0_root,aes(y=-Weight,fill=Plant_part))






data_0 <- data_0 %>%
  group_by(AMF_C,Plant_part)%>%
  mutate(bar=cumsum(Weight))%>%
  summarise(mean=mean(Weight),sd=sd(Weight))
  mutate(Weight=cumsum(Weight))


  ggplot(aes(x=AMF_C,y=Weight,fill=Plant_part))+
  geom_bar(width = 0.6,position = position_stack(vjust = 1),
           stat="identity")+
  geom_errorbar(inherit.aes = F,
                aes(x=AMF_C,ymin=Weight-sd,ymax=Weight+sd),width=0.2,
                position = "identity")
  


data_0_leaf_stem_root <- Height_Biomass %>%
  mutate(Root_dry=-Root_dry)%>%
    select(c("AMF_C","Hg_C","Leaf_dry","Stem_dry","Root_dry"))%>%
    filter(Hg_C=="0") %>%
    pivot_longer(cols = c("Leaf_dry","Stem_dry","Root_dry"),names_to="Plant_part",
                 values_to="Weight") %>%
  mutate(part=c("Above","Above","Below","Above","Above","Below","Above","Above","Below",
                    "Above","Above","Below","Above","Above","Below","Above","Above","Below",
                    "Above","Above","Below","Above","Above","Below"))
data_0_leaf_stem_root
p_shoot_root <- ggplot(data_0_leaf_stem_root,aes(x=AMF_C,y=Weight,fill=Plant_part))+
  geom_bar(stat = "identity",position = position_stack(0.7),width = 0.5)+
  theme_classic()+
  geom_hline(yintercept=0)+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
 
p_shoot_root



data_all_leaf_stem_root <- Height_Biomass %>%
  mutate(Root_dry=-Root_dry)%>%
  select(c("AMF_C","Hg_C","Leaf_dry","Stem_dry","Root_dry"))%>%
  pivot_longer(cols = c("Leaf_dry","Stem_dry","Root_dry"),names_to="Plant_part",
               values_to="Weight")

p_all <- ggplot(data_all_leaf_stem_root,aes(x=AMF_C,y=Weight,fill=Plant_part))+
  facet_grid(~Hg_C,scales = "free")+
  geom_bar(stat = "identity",position = position_stack(0.7),width = 0.5,
           color=aes(AMF_C))+
  theme_classic()+
  theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),panel.spacing = unit(0,'lines'))+
  geom_hline(yintercept=0)
p_all











