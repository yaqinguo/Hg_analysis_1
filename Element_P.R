library(tidyverse)
library(RColorBrewer)
library(ggbreak)
library(multcompView)
library(ggthemes)
library(car)
library(ggtext)
library(emmeans)
element_P<- read_csv("Element_P.csv",col_types = "fcffddddddddddddddd")
table(element_P$AMF,element_P$Hg_C) #to check data is balanced or not
element_P %>% distinct(AMF)
element_P %>% distinct(Hg_C)
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
##P concentration in the root
leveneTest(Root_Concentration ~ AMF*Hg_C,data=element_P) #p=0.5
element_P_25 <- element_P %>%
  select(c(1:4,11)) %>%
  filter(Hg_C=="25")
aov_element_P_25 <- aov(Root_Concentration ~ AMF, data = element_P_25)
Anova(aov_element_P_25,type="III")

element_P_Control <- element_P %>%
  select(c(1:4,11)) %>%
  filter(AMF=="Control")
aov_element_P_Control <- aov(Root_Concentration ~ Hg_C, data = element_P_Control)
Anova(aov_element_P_Control,type="III")

annova_root_concentration <- aov(Root_Concentration ~ AMF*Hg_C,data=element_P)
summary(annova_root_concentration) #only AMF (p=3.376e-10 ***)
shapiro.test(residuals(annova_root_concentration)) #p=0.57
tukey_rootconcentration<- TukeyHSD(annova_root_concentration)
multcompLetters4(annova_root_concentration,tukey_rootconcentration)
# tukey_rootconcentration <- emmeans(annova_root_concentration,specs = pairwise ~ AMF+Hg_C) %>% 
#   multcomp::cld(Letters=letters)
# tukey_rootconcentration
root_concentration <- data_summary(element_P,varname = "Root_Concentration",
                           groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
root_concentration$AMF <- as.factor(root_concentration$AMF)
root_concentration$Hg_C <- as.factor(root_concentration$Hg_C)
write_csv(root_concentration,"P_root_concentration_summary.csv")
y_expression.root <- expression(Root ~ P ~ Concnetration ~ (µg ~ g^-1))
P_rootconcentration <- ggplot(root_concentration,aes(x=Hg_C,y=Root_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Root_Concentration-sd,ymax=Root_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Root_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(axis.title.x = element_blank(),
    legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y=y_expression.root)+
  geom_richtext(x=2,y=3500,label="*P*(AMF)<0.001<br>*P*(Hg)=0.86<br>*P*(AMFXHg)=0.91",fill=NA,label.color=NA,size=3)+
  scale_y_continuous(expand = c(0,0),limits = c(0,3800))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  geom_text(x=3.5,y=3750,label="(B)",size=3)+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))
P_rootconcentration
ggsave("P_rootconcentration_P_bw_240322.tiff",width = 5,heigh=4)
##root content
leveneTest(Root_Content ~ AMF*Hg_C,data=element_P) #p=0.66
shapiro.test(element_P$Root_Content) #p=0.01 
annova_root_content <- lm(Root_Content ~ AMF+Hg_C,data=element_P)
Anova(annova_root_content,type="III") #only AMF (p=7.203e-10 ***)
shapiro.test(residuals(annova_root_content)) #p=0.54
tukey_rootcontent <- emmeans(annova_root_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_rootcontent
root_content <- data_summary(element_P,varname = "Root_Content",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a", "a","a","b","b","b"))
root_content$AMF <- as.factor(root_content$AMF)
root_content$Hg_C <- as.factor(root_content$Hg_C)
write_csv(root_content,"P_root_content_summary.csv")
P_rootcontent <- ggplot(root_content,aes(x=Hg_C,y=Root_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Root_Content-sd,ymax=Root_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Root_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Root P Content (µg)")+
  geom_text(x=3.3, y=1700, label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1800))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_rootcontent
ggsave("P_rootcontent_P.tiff",width = 5,heigh=4)
##P concentration in the leaf
leveneTest(Leaf_Concentration ~ AMF*Hg_C,data=element_P) #p=0.7
shapiro.test(element_P$Leaf_Concentration) #p=0.07
annova_leaf_concentration <- lm(Leaf_Concentration ~ AMF+Hg_C,data=element_P)
Anova(annova_leaf_concentration,type="III") #only AMF (p=3.347e-09 ***)
shapiro.test(residuals(annova_leaf_concentration)) #p=0.11
tukey_leafconcentration <- emmeans(annova_leaf_concentration,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
leaf_concentration <- data_summary(element_P,varname = "Leaf_Concentration",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
leaf_concentration$AMF <- as.factor(leaf_concentration$AMF)
leaf_concentration$Hg_C <- as.factor(leaf_concentration$Hg_C)
write_csv(leaf_concentration,"P_leaf_concentration_summary.csv")
P_leafconcentration <- ggplot(leaf_concentration,aes(x=Hg_C,y=Leaf_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Leaf_Concentration-sd,ymax=Leaf_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Leaf_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Leaf P Concentration (mg/kg)")+
  geom_text(x=3.2, y=2050, label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,2100))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_leafconcentration
ggsave("P_leafconcentration_P.tiff",width = 5,heigh=4)
##P content in the leaf
leveneTest(Leaf_Content ~ AMF*Hg_C,data=element_P) #p=0.6
shapiro.test(element_P$Leaf_Content) #p=0.3 
annova_leaf_content <- lm(Leaf_Content ~ AMF+Hg_C,data=element_P)
Anova(annova_leaf_content,type="III") #only AMF (p=3.384e-08 ***)
shapiro.test(residuals(annova_leaf_content))#p=0.3
tukey_leafcontent <- emmeans(annova_leaf_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_leafcontent
leaf_content <- data_summary(element_P,varname = "Leaf_Content",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
leaf_content$AMF <- as.factor(leaf_content$AMF)
leaf_content$Hg_C <- as.factor(leaf_content$Hg_C)
write_csv(leaf_content,"P_leaf_content_summary.csv")
P_leafcontent <- ggplot(leaf_content,aes(x=Hg_C,y=Leaf_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Leaf_Content-sd,ymax=Leaf_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Leaf_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Leaf P Content (µg)")+
  geom_text(x=3.2,y=1400,label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1500))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_leafcontent
ggsave("P_leafcontent_P.tiff",width = 5,heigh=4)
##P concentration in the stem
leveneTest(Stem_Concentration ~ AMF*Hg_C,data=element_P) #p=0.03
shapiro.test(element_P$Stem_Concentration) #p=0.0001014
annova_stem_concentration <- lm(Stem_Concentration ~ AMF+Hg_C,data=element_P)
Anova(annova_stem_concentration,type="III") # AMF (2.277e-09 ***)
shapiro.test(residuals(annova_stem_concentration))#p=0.01054
tukey_stemconcentration <- emmeans(annova_stem_concentration,specs = pairwise ~AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_stemconcentration
stem_concentration <- data_summary(element_P,varname = "Stem_Concentration",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
stem_concentration$AMF <- as.factor(stem_concentration$AMF)
stem_concentration$Hg_C <- as.factor(stem_concentration$Hg_C)
write_csv(stem_concentration,"P_stem_concentration_summary.csv")
P_stemconcentration <- ggplot(stem_concentration,aes(x=Hg_C,y=Stem_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Stem_Concentration-sd,ymax=Stem_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Stem_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Stem P Concentration (mg/kg)")+
  geom_text(x=3.2,y=4700,label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,4800))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_stemconcentration
ggsave("P_stemconcentration_P.tiff",width = 5,heigh=4)
##P content in the stem
leveneTest(Stem_Content ~ AMF*Hg_C,data=element_P) #p=0.09
shapiro.test(element_P$Stem_Content) #p=0.0001031
annova_stem_content <- lm(Stem_Content ~ AMF+Hg_C,data=element_P)
Anova(annova_stem_content,type="III") #only AMF (1.926e-09 ***)
shapiro.test(residuals(annova_stem_content)) #p=0.006013
tukey_stemcontent <- emmeans(annova_stem_content, specs = pairwise ~ AMF*Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_stemcontent
stem_content <- data_summary(element_P,varname = "Stem_Content",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
stem_content$AMF <- as.factor(stem_content$AMF)
stem_content$Hg_C <- as.factor(stem_content$Hg_C)
write_csv(stem_content,"P_stem_content_summary.csv")
P_stemcontent <- ggplot(stem_content,aes(x=Hg_C,y=Stem_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Stem_Content-sd,ymax=Stem_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Stem_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Stem P Content (µg)")+
  geom_text(x=3.2,y=2900,label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3000))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_stemcontent
ggsave("P_stemcontent_P.tiff",width = 5,heigh=4)
##P concentration in the above including stem and leaf
leveneTest(Above_Concentration ~ AMF*Hg_C,data=element_P) #p=0.05
shoot_element_P_25 <- element_P %>%
  select(c(1:4,17)) %>%
  filter(Hg_C=="25")
aov_shoot_P_25 <- aov(Above_Concentration ~ AMF,shoot_element_P_25)
Anova(aov_shoot_P_25,type="III")

shoot_element_P_Control <- element_P %>%
  select(c(1:4,17)) %>%
  filter(AMF=="Control")
aov_shoot_P_Control <- aov(Above_Concentration ~ Hg_C,shoot_element_P_Control)
Anova(aov_shoot_P_Control,type="III")


annova_above_concentration <- aov(Above_Concentration ~ AMF*Hg_C,data=element_P)
summary(annova_above_concentration) #  AMF (p=2.826e-10***) 
shapiro.test(residuals(annova_above_concentration)) #P=0.1411
tukey.above_concentration <- TukeyHSD(annova_above_concentration)
multcompLetters4(annova_above_concentration,tukey.above_concentration)

# tukey_aboveconcentration <- emmeans(annova_above_concentration,specs = pairwise ~AMF+Hg_C) %>%
#   multcomp::cld(Letters=letters)
# tukey_aboveconcentration
above_concentration <- data_summary(element_P,varname = "Above_Concentration",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
above_concentration$AMF <- as.factor(above_concentration$AMF)
above_concentration$Hg_C <- as.factor(above_concentration$Hg_C)
write_csv(above_concentration,"P_above_concentration_summary.csv")
y_expression.shoot <- expression(Shoot ~ P ~ Concnetration ~ (µg ~ g^-1))
P_aboveconcentration <- ggplot(above_concentration,aes(x=Hg_C,y=Above_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Above_Concentration-sd,ymax=Above_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Above_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(axis.title.x = element_blank(),
    legend.position = c(0.12,0.95),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y=y_expression.shoot)+
  geom_richtext(x=2,y=3500,label="*P*(AMF)<0.001<br>*P*(Hg)=0.52<br>*P*(AMFXHg)=0.35",fill=NA,label.color=NA,size=3)+
  scale_y_continuous(expand = c(0,0),limits = c(0,3800))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  geom_text(x=3.5,y=3750,label="(A)",size=3)+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))
P_aboveconcentration
ggsave("P_aboveconcentration_P_bw_240322.tiff",width = 5,heigh=4)
library(patchwork)
P_aboveconcentration|P_rootconcentration
ggsave("P_above_root_Pconcentration.tiff",width = 8,height = 6)
#P concentration in the above and root
P_above_root_concentration1 <- element_P %>%
  select(c(AMF,Hg_C,Above_Concentration,Root_Concentration))
aov_above <- aov(Above_Concentration ~ AMF*Hg_C,data = P_above_root_concentration1)
summary(aov_above)
tukey.above <- TukeyHSD(aov_above)
multcompLetters4(aov_above,tukey.above)

aov_root <- aov(Root_Concentration ~ AMF*Hg_C,data = P_above_root_concentration1)
summary(aov_root)
tukey.root <- TukeyHSD(aov_root)
multcompLetters4(aov_root,tukey.root)


P_above_root_concentration2 <- P_above_root_concentration1%>% 
  pivot_longer(cols = c(Above_Concentration,Root_Concentration),names_to="plant_P",
             values_to="P_concentration")%>%
  group_by(AMF,Hg_C,plant_P)%>%
  summarise(P_mean=mean(P_concentration),P_sd=sd(P_concentration))%>%ungroup()
P_above_root_concentration2$plant_P <- factor(P_above_root_concentration2$plant_P,
                                                levels = c("Above_Concentration","Root_Concentration"),
                                                labels = c("Shoot","Root"))

P_above_root_concentration2<-P_above_root_concentration2%>%
  mutate(label=c("a","a",
                 "a","a",
                 "a","a",
                 "b","b",
                 "b","b",
                 "b","b"))

labeldat <- data.frame(plant_P=c("Shoot","Root"),x=c(1.0,1.0),y=c(3500,3500),
                       label=c("*P*(AMF)<0.001<br>*P*(Hg)=0.52<br>*P*(AMFXHg)=0.35",
                               "*P*(AMF)<0.001<br>*P*(Hg)=0.86<br>*P*(AMFXHg)=0.91"))
dat_text <- data.frame(plant_P=c("Shoot","Root"),x=c(3.4,3.4),y=c(3920,3920),label=c("(A)","(B)"))

ggplot(P_above_root_concentration2,aes(x=Hg_C,y=P_mean,fill=AMF))+
  geom_errorbar(aes(x=Hg_C,ymin=P_mean-P_sd,ymax=P_mean+P_sd),position = position_dodge(0.7),width=0.2)+
  geom_bar(stat = "identity",position = position_dodge(0.7),color="black",width = 0.6)+
  geom_text(aes(x=Hg_C,y=P_mean+P_sd,label=label),vjust=-0.3,position = position_dodge(0.7))+
  facet_grid(~factor(plant_P,levels=c("Shoot","Root")))+
  theme_bw()+
  labs(y=P ~ Concentration ~ (µg ~ g ^ -1))+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"),
        panel.spacing.x = unit(0,"line"),panel.background = element_blank(),panel.grid = element_blank(),
        strip.text.x = element_text(face = "bold"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,4000))+
  scale_fill_manual(values = c("white","grey80"))+
  geom_text(data=dat_text, aes(x=x,y=y,label=label,fill=NULL))+
  geom_richtext(data=labeldat, aes(x=x,y=y,label=label,fill=NULL),
                fill=NA,label.color=NA,hjust=0,vjust=0,size=3,family = "Arial",fontface="plain")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))
ggsave("P_P_concentration_above_root_250322.tiff",width = 6,height = 6)


##P content in the above including stem and leaf
leveneTest(Above_Content ~ AMF*Hg_C,data=element_P) #p=0.07
shapiro.test(element_P$Above_Content) #p=0.0005423 
annova_above_content <- lm(Above_Content ~ AMF+Hg_C,data=element_P)
Anova(annova_above_content,type="III") #AMF (p=3.192e-10 ***)
shapiro.test(residuals(annova_above_content))#P=0.05047
tukey_abovecontent <- emmeans(annova_above_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_abovecontent
above_content <- data_summary(element_P,varname = "Above_Content",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
above_content$AMF <- as.factor(above_content$AMF)
above_content$Hg_C <- as.factor(above_content$Hg_C)
write_csv(above_content,"P_above_content_summary.csv")
P_abovecontent <- ggplot(above_content,aes(x=Hg_C,y=Above_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Above_Content-sd,ymax=Above_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Above_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Shoot P Content (µg)")+
  geom_text(x=3.2,y=4100,label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,4200))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_abovecontent
ggsave("P_abovecontent_P.tiff",width = 5,heigh=4)

#Transfer factor
leveneTest(Transfer_Factor ~ AMF*Hg_C,data=element_P) #p=0.74
shapiro.test(element_P$Transfer_Factor) #p=0.96 
annova_TF <- lm(Transfer_Factor ~ AMF+Hg_C,data=element_P)
Anova(annova_TF,type="III") #none
shapiro.test(residuals(annova_TF))#P=0.88
tukey_TF <- emmeans(annova_TF,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_TF
TF <- data_summary(element_P,varname = "Transfer_Factor",
                              groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","a","a","a"))
TF$AMF <- as.factor(TF$AMF)
TF$Hg_C <- as.factor(TF$Hg_C)
write_csv(TF,"P_TF_summary.csv")
P_TF <- ggplot(TF,aes(x=Hg_C,y=Transfer_Factor,fill=AMF))+
  geom_errorbar(aes(ymin=Transfer_Factor-sd,ymax=Transfer_Factor+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Transfer_Factor+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="P Transfer Factor")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.5))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_TF
ggsave("P_TF_P.tiff",width = 5,heigh=4)
