library(tidyverse)
library(RColorBrewer)
library(ggbreak)
library(multcompView)
library(ggthemes)
library(car)
library(ggtext)
library(emmeans)
element_Zn<- read_csv("Element_Zn.csv",col_types = "fcfffddddddddddddddd")
table(element_Zn$AMF,element_Zn$Hg_C) #to check data is balanced or not
element_Zn %>% distinct(AMF)
element_Zn %>% distinct(Hg_C)
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
#
Zn_phytoextraction_efficiency <- element_Zn %>%
  select(c(3:4,Root_dry,Above_Content))%>%
  mutate(phytoextraction_efficiency=(Above_Content/Root_dry))
leveneTest(phytoextraction_efficiency ~ AMF*Hg_C,data=Zn_phytoextraction_efficiency)
aov_Zn_phytoextraction_efficiency <- aov(phytoextraction_efficiency ~ AMF+Hg_C,data=Zn_phytoextraction_efficiency)
Anova(aov_Zn_phytoextraction_efficiency,type="III")
shapiro.test(residuals(aov_Zn_phytoextraction_efficiency))
Zn_phytoextraction_efficiency_data <- data_summary(Zn_phytoextraction_efficiency,varname = "phytoextraction_efficiency",
                                                groupnames = c("AMF","Hg_C"))


p_Zn_phytoextraction_efficiency <- ggplot(Zn_phytoextraction_efficiency_data,
                                          aes(x=Hg_C,y=phytoextraction_efficiency,group=AMF,
                                              fill=AMF))+
  geom_line()+
  geom_point()

p_Zn_phytoextraction_efficiency

Zn_uptake_efficiency <- element_Zn %>%
  select(c(3:4,Root_Content,Above_Content,Root_dry))%>%
  mutate(Zn_uptake_efficiency=(Root_Content+Above_Content)/Root_dry)
leveneTest(Zn_uptake_efficiency ~ AMF*Hg_C,data=Zn_uptake_efficiency)
aov_Zn_uptake_efficiency <- aov(Zn_uptake_efficiency ~ AMF+Hg_C,data=Zn_uptake_efficiency)
Anova(aov_Zn_uptake_efficiency,type="III")
Zn_uptake_efficiency_data <- data_summary(Zn_uptake_efficiency,varname = "Zn_uptake_efficiency",
                                       groupnames = c("AMF","Hg_C")) 
Zn_uptake_efficiency <- ggplot(Zn_uptake_efficiency_data,aes(x=Hg_C,y=Zn_uptake_efficiency,group=AMF,color=AMF))+
  #geom_errorbar(aes(ymin=uptake_efficiency-sd,ymax=uptake_efficiency+sd),width=0.1)+
  geom_line()+
  geom_point(size=4,shape=21)+
  #geom_text(aes(x=Hg_C,y=uptake_efficiency+0.7,label=tukey))+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family = "Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Zn_uptake_efficiency")+
  #geom_text(x=3.2, y=38, label="Hg***, AMF*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,300))+
  scale_color_manual(values=c('#999999','#E69F00'))
Zn_uptake_efficiency
ggsave("Zn_uptake_efficiency.tiff",width = 5,heigh=4)

##concentration in the root
leveneTest(Root_Concentration ~ AMF*Hg_C,data=element_Zn) #p=0.03
 

root_Zn_50 <- element_Zn %>%
  select(c(1:4,Root_Concentration)) %>%
  filter(Hg_C=="50")
# aov_root_Zn_25 <- aov(Root_Concentration ~ AMF, data = root_Zn_25)
# Anova(aov_root_Zn_25,type="III")

t.test(Root_Concentration ~ AMF,data = root_Zn_50,var.equal = TRUE)

root_Zn_Control <- element_Zn %>%
  select(c(1:4,11)) %>%
  filter(AMF=="Control")
aov_root_Zn_Control <- aov(Root_Concentration ~ Hg_C, data = root_Zn_Control)
Anova(aov_root_Zn_Control,type="III")
annova_root_concentration <- aov(Root_Concentration ~ AMF+Hg_C,data=element_Zn)
Anova(annova_root_concentration,type="III") #only AMF (p=0.007655 **)
shapiro.test(residuals(annova_root_concentration)) #p=0.57
tukey_rootconcentration <- emmeans(annova_root_concentration,specs = pairwise ~ AMF+Hg_C) %>% 
  multcomp::cld(Letters=letters)
tukey_rootconcentration
root_concentration <- data_summary(element_Zn,varname = "Root_Concentration",
                           groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("ac","ab","ac","ac","c","ac"))
root_concentration$AMF <- as.factor(root_concentration$AMF)
root_concentration$Hg_C <- as.factor(root_concentration$Hg_C)
write_csv(root_concentration,"Zn_root_concentration_summary.csv")
y_expression <- expression(Root~Zn~Concentration~(mg~kg^-1))
P_rootconcentration <- ggplot(root_concentration,aes(x=Hg_C,y=Root_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Root_Concentration-sd,ymax=Root_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  #geom_text(aes(x=Hg_C,y=Root_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position ="none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y=y_expression)+
  geom_text(x=1, y=240, label="ns")+
  geom_text(x=2, y=240, label="*")+
  geom_text(x=3, y=240, label="ns")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,260))+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))
P_rootconcentration
ggsave("P_rootconcentration_Zn_new_bw.tiff",width = 5,heigh=4)
##root content
leveneTest(Root_Content ~ AMF*Hg_C,data=element_Zn) #p=0.08
shapiro.test(element_Zn$Root_Content) #p=0.01 
annova_root_content <- lm(Root_Content ~ AMF+Hg_C,data=element_Zn)
Anova(annova_root_content,type="III") #only AMF (p=0.002869 **)
shapiro.test(residuals(annova_root_content)) #p=0.42
tukey_rootcontent <- emmeans(annova_root_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_rootcontent
root_content <- data_summary(element_Zn,varname = "Root_Content",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("b","b","b","a", "a","a"))
root_content$AMF <- as.factor(root_content$AMF)
root_content$Hg_C <- as.factor(root_content$Hg_C)
write_csv(root_content,"Zn_root_content_summary.csv")
P_rootcontent <- ggplot(root_content,aes(x=Hg_C,y=Root_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Root_Content-sd,ymax=Root_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Root_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Root Zn Content (µg)")+
  geom_text(x=3.3, y=170, label="AMF**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,180))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_rootcontent
ggsave("P_rootcontent_Zn.tiff",width = 5,heigh=4)
## concentration in the leaf
leveneTest(Leaf_Concentration ~ AMF*Hg_C,data=element_Zn) #p=0.7


leaf_concentration_25 <- element_Zn%>%
  select(c(1:4,Leaf_Concentration))%>%
  filter(Hg_C=="25")

t.test(Leaf_Concentration ~ AMF, data=leaf_concentration_25,var.equal = TRUE)


annova_leaf_concentration <- lm(Leaf_Concentration ~ AMF+Hg_C,data=element_Zn)
Anova(annova_leaf_concentration,type="III") #only AMF (p=0.00209 **)
shapiro.test(residuals(annova_leaf_concentration)) #p=0.11
tukey_leafconcentration <- emmeans(annova_leaf_concentration,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_leafconcentration
leaf_concentration <- data_summary(element_Zn,varname = "Leaf_Concentration",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
leaf_concentration$AMF <- as.factor(leaf_concentration$AMF)
leaf_concentration$Hg_C <- as.factor(leaf_concentration$Hg_C)
write_csv(leaf_concentration,"Zn_leaf_concentration_summary.csv")
P_leafconcentration <- ggplot(leaf_concentration,aes(x=Hg_C,y=Leaf_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Leaf_Concentration-sd,ymax=Leaf_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Leaf_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Leaf Zn Concentration (mg/kg)")+
  geom_text(x=3.3, y=43, label="AMF**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,45))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_leafconcentration
ggsave("P_leafconcentration_Zn.tiff",width = 5,heigh=4)
##content in the leaf
leveneTest(Leaf_Content ~ AMF*Hg_C,data=element_Zn) #p=0.77
shapiro.test(element_Zn$Leaf_Content) #p=0.08
annova_leaf_content <- lm(Leaf_Content ~ AMF+Hg_C,data=element_Zn)
Anova(annova_leaf_content,type="III") #only AMF (p=0.01164 * )
shapiro.test(residuals(annova_leaf_content))#p=0.05
tukey_leafcontent <- emmeans(annova_leaf_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_leafcontent
leaf_content <- data_summary(element_Zn,varname = "Leaf_Content",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
leaf_content$AMF <- as.factor(leaf_content$AMF)
leaf_content$Hg_C <- as.factor(leaf_content$Hg_C)
write_csv(leaf_content,"Zn_leaf_content_summary.csv")
P_leafcontent <- ggplot(leaf_content,aes(x=Hg_C,y=Leaf_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Leaf_Content-sd,ymax=Leaf_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Leaf_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Leaf Zn Content (µg)")+
  geom_text(x=3.2,y=1400,label="AMF*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,30))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_leafcontent
ggsave("P_leafcontent_Zn.tiff",width = 5,heigh=4)
##concentration in the stem
leveneTest(Stem_Concentration ~ AMF*Hg_C,data=element_Zn) #p=0.53

stem_concentration_25 <- element_Zn%>%
  select(c(1:4,Stem_Concentration))%>%
  filter(Hg_C=="25")
t.test(Stem_Concentration ~ AMF,data = stem_concentration_25,var.equal=T)

annova_stem_concentration <- lm(Stem_Concentration ~ AMF+Hg_C,data=element_Zn)
Anova(annova_stem_concentration,type="III") # AMF (0.000543 ***)
shapiro.test(residuals(annova_stem_concentration))#p=0.2
tukey_stemconcentration <- emmeans(annova_stem_concentration,specs = pairwise ~AMF*Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_stemconcentration
stem_concentration <- data_summary(element_Zn,varname = "Stem_Concentration",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
stem_concentration$AMF <- as.factor(stem_concentration$AMF)
stem_concentration$Hg_C <- as.factor(stem_concentration$Hg_C)
write_csv(stem_concentration,"Zn_stem_concentration_summary.csv")
P_stemconcentration <- ggplot(stem_concentration,aes(x=Hg_C,y=Stem_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Stem_Concentration-sd,ymax=Stem_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Stem_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Stem Zn Concentration (mg/kg)")+
  geom_text(x=3.2,y=80,label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,85))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_stemconcentration
ggsave("P_stemconcentration_Zn.tiff",width = 5,heigh=4)
##content in the stem
leveneTest(Stem_Content ~ AMF*Hg_C,data=element_Zn) #p=0.055
shapiro.test(element_Zn$Stem_Content) #p=0.02
annova_stem_content <- lm(Stem_Content ~ AMF+Hg_C,data=element_Zn)
Anova(annova_stem_content,type="III") #only AMF (0.0008628 ***)
shapiro.test(residuals(annova_stem_content)) #p=0.2798
tukey_stemcontent <- emmeans(annova_stem_content, specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_stemcontent
stem_content <- data_summary(element_Zn,varname = "Stem_Content",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
stem_content$AMF <- as.factor(stem_content$AMF)
stem_content$Hg_C <- as.factor(stem_content$Hg_C)
write_csv(stem_content,"Zn_stem_content_summary.csv")
P_stemcontent <- ggplot(stem_content,aes(x=Hg_C,y=Stem_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Stem_Content-sd,ymax=Stem_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Stem_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Stem Zn Content (µg)")+
  geom_text(x=3.2,y=48,label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,50))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_stemcontent
ggsave("P_stemcontent_Zn.tiff",width = 5,heigh=4)
##concentration in the above including stem and leaf
leveneTest(Above_Concentration ~ AMF*Hg_C,data=element_Zn) #p=0.4


shoot_Zn_25 <- element_Zn %>%
  select(c(1:4,17)) %>%
  filter(Hg_C=="25")
aov_shoot_Zn_25 <- aov(Above_Concentration ~ AMF, data = shoot_Zn_25)
Anova(aov_shoot_Zn_25,type="III")

shoot_Zn_AMF <- element_Zn %>%
  select(c(1:4,17)) %>%
  filter(AMF=="AMF")
aov_shoot_Zn_AMF <- aov(Above_Concentration ~ Hg_C, data = shoot_Zn_AMF)
Anova(aov_shoot_Zn_AMF,type="III")

annova_above_concentration <- aov(Above_Concentration ~ AMF+Hg_C,data=element_Zn)
Anova(annova_above_concentration,type="III") #  AMF (p=0.0004683 ***) 
shapiro.test(residuals(annova_above_concentration)) #P=0.56
tukey_aboveconcentration <- emmeans(annova_above_concentration,specs = pairwise ~AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_aboveconcentration
above_concentration <- data_summary(element_Zn,varname = "Above_Concentration",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","ab","ab","b"))
above_concentration$AMF <- as.factor(above_concentration$AMF)
above_concentration$Hg_C <- as.factor(above_concentration$Hg_C)
write_csv(above_concentration,"Zn_above_concentration_summary.csv")
y_expression <- expression(Shoot~Zn~Concentration~(mg~kg^-1))
P_aboveconcentration <- ggplot(above_concentration,aes(x=Hg_C,y=Above_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Above_Concentration-sd,ymax=Above_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  #geom_text(aes(x=Hg_C,y=Above_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position ="none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y=y_expression)+
  geom_text(x=1,y=58,label="ns")+
  geom_text(x=2,y=58,label="ns")+
  geom_text(x=3,y=58,label="*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,62))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))
P_aboveconcentration
ggsave("P_aboveconcentration_Zn_new_bw.tiff",width = 5,heigh=4)
##content in the above including stem and leaf
leveneTest(Above_Content ~ AMF*Hg_C,data=element_Zn) #p=0.37



annova_above_content <- aov(Above_Content ~ AMF+Hg_C,data=element_Zn)
Anova(annova_above_content,type="III") #AMF (p=3.192e-10 ***)
shapiro.test(residuals(annova_above_content))#P=0.56
tukey_abovecontent <- emmeans(annova_above_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
tukey_abovecontent
above_content <- data_summary(element_Zn,varname = "Above_Content",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","b","b"))
above_content$AMF <- as.factor(above_content$AMF)
above_content$Hg_C <- as.factor(above_content$Hg_C)
write_csv(above_content,"Zn_above_content_summary.csv")
P_abovecontent <- ggplot(above_content,aes(x=Hg_C,y=Above_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Above_Content-sd,ymax=Above_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Above_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Shoot Zn Content (µg)")+
  geom_text(x=3.2,y=72,label="AMF**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,75))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_abovecontent
ggsave("P_abovecontent_Zn.tiff",width = 5,heigh=4)

##Transfer factor
leveneTest(Transfer_Factor ~ AMF*Hg_C,data=element_Zn) #p=0.45

TF_Zn_25 <- element_Zn %>%
  select(c(1:6,19)) %>%
  filter(Hg_C=="25")
aov_TF_Zn_25 <- aov(Transfer_Factor ~ AMF, data = TF_Zn_25)
Anova(aov_TF_Zn_25,type="III")

TF_Zn_Control <- element_Zn %>%
  select(c(1:6,19)) %>%
  filter(AMF=="Control")
aov_TF_Zn_Control <- aov(Transfer_Factor ~ Hg_C, data = TF_Zn_Control)
Anova(aov_TF_Zn_Control,type="III")

annova_TF <- aov(Transfer_Factor ~ AMF*Hg_C,data=element_Zn)
Anova(annova_TF,type="III") #AMF(P=2.482e-05 ***)
shapiro.test(residuals(annova_TF))#P=0.1
tukey.TF <- TukeyHSD(annova_TF)
cld.TF <- multcompLetters4(annova_TF,tukey.TF)
TF <- data_summary(element_Zn,varname = "Transfer_Factor",
                              groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("abc","bc","c","ab","abc","a"))
TF$AMF <- as.factor(TF$AMF)
TF$Hg_C <- as.factor(TF$Hg_C)
write_csv(TF,"Zn_TF_summary.csv")

labeldat <- data.frame(plant_Zn=c("Shoot","Root"),x=c(1.0,1.0),y=c(250,250),
                       label=c("*P*(AMF)<0.001<br>*P*(Hg)=0.17<br>*P*(AMFXHg)=0.46",
                               "*P*(AMF)<0.01<br>*P*(Hg)=0.57<br>*P*(AMFXHg)=0.98"))

P_TF <- ggplot(TF,aes(x=Hg_C,y=Transfer_Factor,group=AMF,fill=AMF))+
  geom_errorbar(aes(ymin=Transfer_Factor-sd,ymax=Transfer_Factor+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Transfer_Factor+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  geom_richtext(x=2,y=0.45,label="*P*(AMF)=0.02<br>*P*(Hg)=0.83<br>*P*(AMFXHg)=0.70", label.color=NA,fill=NA,hjust=0,size=3)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y="Zn transfer factor")+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.5))+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))
P_TF
ggsave("P_TF_Zn_new_bw.tiff",width = 5,heigh=4)

##Zn concentration in the above and root
Zn_above_root_concentration1 <- element_Zn %>%
  select(c(AMF,Hg_C,treatment,Above_Concentration,Root_Concentration))
  

aov_above <- aov(Above_Concentration~AMF*Hg_C,data = Zn_above_root_concentration1)
summary(aov_above)
tukey_above <- TukeyHSD(aov_above)
cld_above <-multcompLetters4(aov_above,tukey_above)
cld_above

aov_root <- aov(Root_Concentration~AMF*Hg_C,data = Zn_above_root_concentration1)
summary(aov_root)
tukey_root <- TukeyHSD(aov_root)
cld_root <-multcompLetters4(aov_root,tukey_root)
cld_root

Zn_above_root_concentration2<- Zn_above_root_concentration1%>%
  pivot_longer(cols = c(Above_Concentration,Root_Concentration),names_to="plant_Zn",
               values_to="Zn_concentration")%>%
  group_by(AMF,Hg_C,treatment,plant_Zn)%>%
  summarise(Zn_mean=mean(Zn_concentration),Zn_sd=sd(Zn_concentration))%>%ungroup()
Zn_above_root_concentration2$treatment <- factor(Zn_above_root_concentration2$treatment,
                                                levels = c("CK0","CK25","CK50","AMF0","AMF25","AMF50")) 
Zn_above_root_concentration2$plant_Zn <- factor(Zn_above_root_concentration2$plant_Zn,
                                               levels = c("Above_Concentration","Root_Concentration"),
                                               labels = c("Shoot","Root"))
Zn_above_root_concentration2<-Zn_above_root_concentration2%>%
  mutate(label=c("ab","a",
                 "b","a",
                 "b","a",
                 "ab","a",
                 "ab","a",
                 "a","a"))

labeldat <- data.frame(plant_Zn=c("Shoot","Root"),x=c(1.0,1.0),y=c(250,250),
                       label=c("*P*(AMF)<0.001<br>*P*(Hg)=0.17<br>*P*(AMFXHg)=0.46",
                               "*P*(AMF)<0.01<br>*P*(Hg)=0.57<br>*P*(AMFXHg)=0.98"))
dat_text <- data.frame(plant_Zn=c("Shoot","Root"),x=c(3.4,3.4),y=c(290,290),label=c("(A)","(B)"))

ggplot(Zn_above_root_concentration2,aes(x=Hg_C,y=Zn_mean,fill=AMF))+
  geom_errorbar(aes(x=Hg_C,ymin=Zn_mean-Zn_sd,ymax=Zn_mean+Zn_sd),position = position_dodge(0.7),width=0.2)+
  geom_bar(stat = "identity",position = position_dodge(0.7),color="black",width = 0.6)+
  geom_text(aes(x=Hg_C,y=Zn_mean+Zn_sd,label=label),vjust=-0.3,position = position_dodge(0.7))+
  facet_grid(~factor(plant_Zn,levels=c("Shoot","Root")))+
  theme_bw()+
  labs(y=Zn ~ Concentration ~ (µg ~ g ^ -1))+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"),
        panel.spacing.x = unit(0,"line"),panel.background = element_blank(),panel.grid = element_blank(),
        strip.text.x = element_text(face = "bold"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,300))+
  scale_fill_manual(values = c("white","grey80"))+
  scale_y_break(c(60,95),scales = 1,expand = F)+
  geom_text(data=dat_text, aes(x=x,y=y,label=label,fill=NULL))+
  geom_richtext(data=labeldat, aes(x=x,y=y,label=label,fill=NULL),
                fill=NA,label.color=NA,hjust=0,vjust=0,size=3,family = "Arial",fontface="plain")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))
ggsave("P_Zn_concentration_above_root_250322.tiff",width = 6,height = 6)

#here i plot leaf,stem and root together

Zn_concentration_1 <- element_Zn%>%
  select(c(treatment,AMF,Hg_C,Leaf_Concentration,Root_Concentration,Stem_Concentration))
Zn_concentration_2 <- Zn_concentration_1%>%
  pivot_longer(cols=c("Leaf_Concentration","Root_Concentration","Stem_Concentration"),names_to="plant_part",
               values_to="Zn_concentration")
Zn_concentration_2$treatment=factor(Zn_concentration_2$treatment,
                                      levels = c("CK0","AMF0","CK25","AMF25","CK50","AMF50"))
Zn_concentration_2$plant_part=factor(Zn_concentration_2$plant_part,
                                     levels =c("Leaf_Concentration","Stem_Concentration","Root_Concentration") )
Zn_concentration_2 <- Zn_concentration_2%>%
  group_by(treatment,AMF,Hg_C,plant_part)%>%
  summarise(Zn_mean=mean(Zn_concentration),Zn_sd=sd(Zn_concentration))%>%ungroup()

Zn_concentration_2<- Zn_concentration_2%>%
  mutate(label=c("a","a","a",
                 "a","b","a",
                 "a","a","b",
                 "a","a","a",
                 "a","a","a",
                 "b","b","a"))
p_root_stem_leaf<-ggplot(Zn_concentration_2,aes(x=treatment,y=Zn_mean,fill=plant_part))+
  geom_errorbar(aes(ymin=Zn_mean-Zn_sd,ymax=Zn_mean+Zn_sd),width=0.2,
                position = position_dodge(0.9))+
  geom_text(aes(x=treatment,y=Zn_mean+Zn_sd,label=label),position = position_dodge(0.9),vjust=-0.3)+
  geom_bar(stat = "identity",position = position_dodge(0.9),color="black")+
  theme_few()+
  scale_x_discrete(label=c("Hg 0","Hg 0*","Hg 25","Hg 25*","Hg 50","Hg 50*"))+
  scale_y_break(c(80,118),scales = 1,expand = F)+
  geom_vline(xintercept = 2.5, color = "blue",size=0.2)+
  geom_vline(xintercept = 4.5, color = "blue",size=0.2)+
  ylab("Zn Concentration (µg/g)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,260))+
  scale_fill_manual(labels=c("leaf","stem","root"),values =c("forestgreen","gold4","gold3"))+
  theme(legend.title = element_blank(),legend.position = "top",
        axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"))
p_root_stem_leaf
ggsave("root_stem_leaf_Zn_concentration.tiff",width = 4,height = 5)










