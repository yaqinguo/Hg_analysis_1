library(tidyverse)
library(RColorBrewer)
library(ggbreak)
library(multcompView)
library(ggthemes)
library(car)
library(emmeans)
library(gam)
library(mgcv)
Whole_Hg<- read_csv("Whole_Hg.csv",col_types = "ffcffffddddddddddddddd")
table(Whole_Hg$AMF,Whole_Hg$Hg_C) #to check data is balanced or not
Whole_Hg$AMF <- as.numeric(Whole_Hg$AMF)
Whole_Hg$Hg_C<- as.numeric(Whole_Hg$Hg_C)
Whole_Hg$Column <- as.numeric(Whole_Hg$Column)
Whole_Hg$Row <- as.numeric(Whole_Hg$Row)
mod1<- gam(Sand ~ Row+Column+Hg_C+AMF,data=Whole_Hg)
summary(mod1)
mod2<- gam(Sand ~ Row+Column+s(Hg_C,1)+s(AMF,1),data=Whole_Hg)
summary(mod2)
pisa = read.csv('data/pisasci2006.csv')

Whole_Hg %>% distinct(AMF)
Whole_Hg %>% distinct(Hg_C)
Whole_Hg %>% distinct(Column)
Whole_Hg %>% distinct(Row)
Whole_Hg %>% distinct(rep)
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
##Hg concentration in substrate
leveneTest(Sand ~ AMF*Hg_C,data=Whole_Hg) #p=0.05
shapiro.test(Whole_Hg$Sand) #p=0.03
annova_sand_concentration <- lm(Sand ~ AMF*Hg_C,data=Whole_Hg)
Anova(annova_sand_concentration,type="III")
annova_sand_concentration <- lm(Sand ~ rep+Column + Row +AMF+Hg_C,data=Whole_Hg)
summary(annova_sand_concentration) #only Hg (p=2.949e-09 ***)
shapiro.test(residuals(annova_sand_concentration)) #p=0.00129
# par(mfrow=c(2,2))
# plot(annova_sand_concentration)
tukey.cld_sandconcentration <- emmeans(annova_sand_concentration,specs = pairwise ~ AMF+Hg_C) %>% 
  multcomp::cld(Letters=letters)
sand_concentration <- data_summary(Whole_Hg,varname = "Sand",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","b","c","a","b","c"))
sand_concentration$AMF <- as.factor(sand_concentration$AMF)
sand_concentration$Hg_C <- as.factor(sand_concentration$Hg_C)
write_csv(sand_concentration,"Hg_sand_concentration_summary.csv")
P_sandconcentration <- ggplot(sand_concentration,aes(x=Hg_C,y=Sand,fill=AMF))+
  geom_errorbar(aes(ymin=Sand-sd,ymax=Sand+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Sand+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Substrate Hg Concentration (mg/kg)")+
  geom_text(x=3.2, y=0.38, label="Hg***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.4))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_sandconcentration 
ggsave("P_sandconcentration_Hg.tiff",width = 5,heigh=4)

##Hg content in substrate
leveneTest(Sand_Content ~ AMF*Hg_C,data=Whole_Hg) #p=0.04
shapiro.test(Whole_Hg$Sand_Content) #p=0.03 
annova_sand_content <- lm(Sand_Content ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_sand_content,type="III") #only Hg (p=3.868e-09 ***)
shapiro.test(residuals(annova_sand_content))#p=0.0004531

annova_sand_content <- lm(Sand_Content ~ AMF+Hg_C,data=Whole_Hg)
mixedmodel <- lmer(Sand_Content~ AMF*Hg_C + AMF + Hg_C + Block + (1|Block:Column) + (1|Block:Row), data=Whole_Hg)
summary(mixedmodel)
shapiro.test(residuals(mixedmodel))
plot(residuals(mixedmodel))
tukey.cld_sandcontent <- emmeans(annova_sand_content,specs = pairwise ~ AMF+Hg_C) %>% 
  multcomp::cld(Letters=letters)
sand_content <- data_summary(Whole_Hg,varname = "Sand_Content",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","b","c","a","b","c"))
sand_content$AMF <- as.factor(sand_content$AMF)
sand_content$Hg_C <- as.factor(sand_content$Hg_C)
write_csv(sand_content,"Hg_sand_content_summary.csv")
P_sandcontent <- ggplot(sand_content,aes(x=Hg_C,y=Sand_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Sand_Content-sd,ymax=Sand_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Sand_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Sand Hg Content (µg)")+
  geom_text(x=3.3, y=19.0, label="Hg***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,20))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_sandcontent 
ggsave("P_sandcontent_Hg.tiff",width = 5,heigh=4)
##Hg concentration in the root
leveneTest(Root_Concentration ~ AMF*Hg_C,data=Whole_Hg) #p=0.4
shapiro.test(Whole_Hg$Root_Concentration) #p=0.0045 
annova_root_concentration <- lm(Root_Concentration ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_root_concentration,type="III") #only Hg (p=2.536e-11 ***)
shapiro.test(residuals(annova_root_concentration))#P=0.5
tukey_rootconcentration <- emmeans(annova_root_concentration,specs = pairwise ~ AMF+Hg_C) %>% 
  multcomp::cld(Letters=letters)
# ##because data is normal, then the following the method is used
# root_concentration <- glm(Root_Concentration ~ AMF+Hg_C,Whole_Hg,family = "quasipoisson")
# summary(root_concentration)
root_concentration <- data_summary(Whole_Hg,varname = "Root_Concentration",
                           groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("ab","c","d","a","bc","d"))
root_concentration$AMF <- as.factor(root_concentration$AMF)
root_concentration$Hg_C <- as.factor(root_concentration$Hg_C)
write_csv(root_concentration,"Hg_root_concentration_summary.csv")
P_rootconcentration <- ggplot(root_concentration,aes(x=Hg_C,y=Root_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Root_Concentration-sd,ymax=Root_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Root_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Root Hg Concentration (mg/kg)")+
  geom_text(x=3.3, y=37.5, label="Hg***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,40))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_rootconcentration
ggsave("P_rootconcentration_Hg.tiff",width = 5,heigh=4)
##root content
leveneTest(Root_Content ~ AMF*Hg_C,data=Whole_Hg) #p=0.4
shapiro.test(Whole_Hg$Root_Content) #p=0.0047 
annova_root_content <- lm(Root_Content ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_root_content,type="III") #only Hg (p=3.91e-11 ***)
shapiro.test(residuals(annova_root_content))#P=0.8
tukey_rootcontent <- emmeans(annova_root_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
root_content <- data_summary(Whole_Hg,varname = "Root_Content",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("ab", "c","d","a","bc","d"))
root_content$AMF <- as.factor(root_content$AMF)
root_content$Hg_C <- as.factor(root_content$Hg_C)
write_csv(root_content,"Hg_root_content_summary.csv")
P_rootcontent <- ggplot(root_content,aes(x=Hg_C,y=Root_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Root_Content-sd,ymax=Root_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Root_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Root Hg Content (µg)")+
  geom_text(x=3.3, y=28.5, label="Hg***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,30))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_rootcontent
ggsave("P_rootcontent_Hg.tiff",width = 5,heigh=4)
##Hg concentration in the leaf
leveneTest(Leaf_Concentration ~ AMF*Hg_C,data=Whole_Hg) #p=0.8
shapiro.test(Whole_Hg$Leaf_Concentration) #p=0.13 
annova_leaf_concentration <- lm(Leaf_Concentration ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_leaf_concentration,type="III") #only AMF (p=2.929e-07)
shapiro.test(residuals(annova_leaf_concentration))#P=0.01
tukey_leafconcentration <- emmeans(annova_leaf_concentration,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
leaf_concentration <- data_summary(Whole_Hg,varname = "Leaf_Concentration",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("b","b","b","a","a","a"))
leaf_concentration$AMF <- as.factor(leaf_concentration$AMF)
leaf_concentration$Hg_C <- as.factor(leaf_concentration$Hg_C)
write_csv(leaf_concentration,"Hg_leaf_concentration_summary.csv")
P_leafconcentration <- ggplot(leaf_concentration,aes(x=Hg_C,y=Leaf_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Leaf_Concentration-sd,ymax=Leaf_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Leaf_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Leaf Hg Concentration (mg/kg)")+
  geom_text(x=3.3, y=2.8, label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_leafconcentration
ggsave("P_leafconcentration_Hg.tiff",width = 5,heigh=4)
##Hg content in the leaf
leveneTest(Leaf_Content ~ AMF*Hg_C,data=Whole_Hg) #p=0.8
shapiro.test(Whole_Hg$Leaf_Content) #p=0.10 
annova_leaf_content <- lm(Leaf_Content ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_leaf_content,type="III") #only AMF (p=1.611e-06 ***)
shapiro.test(residuals(annova_leaf_content))#P=0.0037
tukey_leafcontent <- emmeans(annova_leaf_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
leaf_content <- data_summary(Whole_Hg,varname = "Leaf_Content",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("b","b","b","a","a","a"))
leaf_content$AMF <- as.factor(leaf_content$AMF)
leaf_content$Hg_C <- as.factor(leaf_content$Hg_C)
write_csv(leaf_content,"Hg_leaf_content_summary.csv")
P_leafcontent <- ggplot(leaf_content,aes(x=Hg_C,y=Leaf_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Leaf_Content-sd,ymax=Leaf_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Leaf_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Leaf Hg Content (µg)")+
  geom_text(x=3.2,y=2.8,label="AMF***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_leafcontent
ggsave("P_leafcontent_Hg.tiff",width = 5,heigh=4)
##Hg concentration in the stem
leveneTest(Stem_Concentration ~ AMF*Hg_C,data=Whole_Hg) #p=0.3
shapiro.test(Whole_Hg$Stem_Concentration) #p=5.391e-05 
annova_stem_concentration <- lm(Stem_Concentration ~ AMF*Hg_C,data=Whole_Hg)
Anova(annova_stem_concentration,type="III") # AMF*Hg (0.005789 **)
shapiro.test(residuals(annova_stem_concentration))#P=0.0005
tukey_stemconcentration <- emmeans(annova_stem_concentration,specs = pairwise ~AMF*Hg_C) %>%
  multcomp::cld(Letters=letters)
stem_concentration <- data_summary(Whole_Hg,varname = "Stem_Concentration",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","a","a","b"))
stem_concentration$AMF <- as.factor(stem_concentration$AMF)
stem_concentration$Hg_C <- as.factor(stem_concentration$Hg_C)
write_csv(stem_concentration,"Hg_stem_concentration_summary.csv")
P_stemconcentration <- ggplot(stem_concentration,aes(x=Hg_C,y=Stem_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Stem_Concentration-sd,ymax=Stem_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Stem_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Stem Hg Concentration (mg/kg)")+
  geom_text(x=3.2,y=1.88,label="AMF X Hg**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_stemconcentration
ggsave("P_stemconcentration_Hg.tiff",width = 5,heigh=4)
##Hg content in the stem
leveneTest(Stem_Content ~ AMF*Hg_C,data=Whole_Hg) #p=0.3
shapiro.test(Whole_Hg$Stem_Content) #p=9.567e-05 
annova_stem_content <- lm(Stem_Content ~ AMF*Hg_C,data=Whole_Hg)
Anova(annova_stem_content,type="III") #only AMF*Hg (0.006259 **)
shapiro.test(residuals(annova_stem_content))#P=0.00088
tukey_stemcontent <- emmeans(annova_stem_content, specs = pairwise ~ AMF*Hg_C) %>%
  multcomp::cld(Letters=letters)
stem_content <- data_summary(Whole_Hg,varname = "Stem_Content",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","a","a","b"))
stem_content$AMF <- as.factor(stem_content$AMF)
stem_content$Hg_C <- as.factor(stem_content$Hg_C)
write_csv(stem_content,"Hg_stem_content_summary.csv")
P_stemcontent <- ggplot(stem_content,aes(x=Hg_C,y=Stem_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Stem_Content-sd,ymax=Stem_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Stem_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Stem Hg Content (µg)")+
  geom_text(x=3.2,y=0.95,label="AMF X Hg**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_stemcontent
ggsave("P_stemcontent_Hg.tiff",width = 5,heigh=4)
##Hg concentration in the above including stem and leaf
leveneTest(Above_Concentration ~ AMF*Hg_C,data=Whole_Hg) #p=0.995
shapiro.test(Whole_Hg$Above_Concentration) #p=0.2051 
annova_above_concentration <- lm(Above_Concentration ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_above_concentration,type="III") # both significant AMF (p=9.770e-06 ***) and Hg (p=0.008409 **)
shapiro.test(residuals(annova_above_concentration))#P=0.5
tukey_aboveconcentration <- emmeans(annova_above_concentration,specs = pairwise ~AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
above_concentration <- data_summary(Whole_Hg,varname = "Above_Concentration",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("bc","cd","d","a","ab","bc"))
above_concentration$AMF <- as.factor(above_concentration$AMF)
above_concentration$Hg_C <- as.factor(above_concentration$Hg_C)
write_csv(above_concentration,"Hg_above_concentration_summary.csv")
P_aboveconcentration <- ggplot(above_concentration,aes(x=Hg_C,y=Above_Concentration,fill=AMF))+
  geom_errorbar(aes(ymin=Above_Concentration-sd,ymax=Above_Concentration+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Above_Concentration+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Above Hg Concentration (mg/kg)")+
  geom_text(x=3.2,y=1.88,label="AMF***, Hg**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_aboveconcentration
ggsave("P_aboveconcentration_Hg.tiff",width = 5,heigh=4)
##Hg content in the above including stem and leaf
leveneTest(Above_Content ~ AMF*Hg_C,data=Whole_Hg) #p=0.94
shapiro.test(Whole_Hg$Above_Content) #p=0.3213 
annova_above_content <- lm(Above_Content ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_above_content,type="III") #both significant AMF (p=1.867e-05 ***) and Hg (0.03367 *)
shapiro.test(residuals(annova_above_content))#P=0.36
tukey_abovecontent <- emmeans(annova_above_content,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
above_content <- data_summary(Whole_Hg,varname = "Above_Content",
                             groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("bcd","cd","d","a","ab","abc"))
above_content$AMF <- as.factor(above_content$AMF)
above_content$Hg_C <- as.factor(above_content$Hg_C)
write_csv(above_content,"Hg_above_content_summary.csv")
P_abovecontent <- ggplot(above_content,aes(x=Hg_C,y=Above_Content,fill=AMF))+
  geom_errorbar(aes(ymin=Above_Content-sd,ymax=Above_Content+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Above_Content+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Shoot Hg Content (µg)")+
  geom_text(x=3.2,y=2.8,label="AMF***, Hg*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_abovecontent
ggsave("P_abovecontent_Hg.tiff",width = 5,heigh=4)
##bioaccumulation
leveneTest(Bioaccumulation ~ AMF*Hg_C,data=Whole_Hg) #p=0.39
shapiro.test(Whole_Hg$Bioaccumulation) #p=0.008298 
annova_Bioaccumulation <- lm(Bioaccumulation ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_Bioaccumulation,type="III") #both AMF (p=0.04371 *) and Hg (p=2.408e-11 ***)
shapiro.test(residuals(annova_Bioaccumulation))#P=0.34
tukey_Bioaccumulation <- emmeans(annova_Bioaccumulation,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
Bioaccumulation <- data_summary(Whole_Hg,varname = "Bioaccumulation",
                              groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("ab","c","d","a","bc","d"))
Bioaccumulation$AMF <- as.factor(Bioaccumulation$AMF)
Bioaccumulation$Hg_C <- as.factor(Bioaccumulation$Hg_C)
write_csv(Bioaccumulation,"Hg_Bioaccumulation_summary.csv")
P_Bioaccumulation <- ggplot(Bioaccumulation,aes(x=Hg_C,y=Bioaccumulation,fill=AMF))+
  geom_errorbar(aes(ymin=Bioaccumulation-sd,ymax=Bioaccumulation+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Bioaccumulation+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Bioaccumulation of Hg")+
  geom_text(x=3.2,y=14,label="AMF*, Hg***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,15))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_Bioaccumulation 
ggsave("P_Bioaccumulation _Hg.tiff",width = 5,heigh=4)
##transfer factor
leveneTest(Transfer_Factor ~ AMF*Hg_C,data=Whole_Hg) #p=4.731e-10 ***
shapiro.test(Whole_Hg$Transfer_Factor) #p=4.641e-07
annova_TF <- lm(Transfer_Factor ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_TF,type="III") #Hg p=3.007e-05 ***
shapiro.test(residuals(annova_TF))#P=0.0005951
tukey_TF <- emmeans(annova_TF,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
TF <- data_summary(Whole_Hg,varname = "Transfer_Factor",
                                groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("b","a","a","b","a","a"))
TF$AMF <- as.factor(TF$AMF)
TF$Hg_C <- as.factor(TF$Hg_C)
write_csv(TF,"Hg_TF_summary.csv")
P_TF<- ggplot(TF,aes(x=Hg_C,y=Transfer_Factor,fill=AMF))+
  geom_errorbar(aes(ymin=Transfer_Factor-sd,ymax=Transfer_Factor+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg_C,y=Transfer_Factor+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Transfer Factor of Hg")+
  geom_text(x=3.2,y=5.6,label="Hg***")+
  scale_y_continuous(expand = c(0,0),limits = c(0,6))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_TF 
ggsave("P_TF_Hg.tiff",width = 5,heigh=4)
library(ggpubr)
ggarrange(P_leafconcentration,P_stemconcentration,P_rootconcentration,ncol = 1,common.legend = T)
