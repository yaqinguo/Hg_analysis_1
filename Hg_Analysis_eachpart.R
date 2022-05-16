library(tidyverse)
library(RColorBrewer)
library(ggbreak)
library(multcompView)
library(ggthemes)
library(car)
library(emmeans)
library(ggpubr)
library(ggtext)
library(patchwork)
Whole_Hg <- read_csv("Whole_Hg.csv",col_types = "ffcfffffddddddddddddddd")
table(Whole_Hg$AMF,Whole_Hg$Hg_C) #to check data is balanced or not
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
# Three_part_Hg_1 <- Whole_Hg %>%
#   select(c("treatment","Leaf_Concentration","Stem_Concentration","Root_Concentration"))
# Three_part_Hg_2 <- Three_part_Hg_1 %>%
#   pivot_longer(cols = c("Leaf_Concentration","Stem_Concentration","Root_Concentration"),
#                names_to="plant_Hg",values_to="Hg_concentration")
# Three_part_Hg_2$plant_Hg = factor(Three_part_Hg_2$plant_Hg,
#                                   levels = c("Leaf_Concentration","Stem_Concentration","Root_Concentration"))
# 
# ggbarplot(Three_part_Hg_2,x="treatment",y="Hg_concentration",add = "mean_sd",fill="plant_Hg",
#           palette = c("forestgreen","gold4","gold3"),width = 0.5)
# ggplot(Three_part_Hg_2,aes(x=plant_Hg,y=Hg_concentration))+
#   facet_wrap(~ treatment)+
#   geom_bar(position = position_dodge(),stat = "identity")+
#   theme_few()+
#   theme(aspect.ratio = 1,axis.title.x = element_blank(),legend.position = "none",
#         axis.ticks.length.y = unit(-0.15,"cm"),axis.ticks.length.x = unit(0,"cm"),axis.line = )+
#   scale_x_discrete(labels=c("Hg 0","Hg 0*","Hg 25","Hg 25*","Hg 50","Hg 50*"))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,40))+
#   labs(y="Hg_Concentration (mg/kg)")
# ggsave("leaf_stem_root.tiff",width = 5,height = 4)

#pie chart
pie_data <- Whole_Hg %>%
  select(c(AMF,Hg_C,Root_Content,Stem_Content,Leaf_Content))%>%
  mutate(Totall_content=Root_Content+Stem_Content+Leaf_Content)
leveneTest(Totall_content ~ AMF*Hg_C, data = pie_data)
aov_total_content <- aov(Totall_content ~ AMF+Hg_C, data = pie_data)
Anova(aov_total_content,type="III")
shapiro.test(residuals(aov_total_content))
TukeyHSD(aov_total_content)
Total_0 <- pie_data %>%
  filter(Hg_C=="0")
t.test(Totall_content ~ AMF, data=Total_0,var.equal=T)
Total_50 <- pie_data %>%
  filter(Hg_C=="50")
t.test(Totall_content ~ AMF, data=Total_50,var.equal=T)
wilcox.test(Totall_content ~ AMF, data=Total_50)



# AMF_0
Pie_AMF_0 <- pie_data %>%
  filter(AMF=="AMF") %>%
  filter(Hg_C=="0") %>%
  pivot_longer(cols = c("Root_Content","Stem_Content","Leaf_Content"),names_to="Plant_part",
               values_to="Content")%>%
  mutate(percent=Content/sum(Content)*100)
summary_data_AMF_0 <- Pie_AMF_0 %>%
  group_by(Plant_part)%>%
  summarise(total=sum(percent))%>%
  arrange(desc(Plant_part))
p_AMF_0 <- ggplot(summary_data_AMF_0,aes(x="",y=total,fill=Plant_part))+
  geom_bar(width=1,stat = "identity",color="white")+
  coord_polar("y",start = 0)+
  theme_void()+
  theme(legend.position = "none",aspect.ratio = 1)
p_AMF_0

# AMF_25
Pie_AMF_25 <- pie_data %>%
  filter(AMF=="AMF") %>%
  filter(Hg_C=="25") %>%
  pivot_longer(cols = c("Root_Content","Stem_Content","Leaf_Content"),names_to="Plant_part",
               values_to="Content")%>%
  mutate(percent=Content/sum(Content)*100)
summary_data_AMF_25 <- Pie_AMF_25 %>%
  group_by(Plant_part)%>%
  summarise(total=sum(percent))%>%
  arrange(desc(Plant_part))
p_AMF_25 <- ggplot(summary_data_AMF_25,aes(x="",y=total,fill=Plant_part))+
  geom_bar(width=1,stat = "identity",color="white")+
  coord_polar("y",start = 0)+
  theme_void()+
  theme(legend.position = "none",aspect.ratio = 1)
p_AMF_25

# AMF_50
Pie_AMF_50 <- pie_data %>%
  filter(AMF=="AMF") %>%
  filter(Hg_C=="50") %>%
  pivot_longer(cols = c("Root_Content","Stem_Content","Leaf_Content"),names_to="Plant_part",
               values_to="Content")%>%
  mutate(percent=Content/sum(Content)*100)
summary_data_AMF_50 <- Pie_AMF_50 %>%
  group_by(Plant_part)%>%
  summarise(total=sum(percent))%>%
  arrange(desc(Plant_part))
p_AMF_50 <- ggplot(summary_data_AMF_50,aes(x="",y=total,fill=Plant_part))+
  geom_bar(width=1,stat = "identity",color="white")+
  coord_polar("y",start = 0)+
  theme_void()+
  theme(legend.position = "none",aspect.ratio = 1)
p_AMF_50

# Control_0
Pie_Control_0 <- pie_data %>%
  filter(AMF=="Control") %>%
  filter(Hg_C=="0") %>%
  pivot_longer(cols = c("Root_Content","Stem_Content","Leaf_Content"),names_to="Plant_part",
               values_to="Content")%>%
  mutate(percent=Content/sum(Content)*100)
summary_data_Control_0 <- Pie_Control_0 %>%
  group_by(Plant_part)%>%
  summarise(total=sum(percent))%>%
  arrange(desc(Plant_part))
p_Control_0 <- ggplot(summary_data_Control_0,aes(x="",y=total,fill=Plant_part))+
  geom_bar(width=1,stat = "identity",color="white")+
  coord_polar("y",start = 0)+
  theme_void()+
  theme(legend.position = "none",aspect.ratio = 1)
p_Control_0

# Control_25
Pie_Control_25 <- pie_data %>%
  filter(AMF=="Control") %>%
  filter(Hg_C=="25") %>%
  pivot_longer(cols = c("Root_Content","Stem_Content","Leaf_Content"),names_to="Plant_part",
               values_to="Content")%>%
  mutate(percent=Content/sum(Content)*100)
summary_data_Control_25 <- Pie_Control_25 %>%
  group_by(Plant_part)%>%
  summarise(total=sum(percent))%>%
  arrange(desc(Plant_part))
p_Control_25 <- ggplot(summary_data_Control_25,aes(x="",y=total,fill=Plant_part))+
  geom_bar(width=1,stat = "identity",color="white")+
  coord_polar("y",start = 0)+
  theme_void()+
  theme(aspect.ratio = 1,legend.position = "none")
p_Control_25
# Control_50
Pie_Control_50 <- pie_data %>%
  filter(AMF=="Control") %>%
  filter(Hg_C=="50") %>%
  pivot_longer(cols = c("Root_Content","Stem_Content","Leaf_Content"),names_to="Plant_part",
               values_to="Content")%>%
  mutate(percent=Content/sum(Content)*100)
summary_data_Control_50 <- Pie_Control_50 %>%
  group_by(Plant_part)%>%
  summarise(total=sum(percent))%>%
  arrange(desc(Plant_part))
p_Control_50 <- ggplot(summary_data_Control_50,aes(x="",y=total,fill=Plant_part))+
  geom_bar(width=1,stat = "identity",color="white")+
  coord_polar("y",start = 0)+
  theme_void()+
  theme(aspect.ratio = 1,legend.position = "none")
p_Control_50 
p_AMF_0+p_AMF_25+p_AMF_50+p_Control_0+p_Control_25+p_Control_50
ggsave("pie_chart_Hg_distrubtion.tiff",width = 5,height = 4)
#here i only care about shoot and root
pie_data <- Whole_Hg %>%
  select(c(AMF,Hg_C,Root_Content,Above_Content))%>%
  pivot_longer(cols = c(Root_Content,Above_Content),names_to="Plant_part",values_to="Content")
pie_data1<- pie_data %>%
  group_by(AMF,Hg_C)%>%
  summarise(percent=(Content/sum(Content))*100,across())%>% ungroup()
pie_0 <- pie_data1%>%
  filter(Plant_part=="Above_Content")%>%
  filter(Hg_C=="0")
leveneTest(percent ~ AMF,data=pie_0)
aov_percent <- aov(percent ~ AMF,data=pie_0)
Anova(aov_percent,type="III")



pie_data2<- pie_data1%>%
  group_by(AMF,Hg_C,Plant_part)%>%
  summarise(total=sum(percent))%>%ungroup()


#phytoextraction efficiency
colnames(Whole_Hg)
phytoextraction_efficiency <- Whole_Hg %>%
  select(c(AMF,Hg_C,Above_Content,Root_dry)) %>%
  mutate(phytoextraction_efficiency=Above_Content/Root_dry)
leveneTest(phytoextraction_efficiency ~ AMF*Hg_C,data = phytoextraction_efficiency)
aov_phytoextraction_efficiency <- aov(phytoextraction_efficiency ~ AMF+Hg_C,data = phytoextraction_efficiency)
Anova(aov_phytoextraction_efficiency,type="III")
shapiro.test(residuals(aov_phytoextraction_efficiency))
tukey_phytoextraction_efficiency <- TukeyHSD(aov_phytoextraction_efficiency)
cld_phytoextraction_efficiency <- multcompLetters4(aov_phytoextraction_efficiency, tukey_phytoextraction_efficiency)
cld_phytoextraction_efficiency 
phytoextraction_efficiency_25 <- phytoextraction_efficiency %>%
  filter(Hg_C=="25")
aov_phytoextraction_efficiency_25 <- aov(phytoextraction_efficiency ~ AMF,data = phytoextraction_efficiency_25)
Anova(aov_phytoextraction_efficiency_25,type="III")
phytoextraction_efficiency_data <- data_summary(phytoextraction_efficiency,varname = "phytoextraction_efficiency",
                                                groupnames = c("AMF","Hg_C"))%>%
  mutate(tukey=c("a","a","a","b","b","a"))
y_expression <- expression(Hg~phytoextraction~efficiency~(µg~g^-1))
P_phytoextraction_efficiency <- ggplot(phytoextraction_efficiency_data,aes(x=Hg_C,y=phytoextraction_efficiency,group=AMF,fill=AMF))+
  #geom_errorbar(aes(ymin=phytoextraction_efficiency-sd,ymax=phytoextraction_efficiency+sd),width=0.1)+
  geom_line()+
  geom_point(size=4,shape=21)+
  #geom_text(aes(x=Hg_C,y=phytoextraction_efficiency+sd,label=tukey))+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position = c(0.1,0.93),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y=y_expression)+
  geom_text(x=1, y=2.85, label="***",size=3,family="Arial")+
  geom_text(x=2,y=2.85,label="*",size=3,family="Arial")+
  geom_text(x=3,y=2.9,label="ns",size=3,family="Arial")+
  scale_y_continuous(expand = c(0,0),limits = c(0,3.5))+
  scale_x_discrete(breaks=c("0","25","50"),labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))
P_phytoextraction_efficiency
ggsave("P_phytoextraction_efficiency_linepot.tiff",width = 5,heigh=4)
#uptake efficiency
uptake_efficiency <- Whole_Hg %>%
  select(c(AMF,Hg_C,Above_Content,Root_Content,Root_dry)) %>%
  mutate(uptake_efficiency=(Above_Content+Root_Content)/Root_dry)
leveneTest(uptake_efficiency ~ AMF*Hg_C,data = uptake_efficiency)
aov_uptake_efficiency <- aov(uptake_efficiency ~ AMF+Hg_C,data = uptake_efficiency)
Anova(aov_uptake_efficiency,type="III")
shapiro.test(residuals(aov_uptake_efficiency))
uptake_efficiency_50 <- uptake_efficiency%>%
  filter(Hg_C=="50")
aov_uptake_efficiency_50 <- aov(uptake_efficiency ~ AMF,data = uptake_efficiency_50)
Anova(aov_uptake_efficiency_50,type="III")
uptake_efficiency_data <- data_summary(uptake_efficiency,varname = "uptake_efficiency",
                                       groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","a","a"))
P_uptake_efficiency <- ggplot(uptake_efficiency_data,aes(x=Hg_C,y=uptake_efficiency,group=AMF,color=AMF))+
  #geom_errorbar(aes(ymin=uptake_efficiency-sd,ymax=uptake_efficiency+sd),width=0.1)+
  geom_line()+
  geom_point(size=4,shape=21)+
  #geom_text(aes(x=Hg_C,y=uptake_efficiency+0.7,label=tukey))+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family = "Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Hg_uptake_efficiency")+
  #geom_text(x=3.2, y=38, label="Hg***, AMF*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,40))+
  scale_color_manual(values=c('#999999','#E69F00'))
P_uptake_efficiency
ggsave("P_uptake_efficiency.tiff",width = 5,heigh=4)
#TF
transfer_efficiency <- Whole_Hg %>%
  select(c(AMF,Hg_C,Above_Content,Root_Content)) %>%
  mutate(transfer_efficiency=Above_Content/Root_Content)
leveneTest(transfer_efficiency ~ AMF*Hg_C,data = transfer_efficiency)
aov_transfer_efficiency <- aov(transfer_efficiency ~ AMF+Hg_C,data = transfer_efficiency)
Anova(aov_transfer_efficiency,type="III")
shapiro.test(residuals(aov_transfer_efficiency))
transfer_efficiency_50 <- transfer_efficiency%>%
  filter(Hg_C=="50")
aov_transfer_efficiency_50 <- aov(transfer_efficiency ~ AMF,data = transfer_efficiency_50)
Anova(aov_transfer_efficiency_50,type="III")
transfer_efficiency_data <- data_summary(transfer_efficiency,varname = "transfer_efficiency",
                                       groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","a","a","b","a","a"))
P_transfer_efficiency <- ggplot(transfer_efficiency_data,aes(x=Hg_C,y=transfer_efficiency,group=AMF,fill=AMF))+
  #geom_errorbar(aes(ymin=uptake_efficiency-sd,ymax=uptake_efficiency+sd),width=0.1)+
  geom_line()+
  geom_point(size=4,shape=21)+
  #geom_text(aes(x=Hg_C,y=uptake_efficiency+0.7,label=tukey))+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family = "Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Hg_transfer_efficiency")+
  #geom_text(x=3.2, y=38, label="Hg***, AMF*")+
  scale_y_continuous(expand = c(0,0),limits = c(0,5))+
  scale_fill_manual(values=c('#999999','#E69F00'))
P_transfer_efficiency
ggsave("P_P_transfer_efficiency_lineplot.tiff",width = 5,heigh=4)
##Hg concentration in substrate
leveneTest(Sand ~ AMF*Hg_C,data=Whole_Hg) #p=0.05
annova_sand_concentration <- aov(Sand ~ AMF*Hg_C,data=Whole_Hg)
summary(annova_sand_concentration)
#Anova(annova_sand_concentration,type="III")
# annova_sand_concentration <- lm(Sand ~ rep+Column + Row +AMF+Hg_C,data=Whole_Hg)
# summary(annova_sand_concentration) #only Hg (p=2.949e-09 ***)
shapiro.test(residuals(annova_sand_concentration)) #p=0.00129
# par(mfrow=c(2,2))
# plot(annova_sand_concentration)
##only consider Hg25 and Hg50
aov_sand <- aov(Sand ~ treatment,data=Whole_Hg)
summary(aov_sand)
tukey_sand <- TukeyHSD(aov_sand)
cld.sand <- multcompLetters4(aov_sand,tukey_sand)

sand_50 <- Whole_Hg %>%
  select(c(1:8,18)) %>%
  filter(Hg_C=="50")
# aov_25_50 <- aov(Sand ~ AMF+Hg_C,data = sand_25_50)  
# summary(aov_25_50)
t.test(Sand ~ AMF, data=sand_50, var.equal = TRUE) 
  
sand_concentration_50 <- Whole_Hg %>%
  select(c(1:8,18)) %>%
  filter(Hg_C=="50")
aov_sand_concentration_50 <- aov(Sand ~ AMF, sand_concentration_50)
Anova(aov_sand_concentration_50,type="III")

sand_concentration_Control <- Whole_Hg %>%
  select(c(1:7,17)) %>%
  filter(AMF=="Control")
aov_sand_concentration_Control <- aov(Sand ~ Hg_C, sand_concentration_Control)
Anova(aov_sand_concentration_Control,type="III")
tukey_sand_concentration_Control <- TukeyHSD(aov_sand_concentration_Control)
cld_sandconcentrtion_Control <- multcompLetters4(aov_sand_concentration_Control, tukey_sand_concentration_Control)
cld_sandconcentrtion_Control 
tukey.cld_sandconcentration <- emmeans(annova_sand_concentration,specs = pairwise ~ AMF+Hg_C) %>% 
  multcomp::cld(Letters=letters)
sand_concentration <- data_summary(Whole_Hg,varname = "Sand",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("","b","ab","","b","a"))
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
  theme(aspect.ratio = 1,
        legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"),axis.title.x = element_blank())+
  labs(y=Substrate ~ Hg ~ Concentration (µg ~ g^-1))+
  geom_richtext(x=1.8,y=0.19,label="*P*(AMF)=0.33<br>*P*(Hg)<0.001<br>*P*(AMFXHg)=0.22",fill=NA,label.color=NA,
                hjust=0,vjust=0,size=3,family = "Arial",fontface="plain")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.23))+
  scale_fill_manual(values=c('white','grey80'),labels=c("Control","AMF"))
P_sandconcentration 
ggsave("P_sandconcentration_Hg_bw_220322.tiff",width = 5,heigh=4)

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
root_concentration_25 <- Whole_Hg %>%
  select(c(1:8,15)) %>%
  filter(Hg_C=="25")
t.test(Root_Concentration ~ AMF, data=root_concentration_25, var.equal = TRUE)
# aov_root_concentration_50 <- aov(Root_Concentration ~ AMF, root_concentration_50)
# Anova(aov_root_concentration_50,type="III")

root_concentration_0_25 <- Whole_Hg %>%
  select(c(1:8,15)) %>%
  filter(Hg_C=="0"&AMF=="Control"|Hg_C=="25"&AMF=="Control")
t.test(Root_Concentration ~ Hg_C, root_concentration_0_25, var.equal = TRUE)

root_concentration_0_50 <- Whole_Hg %>%
  select(c(1:8,15)) %>%
  filter(Hg_C=="0"&AMF=="AMF"|Hg_C=="50"&AMF=="AMF")
t.test(Root_Concentration ~ Hg_C, root_concentration_0_50, var.equal = TRUE)


root_concentration_Control <- Whole_Hg %>%
  select(c(1:7,14)) %>%
  filter(AMF=="Control")
aov_root_concentration_Control <- aov(Root_Concentration ~ Hg_C, root_concentration_Control)
Anova(aov_root_concentration_Control,type="III")
tukey_root_concentration_Control <- TukeyHSD(aov_root_concentration_Control)
cld_rootconcentrtion_Control <- multcompLetters4(aov_root_concentration_Control, tukey_root_concentration_Control)
cld_rootconcentrtion_Control 

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
  mutate(tukey=c("c","b","a","c","b","a"))
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
ggsave("P_rootconcentration_Hg_new.tiff",width = 5,heigh=4)
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
  mutate(tukey=c("c","b","a","c","b","a"))
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
ggsave("P_rootcontent_Hg_new.tiff",width = 5,heigh=4)
##Hg concentration in the leaf
leveneTest(Leaf_Concentration ~ AMF*Hg_C,data=Whole_Hg) #p=0.8
shapiro.test(Whole_Hg$Leaf_Concentration) #p=0.13 

leaf_concentration_25 <- Whole_Hg%>%
  select(1:8,Leaf_Concentration)%>%
  filter(Hg_C=="25")
t.test(Leaf_Concentration ~ AMF, data=leaf_concentration_25, var.equal = TRUE)


leaf_concentration_0_50 <- Whole_Hg%>%
  select(1:8,Leaf_Concentration)%>%
  filter(Hg_C=="0"&AMF=="AMF"|Hg_C=="50"&AMF=="AMF")
t.test(Leaf_Concentration ~ Hg_C, data=leaf_concentration_0_50, var.equal = TRUE)

leaf_concentration_0_25 <- Whole_Hg%>%
  select(1:8,Leaf_Concentration)%>%
  filter(Hg_C=="0"&AMF=="Control"|Hg_C=="25"&AMF=="Control")
t.test(Leaf_Concentration ~ Hg_C, data=leaf_concentration_0_25, var.equal = TRUE)


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

stem_concentration_25 <- Whole_Hg%>%
  select(1:8,Stem_Concentration)%>%
  filter(Hg_C=="25")
t.test(Stem_Concentration ~ AMF, data=stem_concentration_25, var.equal = TRUE)

stem_concentration_0_25 <- Whole_Hg%>%
  select(1:8,Stem_Concentration)%>%
  filter(Hg_C=="0"&AMF=="AMF"|Hg_C=="25"&AMF=="AMF")
t.test(Stem_Concentration ~ Hg_C, data=stem_concentration_0_25, var.equal = TRUE)

stem_concentration_0_50 <- Whole_Hg%>%
  select(1:8,Stem_Concentration)%>%
  filter(Hg_C=="0"&AMF=="Control"|Hg_C=="50"&AMF=="Control")
t.test(Stem_Concentration ~ Hg_C, data=stem_concentration_0_50, var.equal = TRUE)








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

shoot_concentration_0 <- Whole_Hg %>%
  select(c(1:7,20)) %>%
  filter(Hg_C=="0")
aov_shoot_concentration_0 <- aov(Above_Concentration ~ AMF, shoot_concentration_0)
Anova(aov_shoot_concentration_0,type="III")

shoot_concentration_Control <- Whole_Hg %>%
  select(c(1:7,20)) %>%
  filter(AMF=="Control")
aov_shoot_concentration_Control <- aov(Above_Concentration ~ Hg_C, shoot_concentration_Control)
Anova(aov_shoot_concentration_Control,type="III")
tukey_shoot_concentration_Control <- TukeyHSD(aov_shoot_concentration_Control)
cld_shootconcentrtion_Control <- multcompLetters4(aov_shoot_concentration_Control, tukey_shoot_concentration_Control)
cld_shootconcentrtion_Control


annova_above_concentration <- lm(Above_Concentration ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_above_concentration,type="III") # both significant AMF (p=9.770e-06 ***) and Hg (p=0.008409 **)
shapiro.test(residuals(annova_above_concentration))#P=0.5
tukey_aboveconcentration <- emmeans(annova_above_concentration,specs = pairwise ~AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
above_concentration <- data_summary(Whole_Hg,varname = "Above_Concentration",
                                   groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("c","c","ac","b","ab","a"))
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
  labs(x="Hg_Concentration (mg/kg)",y="Shoot Hg Concentration (mg/kg)")+
  geom_text(x=3.2,y=1.88,label="AMF***, Hg**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))+
  scale_fill_manual(values=c('#999999','#E69F00'),labels=c("Control","AMF"))
P_aboveconcentration
ggsave("P_aboveconcentration_Hg_new.tiff",width = 5,heigh=4)
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
  mutate(tukey=c("c","c","ac","b","ab","a"))
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
ggsave("P_abovecontent_Hg_new.tiff",width = 5,heigh=4)
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

TF_50 <- Whole_Hg %>%
  select(c(1:7,22)) %>%
  filter(Hg_C=="50")
aov_TF_50 <- aov(Transfer_Factor ~ AMF, TF_50)
Anova(aov_TF_50,type="III")

TF_AMF <- Whole_Hg %>%
  select(c(1:7,22)) %>%
  filter(AMF=="AMF")
aov_TF_AMF <- aov(Transfer_Factor ~ Hg_C, TF_AMF)
Anova(aov_TF_AMF,type="III")
tukey_TF_AMF <- TukeyHSD(aov_TF_AMF)
cld_TF_AMF <- multcompLetters4(aov_TF_AMF, tukey_TF_AMF)
cld_TF_AMF

annova_TF <- lm(Transfer_Factor ~ AMF+Hg_C,data=Whole_Hg)
Anova(annova_TF,type="III") #Hg p=3.007e-05 ***
shapiro.test(residuals(annova_TF))#P=0.0005951
tukey_TF <- emmeans(annova_TF,specs = pairwise ~ AMF+Hg_C) %>%
  multcomp::cld(Letters=letters)
TF <- data_summary(Whole_Hg,varname = "Transfer_Factor",
                                groupnames = c("AMF","Hg_C")) %>%
  mutate(tukey=c("a","b","b","a","b","b"))
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
ggsave("P_TF_Hg_new.tiff",width = 5,heigh=4)
library(ggpubr)
ggarrange(P_leafconcentration,P_stemconcentration,P_rootconcentration,ncol = 1,common.legend = T)

#Hg concentration in the shoot and root
Hg_root_shoot_concentration1 <- Whole_Hg %>%
  select(c(treatment,AMF,Hg_C,Above_Concentration,Root_Concentration)) 
above <- aov(Above_Concentration~AMF*Hg_C,data=Hg_root_shoot_concentration1)
summary(above)
tukey.above <- TukeyHSD(above)
multcompLetters4(above,tukey.above)
root <- aov(Root_Concentration~AMF*Hg_C,data=Hg_root_shoot_concentration1)
summary(root)
aov_above <- aov(Above_Concentration ~ treatment, data = Hg_root_shoot_concentration1)
summary(aov_above)
tukey_above <- TukeyHSD(aov_above)
cld_above <- multcompLetters4(aov_above, tukey_above)
cld_above

aov_root <- aov(Root_Concentration ~ treatment, data = Hg_root_shoot_concentration1)
summary(aov_root)
tukey_root <- TukeyHSD(aov_root)
cld_root <- multcompLetters4(aov_root, tukey_root)
cld_root

Above_Hg <- Hg_root_shoot_concentration1%>%
  select(c(AMF,Hg_C,treatment,Above_Concentration))%>%
  group_by(AMF,Hg_C,treatment)%>%
  summarise(Above_mean=mean(Above_Concentration),Above_sd=sd(Above_Concentration))%>%ungroup()
Above_Hg$treatment=factor(Above_Hg$treatment,levels = c("CK0","AMF0","CK25","AMF25","CK50","AMF50"))


p_above<-ggplot(Above_Hg,aes(x=Hg_C,y=Above_mean,fill=AMF))+
  geom_errorbar(aes(x=Hg_C,ymin=Above_mean-Above_sd,ymax=Above_mean+Above_sd),position = position_dodge(0.6),width=0.2)+
  geom_bar(stat = "identity",position = position_dodge(0.6),color="black",width = 0.5)+
  #geom_text(aes(x=Hg_C,y=Above_mean+Above_sd,label=label),vjust=-0.3,position = position_dodge(0.6))+
  theme_few()+
  ylab("Shoot Hg Concentration (µg/g)")+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),legend.position ="none",
        axis.ticks.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"),plot.margin=margin(0,0,0,0,unit = "pt"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,38))+
  scale_y_break(c(3,5),scales = 1,expand = F)+
  scale_fill_manual(values = c("white","grey80"))
p_above
Root_Hg <- Hg_root_shoot_concentration1%>%
  select(c(AMF,Hg_C,treatment,Root_Concentration))%>%
  group_by(AMF,Hg_C,treatment)%>%
  summarise(Root_mean=mean(Root_Concentration),Root_sd=sd(Root_Concentration))%>%ungroup()
Root_Hg$treatment=factor(Root_Hg$treatment,levels = c("CK0","AMF0","CK25","AMF25","CK50","AMF50"))

p_root <- ggplot(Root_Hg,aes(x=Hg_C,y=Root_mean,fill=AMF))+
  geom_errorbar(aes(x=Hg_C,ymin=Root_mean-Root_sd,ymax=Root_mean+Root_sd),position = position_dodge(0.6),width=0.2)+
  geom_bar(stat = "identity",position = position_dodge(0.6),color="black",width = 0.5)+
  #geom_text(aes(x=treatment,y=Above_mean+Above_sd,label=label),vjust=-0.3,position = position_dodge(0.9))+
  theme_few()+
  scale_y_continuous(expand = c(0,0),limits = c(0,38),position = "right",sec.axis = dup_axis())+
  scale_y_break(c(3,5),scales = 1,expand = F)+
  ylab("Root Hg Concentration (µg/g)")+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"),plot.margin=margin(0,0,0,0,unit = "pt"),
        axis.title.y.left = element_blank(),axis.text.y.left = element_blank(),
        axis.ticks.y.left=element_blank()
        )+
  scale_fill_manual(values = c("white","grey80"))
p_root
p_above|p_root

Hg_root_shoot_concentration2 <- Hg_root_shoot_concentration1 %>%
  pivot_longer(cols = c(Above_Concentration,Root_Concentration),names_to="plant_Hg",
               values_to="Hg_concentration")%>%
  group_by(AMF,Hg_C,treatment,plant_Hg)%>%
  summarise(mean=mean(Hg_concentration),sd=sd(Hg_concentration))%>%ungroup()
Hg_root_shoot_concentration2$treatment=factor(Hg_root_shoot_concentration2$treatment,
                                          levels = c("CK0","AMF0","CK25","AMF25","CK50","AMF50"))
Hg_root_shoot_concentration2$plant_Hg=factor(Hg_root_shoot_concentration2$plant_Hg,
                                             levels = c("Above_Concentration","Root_Concentration"),
                                             labels = c("Shoot","Root"))   
Hg_root_shoot_concentration2 <- Hg_root_shoot_concentration2%>%
  mutate(label=c("ab","c",
                 "a","b",
                 "a","a",
                 "c","c",
                 "bc","bc",
                 "ab","a"))
labeldat <- data.frame(plant_Hg=c("Shoot","Root"),x=c(1.0,1.0),y=c(30,30),
                       label=c("*P*(AMF)<0.001<br>*P*(Hg)<0.01<br>*P*(AMFXHg)=0.28",
                               "*P*(AMF)=0.13<br>*P*(Hg)<0.001<br>*P*(AMFXHg)=0.14"))
dat_text <- data.frame(plant_Hg=c("Shoot","Root"),x=c(3.4,3.4),y=c(37,37),label=c("(A)","(B)"))
P_above_root <- ggplot(Hg_root_shoot_concentration2,aes(x=Hg_C,y=mean,fill=AMF))+
  geom_errorbar(aes(x=Hg_C,ymin=mean-sd,ymax=mean+sd),position = position_dodge(0.7),width=0.2)+
  geom_bar(stat = "identity",position = position_dodge(0.7),color="black",width = 0.6)+
  geom_text(aes(x=Hg_C,y=mean+sd,label=label),vjust=-0.3,position = position_dodge(0.7))+
  theme_bw(base_size = 10,base_family = "Arial")+
  #geom_text(data=Hg_root_shoot_concentration2%>%filter(plant_Hg=="Shoot"),x=1.5,y=25,label="*P*(AMF)<0.01<br>*P*(Hg)=0.88<br>*P*(AMFXHg)=0.54",size=3)+
  #geom_text(data=Hg_root_shoot_concentration2%>%filter(plant_Hg=="Root"),x=1.5,y=25,label="*P*(AMF)<0.01<br>*P*(Hg)=0.88<br>*P*(AMFXHg)=0.54",size=3)+
  facet_grid(~factor(plant_Hg,levels=c("Shoot","Root")))+
  labs(y=Hg ~ Concentration ~ (µg ~ g^-1))+
  theme(axis.title.x = element_blank(),legend.title = element_blank(),legend.position = "none",
        axis.ticks.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"),
        panel.spacing.x = unit(0,"line"),panel.background = element_blank(),panel.grid = element_blank(),
        strip.text.x = element_text(face = "bold"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,38))+
  scale_fill_manual(values = c("white","grey80"))+
  scale_y_break(c(3,7),scales = 1,expand = F)+
  geom_text(data=dat_text, aes(x=x,y=y,label=label,fill=NULL))+
  geom_richtext(data=labeldat, aes(x=x,y=y,label=label,fill=NULL),
                fill=NA,label.color=NA,hjust=0,vjust=0,size=3,family = "Arial",fontface="plain")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))
P_above_root
#gtable_show_names(P_above_root)
#reposition_legend(P_above_root,"top left",panel = "panel-1-2")
# library(lemon)
# gtable_show_names(P_above_root)
ggsave("P_shoot_root_Hg_concentration_220322.tiff",width = 6,height = 6)


#Hg content in the shoot and root
Hg_root_shoot_content <- Whole_Hg %>%
  select(c(treatment,Above_Content,Root_Content)) %>%
  pivot_longer(cols = c(Above_Content,Root_Content),names_to="plant_Hg",
               values_to="Hg_content")
Hg_root_shoot_content$treatment <- factor(Hg_root_shoot_content$treatment,
                                         levels = c("CK0","AMF0","CK25","AMF25","CK50","AMF50"))
ggplot(Hg_root_shoot_content,aes(x=treatment,y=Hg_content,fill=plant_Hg))+
  geom_bar(stat = "identity",position = position_dodge())+
  scale_x_discrete(labels=c("Hg 0", "Hg 0*","Hg 25","Hg 25*","Hg 50","Hg 50*"))
  

#here i plot leaf,stem and root together
data_concentration_1 <- Whole_Hg%>%
  select(c(treatment,AMF,Hg_C,Leaf_Concentration,Root_Concentration,Stem_Concentration))

data_concentration_2 <- data_concentration_1 %>%
  pivot_longer(cols = c("Leaf_Concentration","Root_Concentration","Stem_Concentration"),names_to="plant_Hg",
               values_to="Hg_concentration")
data_concentration_2$treatment=factor(data_concentration_2$treatment,
  levels = c("CK0","AMF0","CK25","AMF25","CK50","AMF50"))
data_concentration_2$plant_Hg=factor(data_concentration_2$plant_Hg,
                                     levels =c("Leaf_Concentration","Stem_Concentration","Root_Concentration") )
data_concentration_2 <- data_concentration_2%>%
  group_by(treatment,AMF,Hg_C,plant_Hg)%>%
  summarise(Hg_mean=mean(Hg_concentration),Hg_sd=sd(Hg_concentration))%>%ungroup()
data_concentration_2 <- data_concentration_2%>%
  mutate(label=c("b","b","a",
                 "a","a","a",
                 "b","a","a",
                 "a","a","a",
                 "b","a","a",
                 "a","a","a"))
p_root_stem_leaf<-ggplot(data_concentration_2,aes(x=treatment,y=Hg_mean,fill=plant_Hg))+
  geom_errorbar(aes(ymin=Hg_mean-Hg_sd,ymax=Hg_mean+Hg_sd),width=0.2,
                position = position_dodge(0.9))+
  geom_text(aes(x=treatment,y=Hg_mean+Hg_sd,label=label),position = position_dodge(0.9),vjust=-0.3)+
  geom_bar(stat = "identity",position = position_dodge(0.9),color="black")+
  theme_few()+
  scale_x_discrete(label=c("Hg 0","Hg 0*","Hg 25","Hg 25*","Hg 50","Hg 50*"))+
  scale_y_break(c(3,7),scales = 1,expand = F)+
  geom_vline(xintercept = 2.5, color = "blue",size=0.2)+
  geom_vline(xintercept = 4.5, color = "blue",size=0.2)+
  ylab("Hg Concentration (µg/g)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,38))+
  scale_fill_manual(labels=c("leaf","stem","root"),values =c("forestgreen","gold4","gold3"))+
  theme(legend.title = element_blank(),legend.position = "top",
    axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"))
p_root_stem_leaf
ggsave("root_stem_leaf_Hg_concentration.tiff",width = 4,height = 5)

#here is Hg content
data_content_1 <- Whole_Hg%>%
  select(c(treatment,AMF,Hg_C,Leaf_Content,Root_Content,Stem_Content))
data_content_2 <- data_content_1%>%
  pivot_longer(cols=c("Leaf_Content","Root_Content","Stem_Content"),names_to="plant_part",
               values_to="Hg_content")%>%
  group_by(treatment,plant_part)%>%
  summarise(Hg_content_mean=mean(Hg_content),Hg_content_sd=sd(Hg_content))%>%ungroup()
ggplot(data_content_2,aes(x=treatment,y=Hg_content_mean,fill=plant_part))+
  #facet_grid(~AMF)+
  geom_errorbar(aes(ymin=Hg_content_mean-Hg_content_sd,ymax=Hg_content_mean+Hg_content_sd),width=0.2,
                position = position_dodge(0.9))+
  geom_bar(stat = "identity",position = position_dodge(0.9),color="black")+
  theme_few()+
  theme(legend.position = "none",
        axis.title.x = element_blank())+
  scale_x_discrete(label=c("Hg 0","Hg 0*","Hg 25","Hg 25*","Hg 50","Hg 50*"))+
  scale_y_break(c(2,4),scales = 1,expand = F)+
  #scale_y_break(c(10,25),scales = 1)+
  ylab("Hg Content (µg)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,22))+
  scale_fill_manual(values = c("forestgreen","gold4","gold3"))
ggsave("root_stem_leaf_Hg_content.tiff",width = 4,height = 5)


Hg_diff <- read_csv("Hg_diff.csv",col_types = "ffcfffffdddddddddddddddddd")

Hg_diff_0_25 <- Hg_diff%>%
  select(AMF,Hg_C,Stem_Concentration_diff)%>%
  filter(Hg_C=="0"&AMF=="Control"|Hg_C=="25"&AMF=="Control")
t.test(Stem_Concentration_diff ~ Hg_C,data=Hg_diff_0_25,var.equal=T)

Hg_diff_0_25 <- Hg_diff%>%
  select(AMF,Hg_C,Stem_Concentration_diff)%>%
  filter(Hg_C=="0"&AMF=="AMF"|Hg_C=="25"&AMF=="AMF")
t.test(Stem_Concentration_diff ~ Hg_C,data=Hg_diff_0_25,var.equal=T)

Hg_diff_0_50 <- Hg_diff%>%
  select(AMF,Hg_C,Root_Concentration_diff)%>%
  filter(Hg_C=="0"&AMF=="Control"|Hg_C=="50"&AMF=="Control")
t.test(Root_Concentration_diff ~ Hg_C,data=Hg_diff_0_50,var.equal=T)

Hg_diff_0_50 <- Hg_diff%>%
  select(AMF,Hg_C,Root_Concentration_diff)%>%
  filter(Hg_C=="0"&AMF=="AMF"|Hg_C=="50"&AMF=="AMF")
t.test(Root_Concentration_diff ~ Hg_C,data=Hg_diff_0_50,var.equal=T)

Hg_diff1 <- Hg_diff%>%
  select(c(treatment,AMF,Hg_C,Leaf_Concentration_diff,Stem_Concentration_diff,Root_Concentration_diff))%>%
  filter(Hg_C=="25"|Hg_C=="50")

Hg_diff2 <- Hg_diff1%>%
  pivot_longer(cols=c("Leaf_Concentration_diff","Stem_Concentration_diff","Root_Concentration_diff"),
               names_to="plant_diff",values_to="Hg_diff")%>%
  group_by(treatment,plant_diff)%>%
  summarise(Hg_diff_mean=mean(Hg_diff),Hg_diff_sd=sd(Hg_diff))%>%ungroup()
Hg_diff2$treatment=factor(Hg_diff2$treatment,levels = c("CK25","AMF25","CK50","AMF50"))
Hg_diff2$plant_diff=factor(Hg_diff2$plant_diff,levels = c("Leaf_Concentration_diff","Stem_Concentration_diff","Root_Concentration_diff"))

Hg_diff2 <- Hg_diff2%>%
  mutate(label=c("ns","*","*",
                 "ns","***","**",
                 "ns","**","ns",
                 "ns","***","**"))

ggplot(Hg_diff2,aes(x=treatment,y=Hg_diff_mean,fill=plant_diff))+
  geom_errorbar(aes(ymin=0,ymax=Hg_diff_mean+Hg_diff_sd),width=0.2,
                position = position_dodge(0.9))+
  geom_text(aes(x=treatment,y=Hg_diff_mean+Hg_diff_sd,label=label),position = position_dodge(0.9),vjust=-0.3)+
  geom_bar(stat = "identity",position = position_dodge(0.9),color="black")+
  theme_few()+
  scale_x_discrete(label=c("Hg 25","Hg 25*","Hg 50","Hg 50*"))+
  scale_y_break(c(1.3,8),scales = 1,expand = F)+
  #geom_vline(xintercept = 2.5, color = "blue",size=0.2)+
  #geom_vline(xintercept = 4.5, color = "blue",size=0.2)+
  ylab("Hg Concentration difference (µg/g)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,38))+
  scale_fill_manual(labels=c("leaf","stem","root"),values =c("forestgreen","gold4","gold3"))+
  theme(legend.title = element_blank(),legend.position = "top",
        axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.ticks.length.y = unit(-1.2, "mm"))
ggsave("Hg_concentration_difference.tiff",height = 5,width = 4)


#here i want to put above and below 




