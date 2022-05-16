library(tidyverse)
library(ggplot2)
library(ggthemes)
library(multcompView)
library(ggtext)
Ri_gene <- readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",na="NA")
Ri_gene 
Ri_gene$Hg <- as.factor(Ri_gene$Hg)
#ri_tublin
car::leveneTest(Ri_tublin ~ Hg, data=Ri_gene)
shapiro.test(Ri_gene$Ri_tublin)# p-value = 0.0005901
m_ri_tublin <- aov(Ri_tublin ~ Hg, data=Ri_gene)
car::Anova(m_ri_tublin,type="III")
shapiro.test(residuals(m_ri_tublin))#p=0.16
tukey.ri_tublin <- TukeyHSD(m_ri_tublin)
multcompLetters4(m_ri_tublin,tukey.ri_tublin)

P_ri_tublin<- ggplot(Ri_gene,aes(x=Hg, y=Ri_tublin))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(width=0.5)+
  labs(x="Hg_Concentraion (mg/kg)",y="Relative expression of Ri_tublin")+
  theme_few(base_size = 10,base_family = "Arial")
P_ri_tublin <- P_ri_tublin+stat_summary(fun = mean,geom = "point",shape=21,size=2,fill="red")
P_ri_tublin
ggsave("Ri_tublin_geneexpressionration_boxplot.tiff",width = 5,height = 4)
ri_tublin <- Ri_gene %>%
  group_by(Hg)%>%
  summarise(mean_tublin=mean(Ri_tublin,na.rm = T),sd_tublin=sd(Ri_tublin,na.rm = T))
P_ri_tublin<- ggplot(ri_tublin,aes(x=Hg, ymin=mean_tublin,ymax=mean_tublin+sd_tublin))+
  geom_errorbar(aes(x=Hg,y=mean_tublin),width=0.1)+
  geom_bar(aes(x=Hg, y=mean_tublin),stat = "identity",fill="grey80",color="black",width = 0.3)+
  geom_line(aes(x=Hg, y=mean_tublin),group=1,color="Red",linetype=2)+
  geom_point(aes(x=Hg, y=mean_tublin),fill="red",shape=21,size=3)+
  labs(y="Relative expression of *Ri tublin*")+
  geom_text(x=3.4,y=4.85,label="(B)",size=3)+
  geom_richtext(x=3,y=4.5,label="*P*(Hg)=0.18",fill=NA,label.color=NA,size=3)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(plot.margin = margin(0,0,0,0,"pt"),
    axis.title.x = element_blank(),axis.title.y.right = element_markdown())+
  scale_y_continuous(expand = c(0,0),limits = c(0,5),position = "right")+
  scale_x_discrete(breaks=c("0","25","50"),labels=c("Hg 0","Hg 25","Hg 50"))
P_ri_tublin  
ggsave("Ri_tublin_geneexpressionration_barwithlineplot.tiff",width = 5,height = 4)

#MtZIP14
car::leveneTest(MtZIP14 ~ Hg, data=Ri_gene)#p=0.8
shapiro.test(Ri_gene$MtZIP14)#p-value = 0.07
Aov_MtZIP14 <- aov(MtZIP14 ~ Hg, data=Ri_gene)
car::Anova(Aov_MtZIP14,type="III")
shapiro.test(residuals(Aov_MtZIP14))#p=0.228
Tukey_MtZIP14 <- TukeyHSD(Aov_MtZIP14)
cld_MtZIP14 <- multcompLetters4(Aov_MtZIP14, Tukey_MtZIP14)
cld_MtZIP14
pairwise.wilcox.test(Ri_gene$MtZIP14,Ri_gene$Hg,p.adjust.method = "holm")
P_MtZIP14<- ggplot(Ri_gene,aes(x=Hg, y=MtZIP14))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(width=0.5)+
  labs(x="Hg_Concentraion (mg/kg)",y="Relative expression of MtZIP14")+
  geom_text(x=1,y=1.12,label = "a",size=4,family="Arial")+
  geom_text(x=2,y=0.38,label = "b",size=4,family="Arial")+
  geom_text(x=3,y=0.23,label = "b",size=4,family="Arial")+
  theme_few(base_size = 10,base_family = "Arial")
P_MtZIP14+stat_summary(fun = mean,geom = "point",shape=21,size=2,fill="red")
ggsave("MtZIP14_geneexpressionration.tiff",width = 5,height = 4)
detach(package:plyr)
ri_MtZIP14 <- Ri_gene %>%
  group_by(Hg)%>%
  summarise(mean_MtZIP14=mean(MtZIP14,na.rm = T),sd_MtZIP14=sd(MtZIP14,na.rm = T))
P_ri_MtZIP14<- ggplot(ri_MtZIP14,aes(x=Hg, ymin=mean_MtZIP14,ymax=mean_MtZIP14+sd_MtZIP14))+
  geom_errorbar(aes(x=Hg,y=mean_MtZIP14),width=0.1)+
  geom_bar(aes(x=Hg, y=mean_MtZIP14),stat = "identity",fill="grey90",color="black",width = 0.5)+
  #geom_line(aes(x=Hg, y=mean_MtZIP14),group=1,color="Red",linetype=2)+
  #geom_point(aes(x=Hg, y=mean_MtZIP14),fill="red",shape=21,size=3)+
  labs(y="Relative expression of MtZIP14")+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.3))+
  geom_text(x=1,y=1.17,label="b")+
  geom_text(x=2,y=0.56,label="a")+
  geom_text(x=3,y=0.22,label="a")+
  scale_x_discrete(breaks=c("0","25","50"),labels=c("Hg 0","Hg 25","Hg 50"))
P_ri_MtZIP14  
ggsave("MtZIP14_geneexpressionration_bw.tiff",width = 5,height = 4)
#MtPT4
car::leveneTest(MtPT4 ~ Hg, data=Ri_gene)#p=0.1757
shapiro.test(Ri_gene$MtPT4)# p-value = 0.2139
# MtPT4_data <- Ri_gene %>%
#   select(c(4,6)) %>%
#   group_by(Hg) %>%
#   summarise_all(funs(mean,sd),na.rm=T)
m_MtPT4 <- aov(MtPT4 ~ Hg, data = Ri_gene)
car::Anova(m_MtPT4,type="III")
TukeyHSD(m_MtPT4)
P_MtPT4<- ggplot(Ri_gene,aes(x=Hg, y=MtPT4))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(width=0.5)+
  labs(x="Hg_Concentraion (mg/kg)",y="Relative expression of MtPT4")+
  theme_few(base_size = 10,base_family = "Arial")
P_MtPT4+stat_summary(fun = mean,geom = "point",shape=21,size=2,fill="red")
ggsave("MtPT4_geneexpressionration.tiff",width = 5,height = 4)

ri_MtPT4 <- Ri_gene %>%
  group_by(Hg)%>%
  summarise(mean_MtPT4=mean(MtPT4,na.rm = T),sd_MtPT4=sd(MtPT4,na.rm = T))

P_ri_MtPT4<- ggplot(ri_MtPT4,aes(x=Hg, ymin=mean_MtPT4,ymax=mean_MtPT4+sd_MtPT4))+
  geom_errorbar(aes(x=Hg,y=mean_MtPT4),width=0.1)+
  geom_bar(aes(x=Hg, y=mean_MtPT4),stat = "identity",fill="grey90",color="black",width = 0.5)+
  #geom_line(aes(x=Hg, y=mean_MtZIP14),group=1,color="Red",linetype=2)+
  #geom_point(aes(x=Hg, y=mean_MtZIP14),fill="red",shape=21,size=3)+
  labs(y="Relative expression of MtPT4")+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.8))+
  scale_x_discrete(breaks=c("0","25","50"),labels=c("Hg 0","Hg 25","Hg 50"))
P_ri_MtPT4  
ggsave("MtPT4_geneexpressionration_bw.tiff",width = 5,height = 4)


##using the following code could get thr correclation, but none of them have significant correclation
library("ggpubr")
ggscatter(Ri_gene, x = "MtPT4", y = "MtZIP14", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "", ylab = "")
#colonization
Colonization <- readxl::read_xlsx("Colonization_AMF_R.xlsx",na="NA")
Colonization
Colonization$Hg <- as.factor(Colonization$Hg)
#comapre each with zero
with(Colonization,shapiro.test(Percent[Hg=="0"]))
with(Colonization,shapiro.test(Percent[Hg=="25"]))
with(Colonization,shapiro.test(Percent[Hg=="50"]))
data1 <- Colonization %>%
  filter(Hg=="0"|Hg=="25")
var.test(Percent ~ Hg,data = data1)
data2 <- Colonization %>%
  filter(Hg=="0"|Hg=="50")
var.test(Percent ~ Hg,data = data2)
res1 <- t.test(Percent ~ Hg,data = data1)
res1
res2 <- t.test(Percent ~ Hg,data = data2)
res2

car::leveneTest(Percent ~ Hg, data=Colonization)#p=0.795
shapiro.test(Colonization$Percent)#0.612
aov_colonization <- aov(Percent ~ Hg, data=Colonization)
car::Anova(aov_colonization,type="III")
shapiro.test(residuals(aov_colonization))
Tukey_colonization <- TukeyHSD(aov_colonization)
cld_colonization <- multcompLetters4(aov_colonization, Tukey_colonization)
P_colonization <- ggplot(Colonization,aes(x=Hg, y=Percent))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(width=0.5)+
  labs(x="Hg_Concentraion (mg/kg)",y="Mycorrhizal Colonization (%)")+
  geom_text(x=1,y=59,label = "ab",size=4,family="Arial")+
  geom_text(x=2,y=58,label = "a",size=4,family="Arial")+
  geom_text(x=3,y=38,label = "b",size=4,family="Arial")+
  theme_few(base_size = 10,base_family = "Arial")
P_colonization <- P_colonization+stat_summary(fun = mean,geom = "point",shape=21,size=2,fill="red")
P_colonization
ggsave("AMF_colonization.tiff",width = 5,height = 4)
ggarrange(P_colonization,P_ri_tublin,labels = c("a","b"))
ggsave("P_AMF_colonization_ritublin .tiff",width = 8,height = 4)
colonization <-Colonization %>%
  group_by(Hg)%>%
  summarise(mean_colonization=mean(Percent,na.rm = T),sd=sd(Percent,na.rm = T))
P_colonization<- ggplot(colonization,aes(x=Hg, y=mean_colonization))+
  geom_errorbar(aes(x=Hg,ymin=mean_colonization-sd-0.1,ymax=mean_colonization+sd+0.1),width=0.1)+
  geom_bar(aes(x=Hg, y=mean_colonization),stat = "identity",fill="grey80",color="black",width = 0.3)+
  geom_line(aes(x=Hg, y=mean_colonization),group=1,color="Red",linetype=2)+
  geom_point(fill="red",shape=21,size=3)+
  geom_richtext(x=3,y=72,label="*P*(Hg)=0.03",fill=NA,label.color=NA,size=3)+
  labs(y="Mycorrhizal colonization (%)")+
  geom_text(aes(x=Hg,y=mean_colonization+sd+1.5,label =c("ab","a","b")))+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,
        plot.margin = margin(0,0,0,0,unit = "pt"),
        axis.title.x = element_blank())+
  geom_text(x=3.4,y=78,label="(A)",size=3)+
  scale_y_continuous(expand = c(0,0),limits = c(0,80))+
  scale_x_discrete(breaks=c("0","25","50"),labels=c("Hg 0","Hg 25","Hg 50"))
P_colonization 
ggsave("AMF_colonization_barwithlineplot.tiff",width = 5,height = 4)
library(patchwork)
P_colonization|P_ri_tublin
ggsave("P_colonization_ri_tublin.tiff",width=8,height = 6)
