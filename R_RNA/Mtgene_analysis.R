library(tidyverse)
library(emmeans)
library(ggthemes)
library(car)
library(ggtext)
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
ZIP2<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "ZIP2")
ZIP2$AMF <- as.factor(ZIP2$AMF)
ZIP2$Hg <- as.factor(ZIP2$Hg)

car::leveneTest(ZIP2 ~ Hg*AMF,data=ZIP2)#p=0.867
m_ZIP2 <- aov(ZIP2 ~ Hg*AMF,data=ZIP2)
summary(m_ZIP2)
shapiro.test(residuals(m_ZIP2))#p-value = 0.1848
tukey.ZIP2 <- TukeyHSD(m_ZIP2)
multcompView::multcompLetters4(m_ZIP2,tukey.ZIP2)
ZIP2_data <- data_summary(ZIP2,varname = "ZIP2",
                                   groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("a","ab","ab","abc","c","bc"))
ZIP2_data$AMF <- factor(ZIP2_data$AMF,levels = c("Control","AMF"))
ZIP2_data$Hg <- as.factor(ZIP2_data$Hg)
write_csv(ZIP2_data,"ZIP2_summary.csv")
P_ZIP2 <- ggplot(ZIP2_data,aes(x=Hg,y=ZIP2,fill=AMF))+
  geom_errorbar(aes(ymin=ZIP2-sd,ymax=ZIP2+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.6,color="black")+
  geom_text(aes(x=Hg,y=ZIP2+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(axis.title.x = element_blank(),axis.title.y = element_markdown(),
        plot.margin = margin(0,0,0,0,unit = "pt"),
    legend.position ="none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y="Relative expression of *ZIP2*")+
  geom_richtext(x=1,y=4.8, label="*P*(AMF)<0.001<br>*P*(Hg)=0.60<br>*P*(AMFXHg)=0.90",fill=NA,label.color=NA,
                hjust=0,vjust=0,size=3,family = "Arial",fontface="plain")+
  geom_text(x=3.4,y=5.4,label="(A)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,5.5))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values=c('white','grey80'))
P_ZIP2 
ggsave("P_ZIP2_new_bw_250322.tiff",width = 5,heigh=4)
#ZIP6
ZIP6<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "ZIP6")
ZIP6$AMF <- as.factor(ZIP6$AMF)
ZIP6$Hg <- as.factor(ZIP6$Hg)

leveneTest(ZIP6 ~ Hg*AMF,data=ZIP6)#p=0.04
m_ZIP6 <- aov(ZIP6 ~ Hg*AMF,data=ZIP6)#1.409e-09 ***
summary(m_ZIP6)
shapiro.test(residuals(m_ZIP6))#p-value = 0.5342
tukey.ZIP6 <- TukeyHSD(m_ZIP6)
multcompView::multcompLetters4(m_ZIP6,tukey.ZIP6)
ZIP6_data <- data_summary(ZIP6,varname = "ZIP6",
                          groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("b","b","b","a","a","a"))
ZIP6_data$AMF <- factor(ZIP6_data$AMF,levels = c("Control","AMF"))
ZIP6_data$Hg <- as.factor(ZIP6_data$Hg)
write_csv(ZIP6_data,"ZIP6_summary.csv")
P_ZIP6 <- ggplot(ZIP6_data,aes(x=Hg,y=ZIP6,fill=AMF))+
  geom_errorbar(aes(ymin=ZIP6-sd,ymax=ZIP6+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.6,color="black")+
  geom_text(aes(x=Hg,y=ZIP6+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(axis.title.x = element_blank(),axis.title.y.right = element_markdown(angle = 90),
        plot.margin = margin(0,0,0,0,unit = "pt"),
    legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y="Relative expression of *ZIP6*")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  geom_richtext(x=1,y=48, label="*P*(AMF)<0.001<br>*P*(Hg)=0.26<br>*P*(AMFXHg)=0.71",fill=NA,label.color=NA,
                hjust=0,vjust=0,size=3,family = "Arial",fontface="plain")+
  geom_text(x=3.4,y=54,label="(B)")+
  scale_y_continuous(expand = c(0,0),limits = c(0,55),position = "right")+
  scale_fill_manual(values=c('white','grey80'))
P_ZIP6
ggsave("P_ZIP6_bw_250322.tiff",width = 5,heigh=4)
library(patchwork)
P_ZIP2|P_ZIP6
ggsave("p_ZIP2_ZIP6.tiff",width = 6,height = 6)

#MtNAS1
MtNAS1<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "MtNAS1")
MtNAS1$AMF <- as.factor(MtNAS1$AMF)
MtNAS1$Hg <- as.factor(MtNAS1$Hg)

MtNAS1_50 <- MtNAS1 %>%
  filter(Hg=="50")
aov_MtNAS1_50 <- aov(MtNAS1 ~ AMF,MtNAS1_50)
Anova(aov_MtNAS1_50,type="III")

MtNAS1_Control <- MtNAS1 %>%
  filter(AMF=="Control")
aov_MtNAS1_Control <- aov(MtNAS1 ~ Hg,MtNAS1_Control)
Anova(aov_MtNAS1_Control,type="III")


car::leveneTest(MtNAS1 ~ Hg*AMF,data=MtNAS1)#p=0.61
shapiro.test(MtNAS1$MtNAS1)#p-value = 0.17
m_MtNAS1 <- aov(MtNAS1 ~ Hg+AMF,data=MtNAS1)
Anova(m_MtNAS1,type="III")#AMF 0.001392 ** 
shapiro.test(residuals(m_MtNAS1))#p=0.041
tukey.cld_MtNAS1 <- emmeans(m_MtNAS1,specs = pairwise ~ AMF+Hg) %>% 
  multcomp::cld(Letters=letters)
tukey.cld_MtNAS1
MtNAS1_data <- data_summary(MtNAS1,varname = "MtNAS1",
                          groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("bc","bc","c","ab","ab","ab"))
MtNAS1_data$AMF <- factor(MtNAS1_data$AMF,levels = c("Control","AMF"))
MtNAS1_data$Hg <- as.factor(MtNAS1_data$Hg)
write_csv(MtNAS1_data,"MtNAS1_summary.csv")
P_MtNAS1 <- ggplot(MtNAS1_data,aes(x=Hg,y=MtNAS1,fill=AMF))+
  geom_errorbar(aes(ymin=MtNAS1-sd,ymax=MtNAS1+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg,y=MtNAS1+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(x="Hg_Concentration (mg/kg)",y="Relative expression of MtNAS1")+
  geom_text(x=3.2, y=3.9, label="AMF**")+
  scale_y_continuous(expand = c(0,0),limits = c(0,4))+
  scale_fill_manual(values=c('#999999','#E69F00'))
P_MtNAS1
ggsave("P_MtNAS1_new.tiff",width = 5,heigh=4)
#MtPCs
MtPCs<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "MtPCs")
MtPCs$AMF <- as.factor(MtPCs$AMF)
MtPCs$Hg <- as.factor(MtPCs$Hg)

MtPCs_25 <- MtPCs %>%
  filter(Hg=="25")
aov_MtPCs_25 <- aov(MtPCs ~ AMF, MtPCs_25)
Anova(aov_MtPCs_25,type="III")

MtPCs_Control <- MtPCs %>%
  filter(AMF=="Control")
aov_MtPCs_Control <- aov(MtPCs ~ Hg, MtPCs_Control)
Anova(aov_MtPCs_Control,type="III")

car::leveneTest(MtPCs ~ Hg*AMF,data=MtPCs)#p=0.84
m_MtPCs <- aov(MtPCs ~ Hg+AMF,data=MtPCs)
Anova(m_MtPCs,type="III") #no effect
shapiro.test(residuals(m_MtPCs))#p=0.6
MtPCs_data <- data_summary(MtPCs,varname = "MtPCs",
                            groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("a","a","a","a","a","a"))
MtPCs_data$AMF <- factor(MtPCs_data$AMF,levels = c("Control","AMF"))
MtPCs_data$Hg <- as.factor(MtPCs_data$Hg)
write_csv(MtPCs_data,"MtPCs_summary.csv")
P_MtPCs <- ggplot(MtPCs_data,aes(x=Hg,y=MtPCs,fill=AMF))+
  geom_errorbar(aes(ymin=MtPCs-sd,ymax=MtPCs+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg,y=MtPCs+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(
    y="Relative expression of MtPCs")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.6))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values=c('white','grey80'))
P_MtPCs
ggsave("P_MtPCs_bw.tiff",width = 5,heigh=4)
#MtRECS
MtRECS<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "MtRECS")
MtRECS$AMF <- as.factor(MtRECS$AMF)
MtRECS$Hg <- as.factor(MtRECS$Hg)
car::leveneTest(MtRES ~ Hg*AMF,data=MtRECS)#p=0.7531
shapiro.test(MtRECS$MtRES)#p-value = 0.5614
m_MtRECS <- aov(MtRES ~ Hg*AMF,data=MtRECS)
Anova(m_MtRECS,type="III")
shapiro.test(residuals(m_MtRECS))#p-value = 0.2136 
tukey.cld_MtRECS <- emmeans(m_MtRECS,specs = pairwise ~ AMF+Hg) %>% 
  multcomp::cld(Letters=letters)
tukey.cld_MtRECS
MtRECS_data <- data_summary(MtRECS,varname = "MtRES",
                            groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("ab","a","a","ab","b","a"))
MtRECS_data$AMF <- factor(MtRECS_data$AMF,levels = c("Control","AMF"))
MtRECS_data$Hg <- as.factor(MtRECS_data$Hg)
write_csv(MtRECS_data,"MtRECS_summary.csv")
P_MtRECS <- ggplot(MtRECS_data,aes(x=Hg,y=MtRES,fill=AMF))+
  geom_errorbar(aes(ymin=MtRES-sd,ymax=MtRES+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  geom_text(aes(x=Hg,y=MtRES+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(legend.position = c(0.1,0.92),legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  geom_text(x=3.2, y=1.9, label="Hg X AMF*")+
  labs(x="Hg_Concentration (mg/kg)",y="Relative expression of MtRES")+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))+
  scale_fill_manual(values=c('#999999','#E69F00'))
P_MtRECS
ggsave("P_MtRECS.tiff",width = 5,heigh=4)
#MtSODa
MtSODa<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "MtSODa")
MtSODa$AMF <- as.factor(MtSODa$AMF)
MtSODa$Hg <- as.factor(MtSODa$Hg)

MtSODa_25 <- MtSODa %>%
  filter(Hg=="25")
aov_MtSODa_25 <- aov(MtSODa ~ AMF,data = MtSODa_25 )
Anova(aov_MtSODa_25,type="III")

MtSODa_Control <- MtSODa %>%
  filter(AMF=="Control")
aov_MtSODa_Control <- aov(MtSODa ~ Hg,data = MtSODa_Control )
Anova(aov_MtSODa_Control,type="III")


car::leveneTest(MtSODa ~ Hg*AMF,data=MtSODa)#p=0.94
shapiro.test(MtSODa$MtSODa)#p-value = 0.35
m_MtSODa <- aov(MtSODa ~ Hg+AMF,data=MtSODa)
Anova(m_MtSODa,type="III")#no effect
shapiro.test(residuals(m_MtSODa))#p=0.3778
MtSODa_data <- data_summary(MtSODa,varname = "MtSODa",
                           groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("a","a","a","a","a","a"))
MtSODa_data$AMF <- factor(MtSODa_data$AMF,levels = c("Control","AMF"))
MtSODa_data$Hg <- as.factor(MtSODa_data$Hg)
write_csv(MtSODa_data,"MtSODa_summary.csv")
P_MtSODa <- ggplot(MtSODa_data,aes(x=Hg,y=MtSODa,fill=AMF))+
  geom_errorbar(aes(ymin=MtSODa-sd,ymax=MtSODa+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  #geom_text(aes(x=Hg,y=MtSODa+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  labs(y="Relative expression of MtCuZnSODa")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.5))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_fill_manual(values=c('white','grey80'))
P_MtSODa
ggsave("P_MtSODa_bw.tiff",width = 5,heigh=4)
#MtSODb
MtSODb<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "MtSODb")
MtSODb$AMF <- as.factor(MtSODb$AMF)
MtSODb$Hg <- as.factor(MtSODb$Hg)
MtSODb_25 <- MtSODb %>%
  filter(Hg=="25")
aov_MtSODa_25 <- aov(MtSODb ~ AMF, MtSODb_25)
Anova(aov_MtSODa_25,type="III")

MtSODb_AMF <- MtSODb %>%
  filter(AMF=="AMF")
aov_MtSODa_AMF <- aov(MtSODb ~ Hg, MtSODb_AMF)
Anova(aov_MtSODa_AMF,type="III")

car::leveneTest(MtSODb ~ Hg*AMF,data=MtSODb)#p=0.46
shapiro.test(MtSODb$MtSODb)#p=0.69
m_MtSODb <- aov(MtSODb ~ Hg+AMF,data=MtSODb)
Anova(m_MtSODb,type="III")
shapiro.test(residuals(m_MtSODb))#p=0.65
tukey.cld_MtSODb <- emmeans(m_MtSODb,specs = pairwise ~ AMF+Hg) %>% 
  multcomp::cld(Letters=letters)
tukey.cld_MtSODb
MtSODb_data <- data_summary(MtSODb,varname = "MtSODb",
                            groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("a","a","a","a","a","a"))
MtSODb_data$AMF <- factor(MtSODb_data$AMF,levels = c("Control","AMF"))
MtSODb_data$Hg <- as.factor(MtSODb_data$Hg)
write_csv(MtSODb_data,"MtSODb_summary.csv")
P_MtSODb <- ggplot(MtSODb_data,aes(x=Hg,y=MtSODb,fill=AMF))+
  geom_errorbar(aes(ymin=MtSODb-sd,ymax=MtSODb+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  #geom_text(aes(x=Hg,y=MtSODb+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  #geom_text(x=3.2, y=1.9, label="AMF*")+
  labs(y="Relative expression of MtCuZnSODb")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.6))+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 25"))+
  scale_fill_manual(values=c('white',"grey80"))
P_MtSODb
ggsave("P_MtSODb_new_bw.tiff",width = 5,heigh=4)
#MtSODc
MtSODc<-readxl::read_xlsx("Targetgene_figure_O_delete_R.xlsx",sheet = "MtSODc")
MtSODc$AMF <- as.factor(MtSODc$AMF)
MtSODc$Hg <- as.factor(MtSODc$Hg)
MtSODc_50 <- MtSODc %>%
  filter(Hg=="50")
aov_MtSODc_50 <- aov(MtSODc ~ AMF, MtSODc_50)
Anova(aov_MtSODc_50, type="III")

MtSODc_Control <- MtSODc %>%
  filter(AMF=="Control")
aov_MtSODc_Control <- aov(MtSODc ~ Hg, MtSODc_Control)
Anova(aov_MtSODc_Control, type="III")








car::leveneTest(MtSODc ~ Hg*AMF,data=MtSODc)#p=0.93
shapiro.test(MtSODc$MtSODc)#p=0.03
m_MtSODc <- aov(MtSODc ~ Hg+AMF,data=MtSODc)
Anova(m_MtSODc,type="III")
shapiro.test(residuals(m_MtSODc))#p-value = 0.05505
tukey.cld_MtSODc <- emmeans(m_MtSODc,specs = pairwise ~ AMF+Hg) %>% 
  multcomp::cld(Letters=letters)
tukey.cld_MtSODc
MtSODc_data <- data_summary(MtSODc,varname = "MtSODc",
                            groupnames = c("AMF","Hg")) %>%
  mutate(tukey=c("b","b","b","a","a","a"))
MtSODc_data$AMF <- factor(MtSODc_data$AMF,levels = c("Control","AMF"))
MtSODc_data$Hg <- as.factor(MtSODc_data$Hg)
write_csv(MtSODc_data,"MtSODc_summary.csv")
P_MtSODc <- ggplot(MtSODc_data,aes(x=Hg,y=MtSODc,fill=AMF))+
  geom_errorbar(aes(ymin=MtSODc-sd,ymax=MtSODc+sd),width=0.2,
                position = position_dodge(0.7))+
  geom_bar(position = position_dodge(0.7),stat = "identity",
           width = 0.5,color="black")+
  #geom_text(aes(x=Hg,y=MtSODc+sd,label=tukey),position=position_dodge(0.7),vjust = -0.5)+
  theme_few(base_size = 10,base_family = "Arial")+
  theme(aspect.ratio = 1,axis.title.x = element_blank(),
    legend.position = "none",legend.background = element_blank(),legend.title = element_blank(),
        text = element_text(size=10,family="Arial"))+
  geom_text(x=1, y=2.2, label="*")+
  geom_text(x=2, y=2.2, label="*")+
  geom_text(x=3, y=2.2, label="ns")+
  labs(y="Relative expression of MtCuZnSODc")+
  scale_x_discrete(labels=c("Hg 0","Hg 25","Hg 50"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,2.4))+
  scale_fill_manual(values=c('white','grey80'))
P_MtSODc
ggsave("P_MtSODc_bw.tiff",width = 5,heigh=4)
