##Effect of diverse ions concentrations on eDNA detectability as yield of spiked DNA detected.

# Packages
library(readxl)
library(xlsx)
library(ggplot2)
library(tidyverse)
library(drc)
library(MASS)
library(scales)
library(rstatix)
library(rstatix)
library(car)
library(FSA)
library(rcompanion)

#Load Dataset
Data_1 <- read_excel("Metal_ions_qPCR_data.xlsx", sheet = "Calcium")
#Change sheet for each metal ions (Copper, Calcium, Manganese and Iron)

#####
Data_l.1<-Data_1 %>% group_by(sample_trt, biological_rep) %>% get_summary_stats(Cq) 
DATATEST <- Data_l.1 %>% group_by(sample_trt) %>% get_summary_stats(mean)

#Add columns 
Data_2.1<-Data_1 %>% group_by(sample_trt, biological_rep) %>% mutate(moyCq=mean(Cq))

Data_2.2 <-Data_2.1%>% mutate(copy.t= 10^((moyCq -37.6511857678726)/(-3.2915910766764)),
                              copy.max = 150, copy.freq = copy.t/copy.max, Cq.max = 45, Cq.Cqmax= Cq/Cq.max)
Data_3<-Data_2.2 %>% group_by(conc_sed) %>% get_summary_stats(copy.freq)


#Graph
x.break <- sort(unique(Data_2.2$conc_sed))
x.break

######Change x.name for each metal ions######

##Calcium
x.name<- c("110", "115","","150","320", "530", "4,270","20,890")

##Copper
#x.name<- c("0,000","0,004","0,020","0,040","0,080", "0,200", "0,300","0,400")

##Manganese
#x.name<- c("0","4","20","40","200", "400", "4,000","20,000")

##Iron
#x.name<- c("0", "40", "400","2,000", "4,000", "10,000","20,000", "30,000")

x.name

#Change de xlab for each metal ions
Metal_ions_plot<-ggplot() +
  geom_point(data= Data_2.2, aes(x = conc_sed, y = (copy.freq*100)), alpha=0) +
  # geom_line(data=Data3, aes(x=conc_sed, y=mean*100), linetype=5, size=0.7, color="grey50") +
  coord_trans(x="log") +
  geom_point(data=Data_3, aes(x = conc_sed, y=mean*100))+
  geom_errorbar(data=Data_3, aes(x = conc_sed, y=mean*100, ymin=(mean-se)*100, ymax=(mean+se)*100), width=0)+
  scale_y_continuous(breaks=seq(0,10,2), limits=c(0,10), expand=c(0,0)) +
  scale_x_continuous(breaks= x.break, labels=x.name, limits=c(0.05,20.11))+
  xlab("Ca  (mg Kg sediments dry weight)") + ylab("eDNA detectability (%)") +
  theme_classic(base_size = 18) + theme(aspect.ratio = 1) +theme(axis.text.x=element_text(face="plain", angle=90, hjust=0.95,vjust=0.5, size=14))
Metal_ions_plot

#### Statistical analysis ####
Data_4 <-Data_2.2%>% mutate(copy.freq2= (copy.freq+0.001), copy.t2=(log(copy.freq2)))

#Normality and homogenity

shapiro.test(Data_4$copy.t2)

bartlett.test(Data_4$copy.t2~Data_4$sample_trt)
leveneTest(Data_4$copy.t2~Data_4$sample_trt)

results.anova<- aov(Data_4$copy.t2~Data_4$sample_trt)
summary(results.anova)
str(results.anova)
shapiro.test(results.anova$residuals)
plot(results.anova)

anova(results.anova)
TUKEY<-tukey_hsd(results.anova)
summary(TUKEY)
str(TUKEY)

#For all metal ions, non‐normality (Shapiro‐Wilks test) of the data and heteroscedasticity of residuals (Levene test)
#in the analysis of variance (ANOVA) was observed. Consequently, we performed non‐parametric tests (Kruskal‐Wallis)
# to determine differences in DNA detectability. 

# kruskal walis
kruskal.test(Data_2.2$copy.freq~Data_2.2$sample_trt)

#p<0.05 

DT<-dunnTest(Data_2.2$copy.freq~Data_2.2$sample_trt, method = "bh")
str(DT)
DT

Pt<-DT$res
List<-cldList(P.adj~Comparison, data = Pt, threshold = 0.05)






