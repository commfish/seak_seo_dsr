#############################################################################
## Biological data fro yelloweye rockfish in the SEO
## 
## Phil Joy
#############################################################################

{library(tidyverse)
library(ggridges)
library(ggplot2)
library(fishmethods)
library(broom)
library(mosaic)
library(ggpubr)
library(cowplot)
source("r_helper/Port_bio_function.R")
}

YEAR<-2024

#setwd("D:/Groundfish Biometrics/Yelloweye/SAFE reports")

{ Port1<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_1980-1989.csv")
  Port2<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_1990-1999.csv")
  Port3<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_2000-2009.csv")
  Port4<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_2010-2019.csv")
  Port5<-read.csv("Data_processing/Data/SEO_YE_port_sampling_bio_data_2020-2024.csv")
  
  #This report has way more fields now - the original report has project code twice - 
  #edited the oceanAK report to match Ports 1-4 saved in my personal folder in OceanAK.
  
  names(Port5)
  names(Port4)
  
  Port<-rbind(Port1,Port2,Port3, Port4, Port5); nrow(Port)
  #str(Port)
  Port$Year<-as.integer(Port[,1])
  #str(Port)
    #EYAK = EYKT
  unique(Port$Groundfish.Stat.Area)
  unique(Port$Groundfish.Stat.Area.Group)
  unique(Port$Groundfish.Management.Area.Code)
  
  # Port$Groundfish.Management.Area.Code[Port$Groundfish.Management.Area.Code == "EYAK"]<-"EYKT"
  
  statareas<-read.csv("Data_processing/Data/g_stat_area.csv")
  str(statareas)
  unique(statareas$G_MANAGEMENT_AREA_CODE)
  
  SSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEO"])
  CSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "CSEO"])
  NSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEO"])
  EYc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "EYKT"])
  
{  #with(Port, table(Groundfish.Stat.Area,Groundfish.Management.Area.Code))
  #SScodes<-unique(Port$Groundfish.Stat.Area[Port$Groundfish.Management.Area.Code == "SSEO" |
  #                                            substr(Port$Groundfish.Stat.Area.Group,1,4) == "SSEO"])
  #CScodes<-unique(Port$Groundfish.Stat.Area[Port$Groundfish.Management.Area.Code == "CSEO" |
  #                                            substr(Port$Groundfish.Stat.Area.Group,1,4) == "CSEO"])
  #NScodes<-unique(Port$Groundfish.Stat.Area[Port$Groundfish.Management.Area.Code == "NSEO" |
  #                                            substr(Port$Groundfish.Stat.Area.Group,1,4) == "NSEO"])
  #EYcodes<-unique(Port$Groundfish.Stat.Area[Port$Groundfish.Management.Area.Code == "EYKT" |
  #                                            substr(Port$Groundfish.Stat.Area.Group,1,4) == "EYKT"])
  
  #SScodes<-SScodes[!is.na(SScodes)]
  #CScodes<-CScodes[!is.na(CScodes)]
  #NScodes<-NScodes[!is.na(NScodes)]
  #EYcodes<-EYcodes[!is.na(EYcodes)]
  
  #Port %>% mutate(GFMU = ifelse(Groundfish.Management.Area.Code == "",
  #                              ifelse(Groundfish.Stat.Area %in% c(SScodes) |
  #                                       substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("SSEO"),"SSEO",
  #                                     ifelse(Groundfish.Stat.Area %in% c(CScodes) |
  #                                              substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("CSEO"),"CSEO",
  #                                            ifelse(Groundfish.Stat.Area %in% c(NScodes) |
  #                                                     substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("NSEO"),"NSEO",
  #                                                   ifelse(Groundfish.Stat.Area %in% c(EYcodes) |
  #                                                            substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("EYKT"),"EYKT",NA)))),
  #                              Groundfish.Management.Area.Code)) -> Port.xx
  }

  Port %>% mutate(GFMU = ifelse(Groundfish.Management.Area.Code == "",
                                ifelse(Groundfish.Stat.Area %in% c(SSc) |
                                         substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("SSEO"),"SSEO",
                                       ifelse(Groundfish.Stat.Area %in% c(CSc) |
                                                substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("CSEO"),"CSEO",
                                              ifelse(Groundfish.Stat.Area %in% c(NSc) |
                                                       substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("NSEO"),"NSEO",
                                                     ifelse(Groundfish.Stat.Area %in% c(EYc) |
                                                              substr(Port$Groundfish.Stat.Area.Group,1,4) %in% c("EYKT"),"EYKT",NA)))),
                                Groundfish.Management.Area.Code)) -> Port
  
  unknowns<-Port[is.na(Port$GFMU),]
  knowns<-Port[!is.na(Port$GFMU),]
  nrow(unknowns); nrow(knowns)
  unique(unknowns$Groundfish.Management.Area.Code)
  unique(unknowns$Groundfish.Stat.Area.Group)
  unique(unknowns$Groundfish.Stat.Area)
  
  Port<-Port %>% mutate(Sex = 
                          case_when(Sex.Code == 1 ~ "Male",
                                    Sex.Code == 2 ~ "Female")
  )
  Port$Sex<-as.factor(Port$Sex)	
  unique(Port$Sample.Type)
  unique(Port$Project)
  Port.rand<-Port[Port$Sample.Type=="Random",]
  Port.direct<-Port[Port$Project=="Commercial Longline Trip",]
  Port.hal<-Port[Port$Project=="Commercial Halibut Longline",]; nrow(Port.hal)
}

Port<-port.bio(YEAR)
#All random sample from all of SEAK
Port.rand<-Port[Port$Sample.Type=="Random",]

agecomps<-Port[!is.na(Port$Sex) & !is.na(Port$Age),]
nrow(agecomps)
agecomps.rand<-Port.rand[!is.na(Port.rand$Sex) & !is.na(Port.rand$Age),]
nrow(agecomps.rand)

agecomps<-agecomps.rand
#agecomps<-agecomps[agecomps$Project == "Halibut Longline" | agecomps$Project == "Commercial Longline Trip",]
unique(agecomps$Groundfish.Management.Area.Code)
unique(agecomps$GFMU)
#This includes random samples in SEO
agecomps<-agecomps[agecomps$GFMU == "NSEO" |
			agecomps$GFMU == "SSEO" |
			agecomps$GFMU == "CSEO" |
			agecomps$GFMU == "EYKT",]
nrow(agecomps)
hist(agecomps$Age)

ggplot(agecomps, aes(x=Age)) +
  geom_histogram(binwidth=10) +
  facet_wrap(~Year) +
  geom_vline(xintercept=50,col="blue") +
  geom_vline(xintercept=75,col="blue", linetype=2) +
  geom_vline(data=agecomps %>% 
               group_by(Year) %>% 
               summarize(meanage = mean(Age)),
             aes(xintercept = meanage),col="red")

ggplot(agecomps %>% 
         group_by(Year) %>% 
         summarize(meanage = mean(Age)),
       aes(x=meanage)) + geom_histogram()

# mean_age

## get proportion of each age by year...
## will rerun by length bin... 
agecomps<-agecomps %>% 
  group_by(Year,GFMU, Sex) %>%
  mutate(n=n()) %>%
  group_by(Age,Year,GFMU, Sex) %>%
  mutate(count=n()) %>%
  mutate(proportion = count/n)

unique(agecomps$Sex)
EX<-agecomps[agecomps$Year == 2014 & agecomps$GFMU == "CSEO" & agecomps$Sex == "Female",]
view(EX)
plot(EX$Age ~ EX$Sex)
plot(EX$proportion ~ EX$Age)
view(EX[EX$proportion == 1,])
sum(EX$proportion)
hist(EX$proportion)

##make sure proportion sums to 1... 
pr<-unique(EX$Age)
OUT<-data.frame()
i<-1
for (p in pr) {
  SS<-EX[EX$Age == p,]
  OUT[i,"Age"]<-p
  OUT[i,"prop"]<-mean(SS$proportion, na.rm=T)
  i<-i+1
}
sum(OUT$prop, na.rm=T) 

####

str(Port)
lcomps<-Port[!is.na(Port$Sex) & !is.na(Port$Length.Millimeters),]
nrow(lcomps)
lcomps<-Port.rand[!is.na(Port.rand$Sex) & !is.na(Port.rand$Length.Millimeters),]
nrow(lcomps)

lcomps<-lcomps %>% 
  #group_by(Year,Groundfish.Management.Area.Code, Sex) %>%
  group_by(Year,GFMU, Sex) %>%
  mutate(n=n()) %>%
  #group_by(Length.Millimeters,Year,Groundfish.Management.Area.Code, Sex) %>%
  group_by(Length.Millimeters,Year,GFMU, Sex) %>%
  mutate(count=n()) %>%
  mutate(proportion = count/n)
unique(lcomps$Length.Millimeters)

ggplot(lcomps, aes(x=Length.Millimeters)) +
  geom_histogram(binwidth=10, aes(col=Sex, fill = Sex), alpha=0.5) +
  facet_wrap(~Year) +
  geom_vline(aes(xintercept=mean(Length.Millimeters),col=Sex), linetype=2)

lcomps$length<-lcomps$Length.Millimeters
unique(lcomps$length)

ggplot(lcomps, aes(x=length)) +
  geom_histogram(binwidth=10) +
  geom_vline(xintercept=mean(lcomps$length),col="blue") +
  geom_vline(xintercept=quantile(lcomps$length,0.9),col="blue", linetype=2) +
  facet_wrap(~Year, scales="free")

ggplot(lcomps, aes(x=length)) +
  geom_histogram(binwidth=10) +
  geom_vline(xintercept=mean(lcomps$length),col="blue") +
  geom_vline(xintercept=quantile(lcomps$length,0.9),col="blue", linetype=2) +
  facet_wrap(~GFMU, scales="free")

#----------------------------------------------------------------
#Make Length bins

lcomps %>% filter(!c(length < 40)) %>% 
  mutate(length2 = ifelse(length < 41, 41,
                          ifelse(length > 900, 900, length)),
         length_bin = cut(length2, breaks = seq(189.9, 899.9, 10),
                          labels = paste(seq(195, 895, 10)))) %>% 
  select(-length2) -> lendat
str(lendat)
unique(lendat$length)
unique(lcomps$length)
unique(lendat$GFMU)
length(lendat$length[lendat$GFMU == "CSEO"])

lendat.bin<-lendat %>% 
  #group_by(Year,Groundfish.Management.Area.Code, Sex) %>%
  group_by(Year,GFMU, Sex) %>%
  dplyr::mutate(bin.n=n()) %>%
  #group_by(length_bin,Year,Groundfish.Management.Area.Code, Sex) %>%
  group_by(length_bin,Year,GFMU, Sex) %>%
  dplyr::mutate(bin.count=n()) %>%
  mutate(bin.proportion = bin.count/bin.n) 


tst<-lendat %>% filter(GFMU == "CSEO")
unique(tst$length)
unique(tst$Year)
unique(tst$GFMU)

###########################################################################
plus_group<-100
###########################################################################
## age and length compositions for the SEO as a whole


###############################################################################
## Base plots of age and length compositions for management areas: 

unique(Port$GFMU)
gmus <- unique(Port$GFMU)[-c(5, 6)]

for (i in gmus){  
  i<-gmus[1]
  ## AGE BUBBLE PLOTS
  agecompdat <- agecomps %>% 
    filter(Sex %in% c("Female", "Male") &
             GFMU == i &
             Year >= 1990 &
             Age <= plus_group) %>% 
    ungroup()
  
  ggplot(data = agecompdat,
         aes(x = Year, y = Age, size = count)) + #*FLAG* could swap size with proportion_scaled
    geom_point(shape = 21, fill = "black", colour = "black") +
    scale_size(range = c(0, 2)) +
    facet_wrap(~ Sex ) +
    labs(x = "Year", y = "Observed Age") +	#"\nYear", y = "Observed Age\n") +
    guides(size = "none")	+	#FALSE)# +
    theme(text=element_text(size=12),
          axis.text.x = element_text(size=10,angle=45, hjust=1),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=15)) +	#, labels = axisx$labels) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=15)) 

    ggsave(paste0("Figures/Bio_plots_",YEAR,"/",i,"_bubble_agecomp_byyear.png",sep=""),
           dpi=900, height=5, width=7, units="in")
    
    #Length plots
    lendat %>% 
      filter(GFMU == i & Year > 1983) %>% 
      ggplot(aes(length, Year, group = Year, fill = Year)) + 
      geom_density_ridges(aes(point_fill = Year, point_color = Year),
                          alpha = 0.3, scale=1.5, jittered_points = TRUE,
                          rel_min_height = 0.01, point_size=0.5) +
      xlim(300, 900) + 
      xlab("Length (mm)") + 
      ylab(NULL) +
      theme(legend.position = "none") + 
      theme(text=element_text(size=20),
            axis.text.x = element_text(size=10,angle=45, hjust=1),
            axis.text.y = element_text(size=12)) +
      facet_wrap(~ GFMU + Sex) + 
      facet_grid(GFMU ~ Sex)+	
      scale_x_continuous(limits=c(300,900),breaks = c(seq(from=300, to=900, by=50))) +	#, labels = axisx$labels) +
      scale_y_continuous(limits=c(1984,YEAR+2),breaks = c(seq(from=1984, to=(YEAR), by=2)))#scales::pretty_breaks(n=15)) 
    
    ggsave(paste0("Figures/Bio_plots_",YEAR,"/",i,"_ridgejitter_lengthcomp_byyear.png",sep=""),
           dpi=900, height=8, width=5, units="in")
    
    ##SAMPLE SIZE LABEL OPTION 
    lendat %>% 
      filter(GFMU == i & Year > 1983) %>% 
      ggplot(aes(length, Year, group = Year, fill = Year)) + 
      geom_density_ridges(aes(point_fill = Year, point_color = Year),
                          alpha = 0.3, scale=2, #jittered_points = TRUE,
                          rel_min_height = 0.01) +
      xlim(300, 800) + 
      xlab("Length (mm)") + 
      ylab(NULL) +
      theme(legend.position = "none") + 
      theme(text=element_text(size=20),
            axis.text.x = element_text(size=10,angle=45, hjust=1),
            axis.text.y = element_text(size=12)) +
      facet_wrap(~ GFMU + Sex) + 
      facet_grid(GFMU ~ Sex)+	#+
      geom_label(data = lendat%>% 
                   filter(GFMU == i) %>% group_by(Year), 
                 aes(y = Year, x = 775, 
                     label = n), colour="black", fill="white", nudge_y=0.5, size=3)+
      scale_x_continuous(limits=c(300,800),breaks = c(seq(from=300, to=800, by=50))) +	#, labels = axisx$labels) +
      scale_y_continuous(breaks = c(seq(from=1984, to=YEAR, by=2)) )
    
    ggsave(paste0("Figures/Bio_plots_",YEAR,"/",i,"_ridgelabeled_lengthcomp_byyear.png",sep=""),
           dpi=900, height=8, width=5, units="in")
}

#SEO composite plots
str(agecomps)
unique(agecomps$Sample.Type)

agecomps %>% group_by(Year, Sex) %>%
  summarize(samples = n(),
            landings = length(unique(Trip.Number))) -> age_smpl_size

str(age_smpl_size)

ggplot(age_smpl_size,aes(Year,samples)) +
  geom_col(aes(col=Sex, fill=Sex),alpha=0.5, position = position_dodge(width = 0.2)) +
  ylab("Number of age samples") +
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position.inside = c(0.2,0.7)) +
  scale_x_continuous(breaks=seq(1984,2024,2)) ->raw_ss_plot;raw_ss_plot

ggplot(agecomps %>% group_by(Year) %>%
         summarize(samples = n(),
                   landings = length(unique(Trip.Number))),
       aes(Year,landings)) +
  geom_col(aes(),alpha=0.5, position = position_dodge(width = 0.2)) +
  ylab("Number of landings sampled") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_continuous(breaks=seq(1984,2024,2))-> lnds_ss_plot

ggarrange(raw_ss_plot, lnds_ss_plot, nrow=2)

ps<-align_plots(raw_ss_plot,lnds_ss_plot, align = "hv")

ggdraw() + draw_grob(ps[[1]], 0,0.45,1,0.5) +   #grob, x, y, width, height, scale, clip, hjust, vjust, halign, valign
  draw_grob(ps[[2]], 0,0,1,0.5)

ggsave(paste0("Figures/age_ss_byyear_", YEAR, ".png"), dpi=300,  height=5, width=7, units="in")
  
  
ggplot(data = agecomps,
       aes(x = Year, y = Age, size = count)) + #*FLAG* could swap size with proportion_scaled
  geom_point(shape = 21, fill = "black", colour = "black") +
  scale_size(range = c(0, 1)) +
  facet_wrap(~ Sex ) +
  labs(x = "Year", y = "Observed Age") +	#"\nYear", y = "Observed Age\n") +
  guides(size = "none")	+	#FALSE)# +
  theme(text=element_text(size=12),
        axis.text.x = element_text(size=10,angle=45, hjust=1),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=15)) +	#, labels = axisx$labels) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=15)) 

ggsave(paste0("Figures/Bio_plots_",YEAR,"/","SEO","_bubble_agecomp_byyear.png",sep=""),
       dpi=900, height=5, width=7, units="in")

#################### Age bubble plots ########################################
# Code from Adam St. Saviour
# Bubble plot of SF annual AGE distribution
# Count the frequency of each AGE per YEAR
age_distribution <- agecomps %>%
  group_by(Year, Age) %>%
  summarize(count = n())

# bubble plot of AGE distribution
bubble_plot <- ggplot(age_distribution, aes(x = Year, y = Age, size = count)) +
  geom_point(alpha = 0.3, color = 'black') +
  scale_size_continuous(range =  c(1, 20), breaks = c(50,100)) +  # Adjust the size range and Count breaks in legend
  theme_minimal() +
  labs(title = "Yelloweye Rockfish Age Distribution of Commercial Harvest in SEO",
       x = "Year",
       y = "Age",
       size = "Count") +
  theme(plot.title = element_text(hjust = 0.5)); bubble_plot

ggsave(paste0("Figures/Bio_plots_",YEAR,"/","SEO","_bubble_agecomp_byyear2.png",sep=""),
       dpi=900, height=5, width=7, units="in")

age_distribution_sex <- agecomps %>%
  group_by(Year, Age, Sex) %>%
  summarize(count = n())

bubble_plot_bysex <- ggplot(age_distribution_sex, aes(x = Year, y = Age, size = count)) +
  geom_point(alpha = 0.3, color = 'black') +
  scale_size_continuous(range =  c(1, 20), breaks = c(20,40)) +  # Adjust the size range and Count breaks in legend
  facet_grid(~Sex)+
  theme_minimal() +
  labs(title = "Yelloweye Rockfish Age Distribution of Commercial Harvest in SEO by Sex",
       x = "Year",
       y = "Age",
       size = "Count") +
  theme(plot.title = element_text(hjust = 0.5)); bubble_plot_bysex


ggsave(paste0("Figures/Bio_plots_",YEAR,"/","SEO","_bubble_agecomp_byyear2_bysex.png",sep=""),
       dpi=900, height=5, width=7, units="in")


#Length plots - random samples from SEAK - including inside waters
lendat %>% 
  filter(Year > 1983) %>% 
  ggplot(aes(length, Year, group = Year, fill = Year)) + 
  geom_density_ridges(aes(point_fill = Year, point_color = Year),
                      alpha = 0.3, scale=1.5, jittered_points = FALSE,
                      rel_min_height = 0.01, point_size=0.5) +
  xlim(300, 900) + 
  xlab("Length (mm)") + 
  ylab(NULL) +
  theme(legend.position = "none") + 
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=10,angle=45, hjust=1),
        axis.text.y = element_text(size=12)) +
  facet_wrap(~  Sex) + 
  facet_grid( ~ Sex)+	
  scale_x_continuous(limits=c(300,900),breaks = c(seq(from=300, to=900, by=50))) +	#, labels = axisx$labels) +
  scale_y_continuous(limits=c(1984,YEAR+2),breaks = c(seq(from=1984, to=(YEAR), by=2)))#scales::pretty_breaks(n=15)) 

ggsave(paste0("Figures/Bio_plots_",YEAR,"/","SEAK","_ridgejitter_lengthcomp_byyear.png",sep=""),
       dpi=900, height=8, width=5, units="in")

#Same plot as above but filter out NSEI and SSEI

lendat %>% 
  filter(Year > 1983) %>% 
  filter(!Groundfish.Management.Area.Code %in% c("SSEI","NSEI")) %>% 
  ggplot(aes(length, Year, group = Year, fill = Year)) + 
  geom_density_ridges(aes(point_fill = Year, point_color = Year),
                      alpha = 0.3, scale=1.5, jittered_points = FALSE,
                      rel_min_height = 0.01, point_size=0.5) +
  xlim(300, 900) + 
  xlab("Length (mm)") + 
  ylab(NULL) +
  theme(legend.position = "none") + 
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=10,angle=45, hjust=1),
        axis.text.y = element_text(size=12)) +
  facet_wrap(~  Sex) + 
  facet_grid( ~ Sex)+	
  scale_x_continuous(limits=c(300,900),breaks = c(seq(from=300, to=900, by=50))) +	#, labels = axisx$labels) +
  scale_y_continuous(limits=c(1984,YEAR+2),breaks = c(seq(from=1984, to=(YEAR), by=2)))#scales::pretty_breaks(n=15)) 

ggsave(paste0("Figures/Bio_plots_",YEAR,"/","SEO","_ridgejitter_lengthcomp_byyear.png",sep=""),
       dpi=900, height=8, width=5, units="in")

##SAMPLE SIZE LABEL OPTION 
lendat %>% 
  filter(Year > 1983) %>% 
  ggplot(aes(length, Year, group = Year, fill = Year)) + 
  geom_density_ridges(aes(point_fill = Year, point_color = Year),
                      alpha = 0.3, scale=2, #jittered_points = TRUE,
                      rel_min_height = 0.01) +
  xlim(300, 800) + 
  xlab("Length (mm)") + 
  ylab(NULL) +
  theme(legend.position = "none") + 
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=10,angle=45, hjust=1),
        axis.text.y = element_text(size=12)) +
  facet_wrap(~  Sex) + 
  facet_grid( ~ Sex)+	#+
  geom_label(data = lendat %>% 
               group_by(Year), 
             aes(y = Year, x = 775, 
                 label = n), colour="black", fill="white", nudge_y=0.5, size=3)+
  scale_x_continuous(limits=c(300,800),breaks = c(seq(from=300, to=800, by=50))) +	#, labels = axisx$labels) +
  scale_y_continuous(breaks = c(seq(from=1984, to=YEAR, by=2)) )

ggsave(paste0("Figures/Bio_plots_",YEAR,"/","SEO","_ridgelabeled_lengthcomp_byyear.png",sep=""),
       dpi=900, height=8, width=5, units="in")

############################################################################
# Biological data exam for ASA 
colnames(Port.rand)

Bio<-Port.rand %>% 
  mutate(Year = factor(Year),
         year = as.numeric(as.character(Year)),
        Project_cde = factor(Project.Code),
        Stat = factor(Groundfish.Stat.Area),
                    # Station = factor(Station),
        Sex = factor(Sex),
        Mature = Maturity,
        length = Length.Millimeters,
        weight = Weight.Kilograms,
        age=Age) %>% 
  group_by(Year, Stat) %>% 
  mutate(n = length(Age),
         length_mu = mean(length, na.rm = TRUE),
         weight_mu = mean(weight, na.rm = TRUE)) %>% 
  ungroup()

str(Bio)
#------------------------------------------------------------------------------
# male-female length distributions
ggplot(Bio) +
  geom_histogram(aes(x=length, col=Sex, fill=Sex))

#------------------------------------------------------------------------------
#Length-weight relationship

Bio_lw<-Bio %>% filter(Sex %in% c("Female", "Male") &
                         !is.na(length) &
                         !is.na(weight))

lw_allometry <- function(length, a, b) {a * length ^ b}

START <- c(a = 1e-5, b = 3) #Starting values close to Hanselman et al. 2007

fem_fit <- nls(weight ~ lw_allometry(length = length, a, b), 
               data = filter(Bio_lw, Sex == "Female"), start = START)

male_fit <- nls(weight ~ lw_allometry(length = length, a, b), 
                data = filter(Bio_lw, Sex == "Male"), start = START)

all_fit <- nls(weight ~ lw_allometry(length = length, a, b), 
               data = Bio_lw, start = START)

# parameter estimates and plot fit 
beta_m <- tidy(male_fit)$estimate[2]
beta_f <- tidy(fem_fit)$estimate[2]
beta_a <- tidy(all_fit)$estimate[2]

Bio_lw<-bind_rows(Bio_lw %>% filter(Sex %in% "Male") %>% 
                         mutate(Condition = weight/(length^beta_m),
                                Quantile = ntile(Condition,1000)/1000),
                  Bio_lw %>% filter(Sex %in% "Female") %>% 
                         mutate(Condition = weight/(length^beta_f),
                                Quantile = ntile(Condition,1000)/1000))

#filtering out the outliers
Bio_lw_noOL<-Bio_lw[Bio_lw$Quantile > 0.0005 & Bio_lw$Quantile < 0.9995, ]
nrow(Bio_lw_noOL)
nrow(Bio_lw)

# Identify the outliers
outliers <- Bio_lw %>% filter(!(Quantile > 0.0005 & Quantile < 0.9995))

# Write the outliers to a CSV file
write.csv(outliers, "Output/Bio_lw_outliers.csv", row.names = FALSE)

hist(Bio_lw$length, breaks = 50)
max(Bio_lw$length, na.rm=T); 
plot(Bio_lw$weight ~ Bio_lw$length)
points(Bio_lw_noOL$weight ~ Bio_lw_noOL$length, col="blue")

#---------------------------------------------------------------------------
# Empirical weight-at-age ---- outliers left in for this piece
unique(Bio_lw_noOL$Project.Code)
colnames(Bio_lw_noOL)
unique(Bio_lw$Project)
with(Bio_lw, table(Project, Project.Code))
# All years combined
Bio_lw_noOL %>% 
  filter(!is.na(weight) & !is.na(age) & !is.na(Sex)) %>% 
  mutate(Gear = derivedFactor('Longline' = Project.Code %in% c("602","608","621"),
                                'Jig' = Project.Code %in% c("604")),
         Fishery = Project) %>% 
  filter(Age <= plus_group) -> waa

ggplot() +
  geom_point(data = waa, # %>% filter(!is.na(Source)), 
             aes(x = age, y = weight, col = Gear, shape = Sex))

bind_rows(
  waa %>%
    group_by(Gear, Sex, age) %>% 
    dplyr::summarise(weight = mean(weight) %>% round(4)),
  waa %>% 
    group_by(Gear, age) %>% 
    dplyr::summarise(weight = mean(weight) %>% round(4)) %>% 
    dplyr::mutate(Sex = "Combined")) -> emp_waa

ggplot() +
  geom_point(data = emp_waa, # %>% filter(!is.na(Source)), 
             aes(x = age, y = weight, col = Gear, shape = Sex))


str(emp_waa)

write_csv(emp_waa, paste0("Output/empircal_waa_", YEAR, ".csv"))

# Expand to grid to include all age combos and fill in NAs if there are any
# using linear interpolation
rec_age<-min(emp_waa$age)

expand.grid(Gear = unique(emp_waa$Gear),
            Sex = unique(emp_waa$Sex),
            age = seq(rec_age, plus_group, 1))  %>% 
  data.frame()  %>% 
  full_join(emp_waa) %>%
  group_by(Gear, Sex) %>% 
  mutate(weight = zoo::na.approx(weight, maxgap = 20, rule = 4)) -> emp_waa_interpol

ggplot() +
  geom_point(data = emp_waa, # %>% filter(!is.na(Source)), 
             aes(x = age, y = weight, col = Gear, shape = Sex))

write_csv(emp_waa, paste0("Output/empircal_waa_", YEAR, ".csv"))

# Changes in weight-at-age

Bio_lw_noOL %>% 
  select(Year, Project.Code, Sex, age, weight) %>% 
  #filter(year >= 1997) %>% 
  filter(!is.na(weight) & !is.na(age) & !is.na(Sex)) %>% 
  #filter(as.numeric(year) <= YEAR & age >= rec_age) %>% 
  group_by(Year, Sex, age) %>% 
  dplyr::summarise(weight = mean(weight) %>% round(4)) %>%  
  ungroup() %>% 
  mutate(Year = as.character(Year),
         Age = factor(age),
         cohort = as.numeric(Year) - age,
         Cohort = as.factor(cohort))-> df

unique(df$Year)

pal <- ggthemes::canva_pal("Warm and cool")(4) 

# By cohort
df_cohort <- df %>% 
  #filter(cohort >= 2010 & cohort <= YEAR-3 & age >=2 & age <= 5) %>% 
  droplevels()

# Axis ticks for plot (see helper.r tickr() fxn for details)
# axis <- tickr(df_cohort, year, 2)

ggplot(df_cohort, aes(age, weight, colour = Cohort, group = Cohort)) +
  #geom_line(linewidth = 1) +
  geom_point(aes(fill = Cohort), show.legend = FALSE, size = 1) +
  facet_grid(~ Sex) +
  labs(x = "Age", y = "Weight-at-age (grams)\n", colour = "Cohort") +
  guides(colour = guide_legend(ncol = 9)) +
  scale_colour_manual(values = colorRampPalette(pal)(n_distinct(df_cohort$Cohort))) +
  theme(legend.position = "bottom") #+
# scale_x_continuous(breaks = axis$breaks, labels = axis$labels)# -> waa_cohort_plot

ggsave(paste0("Figures/waa_cohort_",YEAR,".png"), dpi = 300, height = 5, width = 7, units = "in")

#df %>% 
#  filter(Age %in% c("2", "3", "4", "5")) %>% 
#  droplevels() -> df
df %>% 
  group_by(Age, Sex) %>% 
  dplyr::summarize(mean_weight = mean(weight, na.rm = TRUE)) -> means

# axis <- tickr(df, year, 5)

ggplot(df, 
       aes(Year, weight, group = Age, colour = Age)) + 
  geom_line() + 
  geom_point() +
  facet_wrap(~ Sex, ncol = 1) +
  geom_hline(data = means, aes(colour = Age, yintercept = mean_weight), alpha = 0.4, linetype = 2) + 
  labs(x = "Year", y = "Weight-at-age (grams)\n", colour = "Age") +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = colorRampPalette(pal)(n_distinct(df$Age))) +
  guides(colour = guide_legend(nrow = 1))# +
# scale_x_continuous(breaks = axis$breaks, labels = axis$labels)

ggsave(paste0("Figures/waa_trends.png"), dpi = 300, height = 5, width = 7, units = "in")


#=========================================================================================
# Length-at-age
# subsets by length, age, sex
Bio_lw_noOL %>% 
  filter(Sex %in% c("Female", "Male") &
           #year >= 1997 & # *FLAG* advent of "modern" survey
           !is.na(length) &
           !is.na(age)) %>% 
  droplevels() -> laa_sub

laa_sub %>% 
  ungroup() %>% 
  filter(Sex == "Female") -> laa_f

laa_sub %>% 
  ungroup() %>% 
  filter(Sex == "Male") -> laa_m

#find some starting values
plot(laa_f$length ~ laa_f$age)
age<-sort(unique(laa_f$age))
Linf<-650
k<-0.05
t0<--5
length<-Linf*(1-exp(-k*(age-t0)))
lines(length~age)

# females
lvb_f <- fishmethods::growth(unit = 1, # length (2 = weight)
                             size = laa_f$length, age = laa_f$age,
                             error = 1, # additive (2 = multiplicative)
                             # starting values from Hanselman et al. 2007 (Appendix C, Table 1)
                             Sinf = Linf, K = k, t0 = t0)
lvb_f

# males
lvb_m <- fishmethods::growth(unit = 1, # length (2 = weight)
                             size = laa_m$length, age = laa_m$age,
                             error = 1, # additive (2 = multiplicative)
                             # starting values from Hanselman et al. 2007 (Appendix C, Table 1)
                             Sinf = Linf, K = k, t0 = t0)
lvb_m

# save param estimates
as_tibble(summary(lvb_f$vout)$parameters[,1:2]) %>% 
  dplyr::mutate(Parameter = c("l_inf", "k", "t0"),
         Sex = "Female") %>% 
  bind_rows(as_tibble(summary(lvb_m$vout)$parameters[,1:2]) %>% 
              mutate(Parameter = c("l_inf", "k", "t0"),
                     Sex = "Male")) %>% mutate(Source = "Commercial fisheries",
                                               Years = paste0(min(laa_sub$year), "-", max(laa_sub$year)),
                                               Region = "SEO",
                                               Function = "Length-based LVB") %>% 
  full_join(laa_sub %>% 
              group_by(Sex) %>% 
              dplyr::summarise(n = n())) -> lvb_pars

# Are there annual trends in length-at-age? First vonB curves by Year. I only
# look at Linf here because vonB parameters are highly correlated. The following
# code could be adapted to look at trends in K. There's also a nifty fxn in this
# package growthlrt() that does likelihood ratio tests to compare multiple
# growth curves
laa_f <- laa_f %>% 
  mutate(iyr = as.integer(as.factor(Year)),
         gyear = paste0(LETTERS[iyr], Year)) %>% #gyears? add sequential letters to 
  # front of year bc growthmultifit needs letter at front for group
  arrange(gyear)

laa_f$year
laa_f$gyear

#l_inf for females in each year... 
multi_f <- fishmethods::growthmultifit(len = laa_f$length, age = laa_f$age, 
                                       group = as.character(laa_f$gyear),
                                       model = 1,
                                       fixed = c(2, 1, 1), # 2 = not fixed, 1 = fixed. order is Linf, K, t0
                                       select = 2,
                                       error = 1, # additive (2 = multiplicative)
                                       # starting values from Hanselman et al. 2007 (Appendix C, Table 1)
                                       Linf = rep(Linf, length(unique(laa_f$iyr))), 
                                       K = rep(k, length(unique(laa_f$iyr))), 
                                       t0 = rep(t0, length(unique(laa_f$iyr))),
                                       plot = FALSE)

laa_m <- laa_m %>% 
  mutate(iyr = as.integer(as.factor(Year)),
         gyear = paste0(LETTERS[iyr], Year)) %>% 
  arrange(gyear)

multi_m <- growthmultifit(len = laa_m$length, age = laa_m$age, 
                          group = as.character(laa_m$gyear),
                          model = 1,
                          fixed = c(2, 1, 1), 
                          select = 2,
                          error = 1, # additive (2 = multiplicative)
                          # starting values from Hanselman et al. 2007 (Appendix C, Table 1)
                          Linf = rep(Linf, length(unique(laa_m$iyr))), 
                          K = rep(k, length(unique(laa_m$iyr))), 
                          t0 = rep(t0, length(unique(laa_m$iyr))),
                          plot = FALSE)

multi <- as_tibble(multi_f$results$parameters[,1:2]) %>% 
  mutate(Parameter = rownames(multi_f$results$parameters),
         Sex = "F") %>% 
  bind_rows(as_tibble(multi_m$results$parameters[,1:2]) %>% 
              mutate(Parameter = rownames(multi_m$results$parameters),
                     Sex = "M")) %>% 
  select(Sex, Parameter, Estimate)

view(multi) #year specific linf for male and female, also one row each for male/fem t0 and K

#get linf for each year sex by adding year specific deviation to year 1
#add col for sex specifc t0 and k, which not varying by year
multi <- multi %>%   #getting linf for each year with a bunch of dplyr steps...
  filter(Parameter %in% c("Linf1", "K1", "t01")) %>% 
  pivot_wider(id_cols = c(Sex), names_from = Parameter, values_from = Estimate) %>% 
  left_join(multi %>% 
              filter(!Parameter %in% c("Linf1", "K1", "t01"))) %>% 
  mutate(Estimate = Linf1 + Estimate) %>% 
  pivot_wider(names_from = Parameter, values_from = Estimate) %>% 
  pivot_longer(cols = -c(K1, t01, Sex), names_to = "Parameter", values_to = "Estimate") %>% 
  bind_cols(as_tibble(c(unique(laa_f$Year), unique(laa_m$Year)))) %>% 
  dplyr::rename(Year = value)

# Deviations by year for Linf vonB 
#dev.off() if not plotting may need to run this...
ggplot() + 
  geom_segment(data = multi %>% 
                 group_by(Sex) %>%
                 mutate(scaled_est = scale(Estimate),
                        mycol = ifelse(scaled_est < 0, "blue", "red")) %>% 
                 ungroup(),
               aes(x = Year, y = 0,
                   xend = Year, yend = scaled_est, 
                   color = mycol), size = 2) +
  scale_colour_grey() +
  geom_hline(yintercept = 0, lty = 2) + 
  guides(colour = FALSE) +
  labs(x = "", y = "Scaled parameter estimates of Linf\n") +
  facet_wrap(~ Sex, ncol = 1)

ggsave(paste0("Figures/trends_Linf_1997_", YEAR, ".png"), 
       dpi=300, height=7, width=6, units="in")

#Save length-at-age predictions

age_pred <- data.frame(age = rec_age:plus_group)

# Survey laa
laa_preds <- age_pred %>% 
  mutate(sex = "Female",
         length = predict(lvb_f$vout, newdata = age_pred)) %>% 
  bind_rows(age_pred %>% 
              mutate(sex = "Male",
                     length = predict(lvb_m$vout, newdata = age_pred))) %>% 
  mutate(Source = "commercial fisheries") 

write_csv(laa_preds, paste0("Output/pred_laa_plsgrp_", plus_group, "_", YEAR, ".csv"))

#--------------------------------------------------------------------------
# Weight-length allometry: W = alpha*L^beta

str(Bio_lw)

Bio_lw_noOL %>% 
  filter(Sex %in% c("Female", "Male") &
           #year >= 1997 & #advent of "modern" survey
           !is.na(length) &
           !is.na(weight)) %>% 
  droplevels() -> allom_sub

# length-weight relationship
lw_allometry <- function(length, a, b) {a * length ^ b}

START <- c(a = 1e-5, b = 3) #Starting values close to Hanselman et al. 2007

fem_fit <- nls(weight ~ lw_allometry(length = length, a, b), 
               data = filter(allom_sub, Sex == "Female"), start = START)

male_fit <- nls(weight ~ lw_allometry(length = length, a, b), 
                data = filter(allom_sub, Sex == "Male"), start = START)

all_fit <- nls(weight ~ lw_allometry(length = length, a, b), 
               data = allom_sub, start = START)

# parameter estimates and plot fit 
beta_m <- tidy(male_fit)$estimate[2]
beta_f <- tidy(fem_fit)$estimate[2]
beta_a <- tidy(all_fit)$estimate[2]

unique(allom_sub$Sex)

bind_rows(tidy(male_fit) %>% mutate(Sex = "Male"),     
          tidy(fem_fit) %>% mutate(Sex = "Female")) %>% 
  bind_rows(tidy(all_fit) %>% mutate(Sex = "Combined")) %>% 
  dplyr::select(Parameter = term, Estimate = estimate, SE = std.error, Sex) %>% 
  mutate(Source = "Commercial fisheries",
         Years = paste0(min(as.numeric(laa_sub$Year)), "-", max(as.numeric(laa_sub$Year))),
         Region = "SEO",
         Function = "Allometric - NLS") %>% 
  full_join(allom_sub %>% 
              group_by(Sex) %>% 
              dplyr::summarise(n = n()) %>% 
              bind_rows(allom_sub %>% 
                          dplyr::summarise(n = n()) %>% 
                          mutate(Sex = "Combined"))) -> allom_pars

ggplot(allom_sub, aes(length, weight, col = Sex, shape = Sex)) +
  geom_jitter(alpha =0.8) + 
  stat_function(fun = lw_allometry, 
                args = as.list(tidy(fem_fit)$estimate),
                col = "#F8766D") + 
  stat_function(fun = lw_allometry, 
                args = as.list(tidy(male_fit)$estimate),
                col = "#00BFC4", lty = 2) + 
  labs(x = "Fork length (cm)", y = "Round weight (kg)", alpha = NULL) +
  # scale_colour_grey() +
  theme(legend.position = c(0.85, 0.2))

ggsave(paste0("Figures/allometry_commfisheries_",min(as.numeric(laa_sub$Year)),"-", YEAR, ".png"),
       dpi=300, height=4, width=6, units="in")

#--------------------------------------------------------------------------
# Weight-at-age:
# subsets by weight, age, sex
Bio_lw_noOL %>% 
  filter(Sex %in% c("Female", "Male") &
           #year >= 1997 & # *FLAG* advent of "modern" survey
           !is.na(age) &
           !is.na(weight)) %>% 
  droplevels() -> waa_sub

waa_sub %>% 
  ungroup() %>% 
  filter(Sex == "Female") -> waa_f

waa_sub %>% 
  ungroup() %>% 
  filter(Sex == "Male") -> waa_m

#find some starting values
plot(waa_f$weight ~ waa_f$age)
age<-sort(unique(waa_f$age))
Winf<-6
k_w<-0.025
t0_w<--2
length<-Winf*(1-exp(-k_w*(age-t0_w)))
lines(length~age, col="red")

#w at age curves
# females
wvb_f <- fishmethods::growth(unit = 2, # 1 = length, 2 = weight
                             size = waa_f$weight, age = waa_f$age,
                             error = 2, # 1 = additive, 2 = multiplicative log(w_i) = log(w_inf) + beta * log(1 - exp * (-k * (age_i - t0))) + error
                             Sinf = Winf, K = k_w, t0 = t0_w,
                             B = allom_pars %>% filter(Sex == "Female" & Parameter == "b") %>% pull(Estimate))
wvb_f # gompertz failed but that wasn't our target so we're good

# males
wvb_m <- fishmethods::growth(unit = 2, 
                             size = waa_m$weight, age = waa_m$age,
                             error = 2, 
                             # starting values from Hanselman et al. 2007 (Appendix C, Table 1)
                             Sinf = Winf, K = k_w, t0 = t0_w,
                             B = allom_pars %>% filter(Sex == "Male" & Parameter == "b") %>% pull(Estimate))
wvb_m

# save param estimates
lvb_pars %>% bind_rows(
  as_tibble(summary(wvb_f$vout)$parameters[,1:2]) %>% 
    mutate(Parameter = c("w_inf", "k", "t0"),
           Sex = "Female") %>% 
    bind_rows(as_tibble(summary(wvb_m$vout)$parameters[,1:2]) %>% 
                mutate(Parameter = c("w_inf", "k", "t0"),
                       Sex = "Male")) %>% mutate(Source = "ADFG longline survey",
                                                 Years = paste0(min(waa_sub$year), "-", max(waa_sub$year)),
                                                 Region = "Chatham Strait",
                                                 Function = "Weight-based LVB") %>% 
    full_join(waa_sub %>% 
                group_by(Sex) %>% 
                summarise(n = n())))  -> lvb_pars

# sexes combined
wvb <- fishmethods::growth(unit = 2, 
                           size = waa_sub$weight, age = waa_sub$age,
                           error = 2, 
                           Sinf = Winf, K = k_w, t0 = t0_w,
                           B = allom_pars %>% filter(Sex == "Combined" & Parameter == "b") %>% pull(Estimate))
wvb # horrible fit. only going through this exercise to supply data inputs to a single sex model if ever desired for comparison

lvb_pars %>% bind_rows(
  as_tibble(summary(wvb$vout)$parameters[,1:2]) %>% 
    mutate(Parameter = c("w_inf", "k", "t0"),
           Sex = "Combined") %>% 
    full_join(waa_sub %>% 
                dplyr::summarise(n = n()) %>% 
                mutate(Sex = "Combined")) %>% 
    mutate(Source = "Commercial fisheries",
           Years = paste0(min(waa_sub$year), "-", max(waa_sub$year)),
           Region = "SEO",
           Function = "Weight-based LVB")) -> lvb_pars


#--------------------------------------------------------------------------
# unfortunately, while the fishmethods package makes it easy to do everything
# else, it does not make it easy to extract model predictions from the growth
# curves when using the multiplicative error structure (the one that assumes
# residuals are lognormally distributed). As such, we have to extract the
# coefficients and do the predictions manually. Part of this is applying a
# log-normal bias correction, which gets us the mean (instead of median)
# predicted value. This correction is E(w_i) = exp(mu_i + 0.5 * sigma^2), where
# mu is the predicted value in log space (i.e. log(w_i) in the formula below):

# log(w_i) = log(w_inf) + beta * log(1 - exp * (-k * (age_i - t0))) + error

coef_f_wvb <- summary(wvb_f$vout)$parameter[,1] #  females
beta_f <- allom_pars %>% filter(Sex == "Female" & Parameter == "b") %>% pull(Estimate) # beta from allometric model (weight-length relationship)

coef_m_wvb <- summary(wvb_m$vout)$parameter[,1] # males
beta_m <- allom_pars %>% filter(Sex == "Male" & Parameter == "b") %>% pull(Estimate) # beta from allometric model (weight-length relationship)

coef_wvb <- summary(wvb$vout)$parameter[,1] # sexes combined
beta_c <- allom_pars %>% filter(Sex == "Combined" & Parameter == "b") %>% pull(Estimate) # beta from allometric model (weight-length relationship)

waa_preds <- age_pred %>% 
  mutate(Sex = "Female",
         round_kg = exp(log(coef_f_wvb[1]) + beta_f * log(1 - exp(-coef_f_wvb[2] * (age_pred$age - coef_f_wvb[3]))) +
                          0.5 * (sigma(wvb_f$vout)^2))) %>% 
  bind_rows(age_pred %>% 
              mutate(Sex = "Male",
                     round_kg = exp(log(coef_m_wvb[1]) + beta_m * log(1 - exp(-coef_m_wvb[2] * (age_pred$age - coef_m_wvb[3]))) +
                                      0.5 * (sigma(wvb_m$vout)^2)))) %>% 
  bind_rows(age_pred %>% 
              mutate(Sex = "Combined",
                     round_kg = exp(log(coef_wvb[1]) + beta_c * log(1 - exp(-coef_wvb[2] * (age_pred$age - coef_wvb[3]))) +
                                      0.5 * (sigma(wvb$vout)^2)))) %>% 
  mutate(Source = "comm fisheries") 
  

write_csv(waa_preds, paste0("Output/pred_waa_plsgrp_", plus_group, "_", YEAR, ".csv"))

# Compare empirical and predicted weight-at-age
ggplot() +
  geom_point(data = emp_waa, # %>% filter(!is.na(Source)), 
             aes(x = age, y = weight, col = Gear, shape = Sex)) +
  geom_line(data = waa_preds,
            aes(x = age, y = round_kg, col = Source, linetype = Sex), size = 1) +
  scale_colour_grey() +
  # expand_limits(y = 0) +
  labs(x = "\nAge", y = "Weight (kg)\n", linetype = "Sex", shape = "Sex")

# ggsave(paste0(YEAR+1,"/figures/compare_empirical_predicted_waa_", YEAR, ".png"), 
#        dpi=300, height=4, width=6, units="in")

#--------------------------------------------------------------------------------
# Maturity: this is based on sablefish model for starters... 

# 0 = immature, 1 = mature. Only conducted for females (because biological
# reference points / harvest strategy are based on female spawning biomass)
Bio_lw_noOL %>% 
  ungroup() %>% 
  filter(Sex == "Female" &
           #year >= 1997 & # *FLAG* advent of "modern" survey
           !is.na(Mature) &
           Mature != c("","Not observed"),
           !is.na(length)) %>% 
  mutate(mat.code = ifelse(Mature %in% c("Immature"),0,
                           ifelse(Mature %in% c("Maturing"),0.5,1))) %>%
  droplevels() -> len_f; nrow(len_f)

Bio_lw_noOL %>% 
  ungroup() %>% 
  filter(Sex == "Male" &
           #year >= 1997 & # *FLAG* advent of "modern" survey
           !is.na(Mature) &
           Mature != c("","Not observed"),
           !is.na(length)) %>% 
  mutate(mat.code = ifelse(Mature %in% c("Immature"),0,
                           ifelse(Mature %in% c("Maturing"),0.5,1))) %>%
  droplevels() -> len_m; nrow(len_m)

# Sample sizes by year
with(len_f, table(year,Mature))
with(len_f, table(year,Maturity))
with(len_f, table(year,mat.code))
unique(len_f$Maturity)

ggplot(len_f) +
  geom_point(aes(x=length, y = mat.code))

# base model
fit_length <- glm(mat.code ~ length, data = len_f %>% filter(mat.code != 0.5),
                  family = binomial)
len <- seq(0, max(len_f$length), 0.05)
(L50 <- round(- coef(fit_length)[1]/coef(fit_length)[2],1))
(kmat <- round(((coef(fit_length)[1] + coef(fit_length)[2]*len) / (len - L50))[1], 2))

plot(data = len_f, mat.code ~ length)
plot(fit_length)
# by year, for comparison
fit_length_year <- glm(mat.code ~ length * Year, 
                       data = len_f %>% filter(mat.code != 0.5), 
                       family = binomial)
plot(fit_length_year)

summary(fit_length_year)
AIC(fit_length, fit_length_year)

# fit_length_year <- glm(Mature ~ length * Year, data = len_f, family = quasibinomial)

# New df for prediction for fit_length
new_len_f_simple <- data.frame(length = seq(0, max(len_f$length), 0.05))

# New df for prediction for fit_length_year
new_len_f <- data.frame(length = rep(seq(0, max(len_f$length), 0.05), n_distinct(len_f$year)),
                        Year = factor(sort(rep(unique(len_f$year), length(seq(0, max(len_f$length), 0.05))), decreasing = FALSE)))

# Get predicted values for fit_length   ; ?broom::augment    
broom::augment(x = fit_length, 
               newdata = new_len_f_simple, 
               type.predict = "response") %>% 
  select(length, fitted = .fitted) -> pred_simple #, se =.se.fit #(earlier versions of broom output .se.fit!)

# Get predicted values by year for fit_length_year          
broom::augment(x = fit_length_year, 
               newdata = new_len_f, 
               type.predict = "response") %>% 
  select(Year, length, fitted = .fitted) -> pred #, se =.se.fit #(earlier versions of broom output .se.fit!)

#Length-based maturity curves - 2019 and 2020 appear to have early
#age-at-maturation relative to other years
ggplot() +
  geom_line(data = pred, 
            aes(x = length, y = fitted, group = Year, colour = as.numeric(as.character(Year)))) +
  geom_line(data = pred_simple, aes(x = length, y = fitted, lty = "All years combined"),
            colour = "black", size = 1) +
  lims(x = c(0, max(len_f$length))) +
  scale_colour_gradientn(colours=rainbow(4)) +
  scale_linetype_manual(values = 2) +
  labs(x = "\nLength (cm)", y = "Probability\n", colour = "Year", lty = NULL) +
  theme(legend.position = c(.8, .4))

ggsave(paste0("Figures/maturity_atlength_byyear_srvfem_", YEAR, ".png"), 
       dpi=300, height=4, width=6, units="in")

# Parameter estimates by year
broom::tidy(fit_length_year) %>% 
  select(param = term,
         est = estimate) -> mature_results
view(mature_results)

#****** FLAG: Of all the samples, only 6 have been observed as immature.
#*#  This data is not suitable for estimating a maturity curve!!! 
#*
#*The rest of this maturity code is copied from nsei_sablefish and is NOT modified 
#* to deal with this data.  Leaving it in here for potential future analysis IF
#* data becomes available...   

{# Note on glm() output b/c I can never remember
# (Intercept) = intercept for Year1997 (or first level of factor) on logit scale
# length = slope for Year1997 on logit scale
# Year1998 = difference in intercept (Year1998 - Year1997)
# length:Year1998 = difference in slope (lengthYear1998 - Year1997)
bind_rows(
  # filter out intercepts, derive estimates by yr  
  mature_results %>% 
    filter(param == "(Intercept)") %>% 
    mutate(param = "Year1997",
           Parameter = "b_0"),
  mature_results %>% 
    filter(!param %in% c("(Intercept)", "length") &
             !grepl(':+', param)) %>%     # filter(!grepl('length:Year\\d{4}', param)) %>% #alternative regex
    mutate(est = est + mature_results$est[mature_results$param == "(Intercept)"],
           Parameter = "b_0"),
  # filter out slopes (contain >1 :), derive estimates by year
  mature_results %>% 
    filter(param == "length") %>% 
    mutate(param = "length:Year1997",
           Parameter = "b_1"),
  mature_results %>% 
    filter(grepl(':+', param)) %>% 
    mutate(est = est + mature_results$est[mature_results$param == "length"],
           Parameter = "b_1")) %>% #View()
  group_by(Parameter) %>% 
  # mutate(scaled_est = scale(est)) %>%
  mutate(scaled_est = (est - mean(est)) / sd(est)) %>% 
  ungroup() %>% 
  mutate(year = rep(1997:max(len_f$year), 2)) -> mature_results; view(mature_results)

# Deviations by year for length-based maturity curve param estimates
ggplot() + 
  geom_segment(data = mature_results %>% 
                 mutate(mycol = ifelse(scaled_est < 0, "blue", "red")),
               aes(x = year, y = 0,
                   xend = year, yend = scaled_est, 
                   color = mycol), size = 2) +
  geom_hline(yintercept = 0, lty = 2) + 
  guides(colour = FALSE) +
  labs(x = "", y = "Scaled parameter estimates") +
  facet_wrap(~ Parameter, ncol = 1)

# Next convert maturity length-based predictions to age using vonB predictions
# from ll survey length-at-age data
age_pred <- seq(0, plus_group, by = 0.01)

vb_pars <- lvb_pars %>% filter(Function == "Length-based LVB" &
                                 Sex == "Female" & 
                                 Source == "ADFG longline survey")
age_pred <- data.frame(age = age_pred,
                       length = round(vb_pars$Estimate[1] * (1 - exp(- vb_pars$Estimate[2] * (age_pred - vb_pars$Estimate[3]))), 1))

# Match lengths back to lengths predicted by vonb for fit_length
pred_simple <- merge(pred_simple, age_pred, by = "length")
which(is.na(pred_simple)) # should be integer(0) head(pred_simple)

# Match lengths back to lengths predicted by vonb for fit_length_year
pred <- merge(pred, age_pred, by = "length") 
which(is.na(pred)) # should be integer(0)
str(pred)
# Get length at 50% maturity (L_50 = -b_0/b_1) and get a_50 from age_pred
# (predicted from vonB)
merge(mature_results %>% 
        select(-scaled_est, -param) %>% 
        pivot_wider(names_from = Parameter, values_from = est) %>%
        mutate(length = - round(b_0 / b_1, 1)), #length at 50% maturity)
      age_pred, by = "length") %>% 
  arrange(year) %>% 
  select(year, l_50 = length, a_50 = age) %>% 
  group_by(year) %>% 
  mutate(a_50 = round(mean(a_50), 1)) %>%
  distinct() %>% ungroup() %>% 
  mutate(mu_a_50 = mean(a_50),
         mu_l_50 = mean(l_50)) -> mat_50_year
data.frame(mat_50_year)
# trends in L50 and a50 by year
ggplot(mat_50_year, aes(x = year, y = l_50)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mu_l_50), lty = 2)

ggplot(mat_50_year, aes(x = year, y = a_50)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = mu_a_50), lty = 2)

merge(pred %>% mutate(year = Year), mat_50_year, by = "year") -> pred

# Age-based maturity curves estimated from length-based maturity and vonB growth
# curve (light blue lines are annual mean predictions, dark blue is the mean)
ggplot() +
  geom_line(data = pred, aes(x = age, y = fitted, group = Year, colour = as.numeric(as.character(Year)))) +  
  geom_line(data = pred_simple, aes(x = age, y = fitted, lty = "All years combined"),
            colour = "black", size = 1) +
  scale_colour_gradientn(colours=rainbow(4)) +
  lims(x = c(0, 20)) +
  scale_linetype_manual(values = 2) +
  labs(x = "\nAge (yr)", y = "Probability\n", colour = "Year", lty = NULL) +
  theme(legend.position = c(.8, .4)) 

ggsave(paste0(YEAR+1,"/figures/maturity_atage_byyear_srvfem_", YEAR, ".png"), 
       dpi=300, height=4, width=6, units="in")

# Comparison with age-based maturity curve
fit_age_year <- glm(Mature ~ age * Year, data = len_f, family = binomial)

# New df for prediction
new_f <- data.frame(age = seq(0, 30, by = 0.01), n_distinct(laa_f$year),
                    Year = factor(sort(rep(unique(laa_f$year), 
                                           length(seq(0, 30, by = 0.01))), 
                                       decreasing = FALSE)))

# Get predicted values by year and take the mean           
broom::augment(x = fit_age_year, 
               newdata = new_f, 
               type.predict = "response") %>% 
  select(Year, age, fitted = .fitted) %>% #, se =.se.fit # prev versions of broom output se
  group_by(age) %>% 
  # Just use mean for illustrative purposes
  mutate(Probability = mean(fitted)) -> pred_age

# Comparison of maturity at age curves. Blue is derived from length-based
# maturity cuve, red is estimated directly from age. Light lines are annual mean
# predictions, dark lines are the means.
ggplot() +
  geom_line(data = pred, aes(x = age, y = fitted, group = Year), col = "lightblue") +  
  geom_line(data = pred_simple, aes(x = age, y = fitted),
            colour = "blue", size = 2, alpha = 0.5) +
  lims(x = c(0, 20)) +
  # scale_linetype_manual(values = 2) +
  labs(x = "\nAge (yr)", y = "Probability\n") +
  theme(legend.position = c(.8, .4)) +
  geom_line(data = pred_age,
            aes(x = age, y = fitted, group = Year), 
            colour = "lightpink", alpha = 0.5) +
  geom_line(data = pred_age,
            aes(x = age, y = Probability), 
            colour = "red", size = 2, alpha = 0.5) +
  labs(x = "Age", y = "Probability")

# Length-based (translated to age; aka the blue one) is more realistic than
# age-based (the red one). Also there is no clear reason to choose the more
# complicated model (fit_length_year) over the simpler model (fit_length)
### pj2023: trend of earlier maturity continues in last two years with the influx
#           of those young year-classes that are dominating the biomass.  
#           Something to try/consider is year specific maturity curves...?

# Maturity at age for YPR and SCAA models
pred_simple %>%  
  filter(age %in% c(rec_age:plus_group)) %>%
  right_join(data.frame(age = rec_age:plus_group)) %>%
  arrange(age) %>% 
  # interpolate fitted probability to fill in any missing values - feel free to
  # revisit rounding if so desired
  mutate(Sex = "Female",
         Source = "LL survey",
         # changed rounding from 2 to 4 in 2020
         probability = round(zoo::na.approx(fitted, maxgap = 2, rule = 2), 4)) %>% 
  select(age, probability) %>% 
  # arrange(age) %>% 
  write_csv(paste0(YEAR+1,"/output/fem_maturityatage_llsrv_plsgrp", plus_group, "_", YEAR, ".csv"))

#Derive age at 50% maturity and kmat (slope of logistic curve)
b0 <- fit_length$coefficients[1]
b1 <- fit_length$coefficients[2]
(L50 <- round(-b0/b1, 1))
(a50 <- age_pred %>% 
    right_join(data.frame(length = L50)) %>% 
    group_by(length) %>% 
    dplyr::summarise(a50 = round(mean(age), 1)) %>% 
    pull(a50))
(kmat <- round(((coef(fit_length)[1] + coef(fit_length)[2]*len) / (len - L50))[1], 2))

data.frame(year_updated = YEAR,
           L50 = L50,
           kmat = kmat,
           a50 = a50) %>% 
  write_csv(paste0(YEAR+1,"/output/maturity_param_", YEAR))  #should this be a csv file? 

#if we wanted to update this to using year specific maturity curves... 
fit_length_year <- glm(Mature ~ length * Year, data = len_f, family = binomial)

b0r <- fit_length_year$coefficients[1]
b1r <- fit_length_year$coefficients[2]
b2r <- c(0,fit_length_year$coefficients[c(3:(2+(YEAR-1997)))])#Year effect
b3r <- c(0,fit_length_year$coefficients[c((3+(YEAR-1997)):(length(fit_length_year$coefficients)))]) #year * length interaction

#need to rerun and round for predicted values... 
new_len_f2 <- data.frame(length = rep(seq(0, 120, 0.0005), n_distinct(len_f$year)),
                         Year = factor(sort(rep(unique(len_f$year), length(seq(0, 120, 0.0005))), decreasing = FALSE)))
broom::augment(x = fit_length_year, 
               newdata = new_len_f2, 
               type.predict = "response") %>% 
  select(Year, length, fitted = .fitted) -> pred2
pred2$length = round(pred2$length,3)
pred2 <- merge(pred2, age_pred, by = "length") 
merge(pred2 %>% mutate(year = Year), mat_50_year, by = "year") -> pred2


round(-(b0r+b2r[1])/(b1r+ b3r[1]), 1) #1997
round(-(b0r+b2r[2])/(b1r+ b3r[2]), 1) 
round(-(b0r+b2r[3])/(b1r+ b3r[3]), 1) 
round(-(b0r+b2r[24])/(b1r+ b3r[24]), 1)
round(-(b0r+b2r[25])/(b1r+ b3r[25]), 1)
round(-(b0r+b2r[26])/(b1r+ b3r[26]), 1) #2022

L50r <- round(-(b0r+b2r)/(b1r+ b3r), 1) #all years

L50df<-data.frame(length = L50r,
                  year = as.factor(replace_na(as.numeric(gsub("[^0-9]","",names(L50r))),1997)))

a50r <- pred2 %>% 
  right_join(data.frame(L50df)) %>% 
  group_by(year,length) %>% 
  dplyr::summarise(a50r = round(mean(age), 1)) %>% 
  pull(a50r)
#******* Need to check if katr1 is correctly calculated? 
kmatr <- round(((b0r + b2r + b3r+ (b1r)*len_r) / (len_r - L50r))[c(1:length(L50r))], 2)
#kmatr2 <- round(((b0r + b2r + b3r* (b1r)*len_r) / (len_r - L50r))[c(1:length(L50r))], 2)
#kmatr3 <- round(((b0r + b2r + (b3r*b1r)*len_r) / (len_r - L50r))[c(1:length(L50r))], 2)

data.frame(year = seq(1997,YEAR,1),
           year_updated = YEAR,
           L50 = L50r,
           kmat = kmatr,
           a50 = a50r) %>% 
  write_csv(paste0(YEAR+1,"/output/maturity_param_byyear_", YEAR))


# # Equation text for plotting values of a_50 and kmat
# a50_txt <- as.character(
#   as.expression(substitute(
#     paste(italic(a[50]), " = ", xx),
#     list(xx = formatC(a50, format = "f", digits = 1)))))
# 
# kmat_txt <- as.character(
#   as.expression(substitute(
#     paste(italic(k[mat]), " = ", xx),
#     list(xx = formatC(kmat, format = "f", digits = 1)))))


}
#----------------------------------------------------------------------------------
# Sex ration
# restrict age range
aa <- rec_age:plus_group

Bio_lw_noOL %>% 
  filter(age %in% aa) %>% 
  ungroup() %>% 
  dplyr::select(Sex, age) %>% 
  dplyr::filter(Sex %in% c("Female", "Male")) %>% 
  na.omit() %>% 
  droplevels() %>% 
  dplyr::count(Sex, age) %>% 
  dplyr::group_by(age) %>% 
  dplyr::mutate(proportion = round(n / sum(n), 2),
         Source = "fishery") %>% 
  filter(Sex == "Female") -> byage
str(byage); view(byage)

min(byage$age)
aa <- min(byage$age):plus_group

byage<-rbind(byage,pot_missing_ages) %>% 
  arrange(factor(Source, levels = c("LL survey","LL fishery","Pot fishery")),
          age)
# get generalized additive model fits and predictions
# survey
library(gamm4)
fitage <- gam(I(Sex == "Female") ~ s(age, k=4), 
                  data = filter(Bio_lw_noOL, age %in% aa, 
                                Sex %in% c("Female", "Male")),
                  family = "quasibinomial")
summary(fitage)
plot(fitage)

predage <- predict(fitage, newdata = data.frame(age = aa),
                       type = "response", se = TRUE)

# combine with the sex_ratio df
# *FLAG* bind_cols goes by col position, make sure survey is first

length(unique(byage$age))

bind_cols(  #1st year pot data does not have all ages... 
  byage,bind_rows(as_tibble(do.call(cbind, predage)) %>% 
                    dplyr::mutate(source_check = "fishery")))-> byage
  #do.call cbinds each vector in the predict() output list 
  #bind_rows(as_tibble(do.call(cbind, predage)) %>% 
  #            mutate(source_check = "fishery"),
  #          as_tibble(do.call(cbind, fsh_ll_predage)) %>% 
  #            mutate(source_check = "LL fishery"),
  #          as_tibble(do.call(cbind, fsh_pot_predage)) %>% 
  #            mutate(source_check = "Pot fishery")) 

view(byage)
view(byage_simp)#not sure what this is

# plot
ggplot(byage, aes(x = age)) +
  geom_line(aes(y = fit, col = Source)) +
  geom_ribbon(aes(ymin = fit - se.fit*2, ymax = fit + se.fit*2, fill = Source),  alpha = 0.2) +
  geom_point(aes(y = proportion, col = Source)) +  
  expand_limits(y = c(0.0, 1)) +
  xlab("\nAge") +
  ylab("Proportion of females\n") +
  geom_hline(yintercept = 0.5, lty = 2, col = "grey") +
  scale_colour_viridis_d(option="D",begin=0,end=0.75) +
  scale_fill_viridis_d(option="D",begin=0,end=0.75) +
  #scale_fill_manual(values = c("black", "grey")) +
  theme(legend.position = c(0.85, 0.3),
        legend.box.background = element_rect(fill = "transparent",
                                             colour = "transparent")) -> byage_plot

# proportion of females by year in the fishery and survey 

Bio_lw_noOL %>% 
  filter(age %in% aa) %>% 
  ungroup() %>% 
  select(Sex, year) %>% 
  filter(Sex %in% c("Female", "Male")) %>% 
  na.omit() %>% 
  droplevels() %>% 
  dplyr::count(Sex, year) %>% 
  group_by(year) %>% 
  mutate(proportion = round(n / sum(n), 2),
         Source = "fishery") %>% 
  filter(Sex == "Female") -> byyear_base
view(byyear_base)

Bio_lw_noOL %>% 
  filter(age %in% aa) %>% 
  ungroup() %>% 
  select(Sex, year) %>% 
  filter(Sex %in% c("Female", "Male")) %>% 
  na.omit() %>% 
  droplevels() %>% 
  dplyr::count(Sex, year) %>% 
  group_by(year) %>% 
  mutate(proportion = round(n / sum(n), 2),
         Source = "fishery") %>% 
  filter(Sex == "Female") -> byyear
view(byyear)

# Save output for YPR analysis
write_csv(byyear, paste0("Output/sexratio_byyear_plsgrp", plus_group, "_", YEAR, ".csv"))

# get generalized additive model fits and predictions
# survey
fityear <- gam(I(Sex == "Female") ~ s(year, k = 4), 
                   data = filter(Bio_lw_noOL, Sex %in% c("Female", "Male")),
                   family = "quasibinomial") #quasi- allows non-integers and avoids warnings...
summary(fityear); plot(fityear)

fsh_yrs <- byyear %>% select(year) %>% range()

fsh_predyear <- predict(fityear, 
                        newdata = data.frame(year = min(fsh_yrs):max(fsh_yrs) ),
                        type = "response", se = TRUE)

# combine with the sex_ratio df
# *FLAG* bind_cols goes by col position, make sure survey is first
length(unique(byyear$year)); unique(byyear$Sex)

#some missing years we need to fill in... 
yearlist<- seq(min(byyear$year),YEAR,1)
missingyears<-yearlist %>% data.frame %>% filter(!(yearlist %in% byyear$year))
missing_years<-cbind(Sex = rep("Female",nrow(missingyears)),
                        year = missingyears,
                        n = rep(0,nrow(missingyears)),
                        proportion = rep(NA,nrow(missingyears)),
                        Source = rep("fishery",nrow(missingyears))) %>%
dplyr::rename(year = ".")
byyear<-rbind(byyear,missing_years) %>% 
  arrange(year)
#  arrange(factor(Source, levels = c("LL survey","LL fishery","Pot fishery")),
#          age)


bind_cols(
  byyear,# %>% filter(!(Source == "Pot fishery")),
  #do.call cbinds each vector in the predict() output list 
  bind_rows(tibble::as_tibble(do.call(cbind, fsh_predyear))%>% mutate(source_check = "fishery"))) %>% 
  mutate(proportion = as.numeric(proportion))-> byyear
str(byyear)


# plots

ggplot(data = byyear, aes(x = year)) +
  geom_line(aes(y = fit, col = Source), size = 1) +
  geom_ribbon(aes(ymin = fit - se.fit*2, ymax = fit + se.fit*2, 
                  fill = Source),  alpha = 0.2) +
  geom_point(data=byyear, aes(y = proportion, col = Source)) +  
  #geom_point(data=byyear %>% filter(Source == "Pot fishery"),
  #           aes(x=year,y=proportion, col = "green")) +
  expand_limits(y = c(0.0, 1)) +
  geom_hline(yintercept = 0.5, lty = 2, col = "grey") +
  xlab("\nYear") +
  ylab("Proportion of females\n") +
  scale_colour_viridis_d(option="D",begin=0,end=0.75) +
  scale_fill_viridis_d(option="D",begin=0,end=0.5) +
  #scale_colour_manual(values = c("black", "grey")) +
  #scale_fill_manual(values = c("black", "grey")) +
  theme(legend.position = "none") -> byyear_plot
byyear_plot

library(cowplot)
plot_grid(byage_plot, byyear_plot, align = c("h"), ncol = 1)
ggsave(paste0("Figures/sex_ratios_", YEAR, ".png"), dpi=300,  height=6, width=7, units="in")

#---------------------------------------------------------------------------
# Age compositions ----
colnames(Bio_lw_noOL)
#quos() uses stand eval in dplyr, eval cols with nonstand eval using !!!
unique(Bio_lw_noOL$Project)
cols <- quos(Source, year, Sex, age) 

bind_rows(Bio_lw_noOL %>% filter(Project == "Commercial Longline Trip") %>%
            mutate(Source = "LL fishery") %>% select(!!!cols), 
          Bio_lw_noOL %>% filter(Project == "Atypical Longline Sample") %>%
            mutate(Source = "LL fishery atyp") %>% select(!!!cols),
          Bio_lw_noOL %>% filter(Project == "Commercial Halibut Longline") %>%
            mutate(Source = "Hal LL fishery") %>% select(!!!cols),
          Bio_lw_noOL %>% filter(Project == "Commercial Jig Trip") %>%
            mutate(Source = "Jig") %>% select(!!!cols)) %>%
  filter(Sex %in% c('Female', 'Male') & !is.na(age)) %>% 
  droplevels() %>% 
  mutate(age = ifelse(age >= plus_group, plus_group, age)) %>% 
  filter(age >= 10) -> all_bio  # Plus group
view(all_bio)

# Age comps (sex-specific)
all_bio %>% 
  dplyr::count(Source, Sex, year, age) %>%
  group_by(Source, Sex, year) %>% 
  mutate(proportion = round( n / sum(n), 5),
         tot_samps = sum(n)) %>% 
  bind_rows(all_bio %>% # Age comps (sexes combined)
              dplyr::count(Source, year, age) %>%
              group_by(Source, year) %>% 
              mutate(proportion = round( n / sum(n), 5),
                     Sex = "Sex combined",
                     tot_samps = sum(n))) -> agecomps  

agecomps %>% group_by(year, Source, Sex) %>%
  dplyr::mutate(tot_samps2 = n()) %>% ungroup() -> agecomps2

range(agecomps$proportion)
range(agecomps$age)
range(agecomps$tot_samps)
agecomps %>% filter(Sex == "Sex combined")

eg<-agecomps2 %>% filter(Source == "Jig" & Sex == "Female" & 
                          year == 2002)
sum(eg$proportion)

# Years with pot bio data
# potsrv_bio %>% 
#   filter(!is.na(age) & !is.na(Sex)) %>% 
#   distinct(year) -> pot_yrs

# complete() was behaving weirdly. Expand to grid to include all age combos
library(padr)
expand.grid(year = unique(agecomps$year), 
            Source = unique(agecomps$Source),
            Sex = unique(agecomps$Sex),
            age = seq(10, plus_group, 1))  %>% 
  data.frame()  %>% 
  full_join(agecomps) %>%
  fill_by_value(n, proportion, value = 0) %>% 
  #fill_by_value(tot_samps) %>%
  group_by(year,Source,Sex) %>%
  mutate(tot_samps = ifelse(!is.na(tot_samps),tot_samps,
                            ifelse(length(unique(tot_samps)) > 1,
                                   tot_samps[!is.na(tot_samps)],0))) %>% ungroup() %>%
  mutate(Age = factor(age),
         proportion = round(proportion, 5)) -> agecomps #%>%
  # Keep only relevant years for each Source
  #filter(c(Source == "LL fishery" & year >= 2002) |
  #         c(Source == "LL survey" & year >= 1997) |
  #         c(Source == "Pot fishery" & year >= 2022)) -> agecomps

# Check that they sum to 1
agecomps %>% 
  group_by(Source, Sex, year) %>% 
  summarise(sum(proportion)) %>% 
  print(n = Inf)

# Sample sizes by source/year/sex
agecomps %>% 
  group_by(Source, year, Sex) %>% 
  dplyr::summarize(n = sum(n)) %>% 
  arrange(year) %>% 
  pivot_wider(names_from = year, values_from = n, values_fill = list(n = 0)) %>% 
  write_csv(paste0("Output/n_agecomps_plsgrp", plus_group, "_", YEAR, ".csv"))

# Age comp matrix
agecomps %>% write_csv(paste0("Output/agecomps_plsgrp", plus_group, "_", YEAR, ".csv"))

# Bargraph for presentation
agecomps %>% 
  filter(#year == YEAR & 
           tot_samps > 50,
           Source == "LL fishery" &
           Sex %in% c("Male", "Female")) %>% 
  ggplot(aes(age, proportion, fill = Sex, col=Sex), alpha = 0.5) +
  geom_bar(stat = "identity",
           position = "dodge") +
  facet_wrap(~year) +
  # position = position_dodge(preserve = "single")) +
  # scale_fill_grey(start = 0.3, end = 0.8) +
  scale_x_continuous(breaks = seq(min(agecomps$age), max(agecomps$age), 10), 
                     labels =  seq(min(agecomps$age), max(agecomps$age), 10)) +
  labs(x = "\nAge", y = "Proportion\n") +
  theme(legend.position = c(0.9, 0.1)) +
  ggtitle(paste0("Longline fishery"))

ggsave(paste0("Figures/agecomp_bargraph_ll_fsh_", YEAR, ".png"), 
       dpi=300, height=3, width=9, units="in")

agecomps %>% 
  filter(#year == YEAR & 
    tot_samps > 50,
    Source == "Hal LL fishery" &
      Sex %in% c("Male", "Female")) %>% 
  ggplot(aes(age, proportion, fill = Sex, col=Sex), alpha = 0.5) +
  geom_bar(stat = "identity",
           position = "dodge") +
  facet_wrap(~year) +
  # position = position_dodge(preserve = "single")) +
  # scale_fill_grey(start = 0.3, end = 0.8) +
  scale_x_continuous(breaks = seq(min(agecomps$age), max(agecomps$age), 10), 
                     labels =  seq(min(agecomps$age), max(agecomps$age), 10)) +
  labs(x = "\nAge", y = "Proportion\n") +
  theme(legend.position = c(0.9, 0.1)) +
  ggtitle(paste0("Halibut longline fishery"))

ggsave(paste0("Figures/agecomp_bargraph_hal_fsh_", YEAR, ".png"), 
       dpi=300, height=3, width=9, units="in")

# All years smoothed by source
agecomps %>% 
  filter(age < plus_group & Sex == "Sex combined" & tot_samps > 50) %>% 
  ggplot(aes(x = age, y = proportion, colour = Source, linetype = Source)) +
  geom_point(size = 1, alpha = 0.1) +
  stat_smooth(size = 1, se = FALSE) +
  scale_colour_grey() +
  # scale_colour_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) +
  scale_y_continuous(limits = c(0, 0.1),
                     breaks = round(seq(min(agecomps$proportion), 0.1, 0.02), 2), 
                     labels =  round(seq(min(agecomps$proportion), 0.1, 0.02), 2)) +
  xlab('\nAge') +
  ylab('Proportion\n') +
  theme(legend.position = c(0.8, 0.8))

ggsave(paste0("Figures/agecomp_bydatasource.png"), 
       dpi=300, height=5, width=5, units="in")

## LL fishery and Halibut longline fishery are very similar.  Jig fishery data is VERY limited
## Not sure what the deal is with the atypical longline fishery.. but limited data and very different
## length distributions seen
## If I was running an ASA now I would combine the long line data and ditch the
## Jig and atypical fishery... 

bind_rows(Bio_lw_noOL %>% 
            filter(Project == "Commercial Longline Trip") %>%
            mutate(Source = "LL fishery") %>% 
            select(!!!cols), 
          Bio_lw_noOL %>% 
            filter(Project == "Commercial Halibut Longline") %>%
            mutate(Source = "LL fishery") %>% 
            select(!!!cols)) %>%
  filter(Sex %in% c('Female', 'Male') & !is.na(age)) %>% 
  droplevels() %>% 
  mutate(age = ifelse(age >= plus_group, plus_group, age)) %>% 
  filter(age >= 10) -> ll_bio  # Plus group
#view(all_bio)
unique(ll_bio$Source)

# Age comps (sex-specific)
ll_bio %>% 
  dplyr::count(Source, Sex, year, age) %>%
  group_by(Source, Sex, year) %>% 
  mutate(proportion = round( n / sum(n), 5),
         tot_samps = sum(n)) %>% 
  bind_rows(ll_bio %>% # Age comps (sexes combined)
              dplyr::count(Source, year, age) %>%
              group_by(Source, year) %>% 
              mutate(proportion = round( n / sum(n), 5),
                     Sex = "Sex combined",
                     tot_samps = sum(n))) -> agecomps_ll  

expand.grid(year = unique(agecomps_ll$year), 
            Source = unique(agecomps_ll$Source),
            Sex = unique(agecomps_ll$Sex),
            age = seq(10, plus_group, 1))  %>% 
  data.frame()  %>% 
  full_join(agecomps_ll) %>%
  fill_by_value(n, proportion, value = 0) %>% 
  #fill_by_value(tot_samps) %>%
  group_by(year,Source,Sex) %>%
  mutate(tot_samps = ifelse(!is.na(tot_samps),tot_samps,
                            ifelse(length(unique(tot_samps)) > 1,
                                   tot_samps[!is.na(tot_samps)],0))) %>% ungroup() %>%
  mutate(Age = factor(age),
         proportion = round(proportion, 5)) -> agecomps_ll

agecomps_ll %>% 
  group_by(Source, year, Sex) %>% 
  dplyr::summarize(n = sum(n)) %>% 
  arrange(year) %>% 
  pivot_wider(names_from = year, values_from = n, values_fill = list(n = 0)) %>% 
  write_csv(paste0("Output/combll_n_agecomps_plsgrp", plus_group, "_", YEAR, ".csv"))

# Age comp matrix
agecomps_ll %>% write_csv(paste0("Output/combll_agecomps_plsgrp", plus_group, "_", YEAR, ".csv"))

agecomps_ll %>% 
  filter(#year == YEAR & 
    tot_samps > 50,
    #Source == "LL fishery" &
      Sex %in% c("Male", "Female")) %>% 
  ggplot(aes(age, proportion, fill = Sex, col=Sex), alpha = 0.5) +
  geom_bar(stat = "identity",
           position = "dodge") +
  facet_wrap(~year) +
  # position = position_dodge(preserve = "single")) +
  # scale_fill_grey(start = 0.3, end = 0.8) +
  scale_x_continuous(breaks = seq(min(agecomps$age), max(agecomps$age), 10), 
                     labels =  seq(min(agecomps$age), max(agecomps$age), 10)) +
  labs(x = "\nAge", y = "Proportion\n") +
  theme(legend.position = c(0.9, 0.1)) +
  ggtitle(paste0("Halibut and general longline fishery"))

ggsave(paste0("Figures/agecomp_bargraph_combll_fsh_", YEAR, ".png"), 
       dpi=300, height=3, width=9, units="in")

# bubble plots filled circles

ggplot(data = agecomps_ll %>% filter(tot_samps > 50),
       aes(x = year, y = age, size = proportion)) + #*FLAG* could swap size with proportion_scaled
  geom_point(shape = 21, fill = "black", colour = "black") +
  scale_size(range = c(0, 4)) +
  facet_wrap(~ Sex) +
  labs(x = "\nYear", y = "Observed Age\n") +
  guides(size = FALSE)+
  theme_classic()

ggsave(paste0("Figures/bubble_combll_agecomp_byyear.png"), dpi=300, height=5, width=7.5, units="in")

#---------------------------------------------------------------------------
# Length compositions
Bio_lw_noOL

histogram(Bio_lw_noOL$length)
range(Bio_lw_noOL$length)

bind_rows(Bio_lw_noOL %>% filter(Project == "Commercial Longline Trip") %>%
            mutate(Source = "LL fishery"), 
          Bio_lw_noOL %>% filter(Project == "Atypical Longline Sample") %>%
            mutate(Source = "LL fishery atyp"),
          Bio_lw_noOL %>% filter(Project == "Commercial Halibut Longline") %>%
            mutate(Source = "Hal LL fishery"),
          Bio_lw_noOL %>% filter(Project == "Commercial Jig Trip") %>%
            mutate(Source = "Jig")) %>% 
  filter(#year >= 1997 &
    Sex %in% c("Female", "Male") &
      !is.na(length),
  ) %>% 
  select(year, Sex, length, Source) %>%
#  filter(!c(length < 40)) %>% 
  mutate(length2 = ifelse(length < 300, 300,
                          ifelse(length > 860, 860, length)),
         length_bin = cut(length2, breaks = seq(299.9, 859.9, 20),
                          labels = paste(seq(310, 860, 20)))) %>% 
  select(-length2) -> lendat2

lendat2 %>% 
  # Length comps by Source, year, and Sex 
  count(Source, Sex, year, length_bin) %>%
  group_by(Source, Sex, year) %>% 
  #  mutate(proportion = round( n / sum(n), 4)) %>% 
  mutate(proportion = n / sum(n),
         tot_samps = sum(n)) %>%
  bind_rows(lendat2 %>% # Sexes combined
              count(Source, year, length_bin) %>%
              group_by(Source, year) %>% 
              mutate(proportion = round( n / sum(n), 4),
                     Sex = "Sex combined",
                     tot_samps = sum(n))) -> lencomps

# complete() was behaving weirdly. Expand to grid to include all length combos
expand.grid(year = unique(lencomps$year), 
            Source = unique(lencomps$Source),
            Sex = unique(lencomps$Sex),
            length_bin = sort(unique(lendat$length_bin)))  %>% 
  data.frame()  %>% 
  full_join(lencomps) %>%
  fill_by_value(n, proportion, value = 0) %>% 
  group_by(year,Source,Sex) %>%
  mutate(tot_samps = ifelse(!is.na(tot_samps),tot_samps,
                            ifelse(length(unique(tot_samps)) > 1,
                                   tot_samps[!is.na(tot_samps)],0))) %>% ungroup() %>%
  mutate(#length_bin = factor(length_bin),
    proportion = round(proportion, 6)) -> lencomps # %>%
  # Keep only relevant years for each Source
  #filter(c(Source == "LL fishery" & year >= 2002) |
  #         c(Source == "LL survey" & year >= 1997) |
  #         c(Source == "Pot fishery" & year >= 2022) #|
         # c(Source == "Pot survey" & year %in% pot_yrs$year)
  #) -> lencomps

# Check that they sum to 1
lencomps %>% 
  group_by(Source, Sex, year) %>% 
  summarise(sum(proportion)) %>% View()

write_csv(lencomps, paste0("Output/lengthcomps_", YEAR, ".csv"))

str(lendat)
lendat2 %>% 
  # Mean length comp for comparison
  count(Source, Sex, length_bin) %>%
  group_by(Source, Sex) %>% 
  mutate(proportion = round( n / sum(n), 4)) %>% 
  arrange(Source, Sex, length_bin)  %>% 
  data.frame() %>%  
  complete(Source, length_bin,
           fill = list(n = 0, proportion = 0)) %>% 
  bind_rows(lendat2 %>% # Sexes combined
              count(Source, length_bin) %>%
              group_by(Source) %>% 
              mutate(proportion = round( n / sum(n), 4),
                     Sex = "Sex combined")) -> mu_lencomps

mu_lencomps %>% 
  group_by(Source, Sex) %>% 
  summarise(sum(proportion))

ll_lencomps <- lencomps %>% 
  filter(Source == "LL fishery")

hal_lencomps <- lencomps %>% 
  filter(Source == "Hal_LL fishery")

# ggridge plots

lendat2 %>% 
  #filter(Source != "Pot survey") %>% 
  #mutate(Source = derivedFactor("Survey" = Source == "LL survey",
  #                              "LL Fishery" = Source == "LL fishery",
  #                              "Pot Fishery" = Source == "Pot fishery",
  #                              .ordered = TRUE)) %>% 
  ggplot(aes(length, year, group = year, fill = year)) + 
  geom_density_ridges(aes(point_fill = year, point_color = year),
                      alpha = 0.3) +
  geom_vline(xintercept = 63, linetype = 4) +
  xlim(300, 900) + 
  xlab("\nLength (cm)") + 
  ylab(NULL) +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  facet_wrap(~ Source)

ggsave(paste0("Figures/lengthcomp_ggridges_", YEAR, ".png"), 
       dpi=300, height=8, width=10, units="in")

# ggride plot for len dat by sex (for TMB inputs)

lendat2 %>% 
  #filter(! c(Source %in% c("Pot survey", "LL fishery", "Pot fishery"))) %>% 
  #mutate(Source = derivedFactor("Survey" = Source == "LL survey")) %>% 
  ggplot(aes(length, year, group = year, fill = year)) + 
  geom_density_ridges(aes(point_fill = year, point_color = year), alpha = 0.3) +
  #geom_vline(xintercept = 61, linetype = 4) + # L50
  xlim(300, 900) + 
  labs(x = "\nLength (cm)", y = "Year\n") +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  facet_wrap(~ Sex) +
  ggtitle("Fishery")

ggsave(paste0(YEAR+1,"/figures/tmb/lencomp_fsh_",YEAR,".png"), 
       dpi=300, height=8, width=10, units="in")

# All years smoothed by source
ggplot() +
  geom_point(data = lencomps %>% 
               filter(Sex == "Sex combined"),
             aes(x = length_bin, y = proportion, 
                 colour = Source),
             size = 1, alpha = 0.2) +
  # stat_smooth(size = 1.5, se = FALSE) +
  geom_line(data = mu_lencomps %>% 
              filter(Sex == "Sex combined"),
            aes(x = length_bin, y = proportion, colour = Source, 
                group = Source, linetype = Source), size = 1) +
  scale_x_discrete(breaks = seq(300, 900, 20),
                   labels = seq(300, 900, 20)) +
  #scale_colour_grey() +
  scale_colour_viridis_d(option = "A", begin=0, end=0.75) +
  xlab('\nFork length (mm)') +
  ylab('Proportion\n') +
  theme(legend.position = c(0.8, 0.8))

ggsave(YEAR+1,"/figures/lengthcomp_bydatasource.png", 
       dpi=300, height=4.5, width=5, units="in")

# length comp figs, requested by AJ Lindley 2018-09-07
lencomps %>% 
  group_by(year, Source, Sex) %>% 
  dplyr::summarize(N = sum(n),
                   label = paste0("n = ", prettyNum(N, big.mark = ","))) %>% 
  ungroup() %>% 
  mutate(length_bin = "91", proportion = 0.18) -> labels 

# For fishery

ggplot(data = lencomps %>% 
         # Last 10 years of data
         filter(year >= YEAR - 10 & 
                  Sex != "Sex combined" &
                  Source == "LL fishery"), 
       aes(x = length_bin, y = proportion)) + 
  geom_bar(stat = "identity", colour = "lightgrey", fill = "lightgrey", width = 0.8) +
  geom_line(data = lencomps %>% 
              # Compare all past years to this year
              filter(year == YEAR & 
                       Sex != "Sex combined" &
                      Source == "LL fishery") %>% 
              select(-year),
            aes(x = length_bin, y = proportion, group = 1),
            colour = "black") +
  geom_text(data = labels %>% 
              filter(year >= YEAR - 10 & 
                       Sex != "Sex combined" &
                       Source == "LL fishery"),
            aes(x = length_bin, y = proportion, label = label),
            size = 3, family = "Times") +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = round(seq(0, 0.2, 0.1), 2),
                     labels =  round(seq(0, 0.2, 0.1), 2)) +
  scale_x_discrete(breaks = seq(41, 99, 6),
                   labels = seq(41, 99, 6)) +
  facet_grid(year ~ Sex) +
  labs(x = "\nFork length (cm)", y = "Proportion-at-length (longline fishery)\n") +
  theme(strip.placement = "outside") 

ggsave(paste0(YEAR+1,"/figures/llfsh_lencomps_", YEAR-10, "_", YEAR, ".png"), 
       dpi=300, height=8, width=6.5, units="in")

# Summary stats output for length comps (requested by AJ Linsley 20180907)
lendat2 %>% 
  #filter(Source %in% c("LL survey", "LL fishery","Pot fishery")) %>% 
  group_by(Source, Sex, year) %>% 
  dplyr::summarize(mean = mean(length),
                   min = min(length),
                   max = max(length)) %>% 
  mutate(variable = "Fork length") -> lensum

# axis <- tickr(lensum, year, 5)

lensum %>% 
  ggplot(aes(x = year, y = mean, colour = Source)) +
  geom_point() +
  geom_line() +
  scale_colour_grey(guide = FALSE) +
  facet_wrap(~ Sex) +
  # scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  labs(x = NULL, y = "Mean\nfork\nlength\n(cm)") +
  theme_classic(axis.title.y = element_text(angle=0)) -> l

#quos() uses stand eval in dplyr, eval cols with nonstand eval using !!!
cols <- quos(Source, year, Sex, age) 

# leaving off here on 7-14 
bind_rows(Bio_lw_noOL %>% filter(Project == "Commercial Longline Trip") %>%
            mutate(Source = "LL fishery")%>% select(!!!cols), 
          Bio_lw_noOL %>% filter(Project == "Atypical Longline Sample") %>%
            mutate(Source = "LL fishery atyp")%>% select(!!!cols),
          Bio_lw_noOL %>% filter(Project == "Commercial Halibut Longline") %>%
            mutate(Source = "Hal LL fishery")%>% select(!!!cols),
          Bio_lw_noOL %>% filter(Project == "Commercial Jig Trip") %>%
            mutate(Source = "Jig") %>% select(!!!cols)) %>%
#bind_rows(
#  fsh_bio %>% filter(Gear == "Longline") %>% 
#    mutate(Source = "LL fishery") %>% select(!!!cols), 
#  fsh_bio %>% filter(Gear == "Pot") %>% 
#    mutate(Source = "Pot fishery") %>% select(!!!cols)) %>% 
#  bind_rows(srv_bio %>% mutate(Source = "LL survey") %>% select(!!!cols)) %>% 
  filter(Sex %in% c('Female', 'Male') & !is.na(age)) %>% 
  group_by(Source, Sex, year) %>% 
  dplyr::summarize(mean = mean(age),
                   min = min(age),
                   max = max(age)) %>% 
  mutate(variable = "Age") -> agesum

agesum %>% 
  mutate(age = round(mean, 0)) -> agesum1

# axisy <- tickr(agesum1, age, 3)

agesum %>% 
  ggplot(aes(x = year, y = mean, colour = Source)) +
  geom_point() +
  geom_line() +
  scale_colour_grey() +
  facet_wrap(~ Sex) +
  # scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  # scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) +
  labs(x = NULL, y = "Mean\nage\n(yrs)") +
  theme(legend.position = "bottom",
        axis.title.y = element_text(angle=0)) -> a

cowplot::plot_grid(l, a, axis = "lrtb", align = "hv", ncol = 1) -> compare_comp_sums
compare_comp_sums
ggsave(paste0("Figures/compare_comp_summaries.png"),
       plot = compare_comp_sums,
       # dpi=300, height=5.5, width=6.5, units="in")
       dpi=300, height=7, width=7, units="in")

bind_rows(agesum, lensum) %>% 
  write_csv(paste0("Output/comps_summary.csv"))

# survey length comps by stat area, requested by A Olson 2021-02-18
unique(Bio_lw_noOL$Stat)

Bio_lw_noOL %>% 
  filter(#year >= 1997 &
           Sex %in% c("Female", "Male") &
           !is.na(length)) %>% 
  select(year, Stat, length) %>% 
 # filter(!c(length < 40)) %>% 
  mutate(length2 = ifelse(length < 300, 300,
                          ifelse(length > 860, 860, length)),
         length_bin = cut(length2, breaks = seq(299.9, 859.9, 20),
                          labels = paste(seq(310, 860, 20)))) %>% 
  select(-length2) -> statlen #%>% 
#  mutate(Stat = derivedFactor("345731" = Stat == "345731",
#                              "345701" = Stat == "345701",
#                              "345631" = Stat == "345631",
#                              "345603" = Stat == "345603",
#                              .ordered = TRUE)) -> statlen

statlen %>% 
  # Length comps by Source, year, and Sex 
  count(Stat, year, length_bin) %>%
  group_by(Stat, year) %>% 
  mutate(proportion = round( n / sum(n), 4)) -> statcomps

# complete() was behaving weirdly. Expand to grid to include all length combos
expand.grid(year = unique(statcomps$year), 
            Stat = unique(statcomps$Stat),
            length_bin = sort(unique(statlen$length_bin)))  %>% 
  data.frame()  %>% 
  full_join(statcomps) %>%
  fill_by_value(n, proportion, value = 0) %>% 
  mutate(proportion = round(proportion, 4)) -> statcomps

# Check that they sum to 1
statcomps %>% 
  group_by(Stat, year) %>% 
  summarise(sum(proportion)) %>% view

statcomps %>% 
  group_by(year, Stat) %>% 
  dplyr::summarize(N = sum(n),
                   label = paste0("n = ", prettyNum(N, big.mark = ","))) %>% 
  ungroup() %>% 
  mutate(length_bin = "91", proportion = 0.18) -> labels 

# For survey
ggplot(data = statcomps %>% 
         # Last 10 years of data
         filter(year >= YEAR - 10), 
       aes(x = length_bin, y = proportion)) + 
  geom_bar(stat = "identity", colour = "lightgrey", fill = "lightgrey", width = 0.8) +
  geom_line(data = statcomps %>% 
              # Compare all past years to this year
              filter(year == YEAR) %>% 
              select(-year),
            aes(x = length_bin, y = proportion, group = 1),
            colour = "black") +
  geom_text(data = labels %>% 
              filter(year >= YEAR - 10),
            aes(x = length_bin, y = proportion, label = label),
            size = 3, family = "Times") +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = round(seq(0, 0.2, 0.1), 2),
                     labels =  round(seq(0, 0.2, 0.1), 2)) +
  scale_x_discrete(breaks = seq(310, 860, 20),
                   labels = seq(310, 860, 20)) +
  facet_grid(year ~ Stat) +
  labs(x = "\nFork length (mm)", y = "Proportion-at-length (commercial fisheries)\n") +
  theme(strip.placement = "outside") 

# ggridge plots

statlen %>% 
  filter(year >= YEAR - 5) %>% 
  ggplot(aes(length, year, group = year, fill = year)) + 
  geom_density_ridges(aes(point_fill = year, point_color = year),
                      alpha = 0.3) +
  xlim(300, 900) + 
  xlab("\nLength (mm)") + 
  ylab(NULL) +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  facet_wrap(~ Stat, ncol = 1)

statlen %>%
  filter(year >= YEAR - 10 & between(length, 40, 99)) %>%
  ggplot(aes(x = length, y = factor(year), group = interaction(year, Stat),
             fill = Stat)) + #, col = Stat)) +
  geom_density_ridges(alpha = 0.5, size = 0.5, col = "white") +
  # scale_fill_manual(values = c("grey30", "#00BFC4", "#da2055", "#daa520" )) +
  # scale_colour_manual(values = c("grey30", "#00BFC4", "#daa520", "#da2055")) +
  labs(x = "Length (cm)", y = NULL, fill = NULL, col = NULL) + #, title = "Dusky rockfish") +
  theme_light() +
  scale_x_continuous(limits = c(300, 900)) +
  theme(legend.position = "top")

ggsave(paste0("Figures/lengthcomp_statarea_", YEAR, ".png"), 
       dpi=300, height=8, width=10, units="in")


























#--------------------------------------------------------------------------
## END END END END

############################################################################

############################################################################
############################################################################
## Plotting Scrap
############################################################################
############################################################################
## INDIVIDUAL PLOTS
## Age Bubble plots

# bubble plots filled circles

# survey
agecompdat <- agecomps %>% 
  filter(Sex %in% c("Female", "Male") &
           #          Source %in% c("LL survey") &
           GFMU == "SSEO" &
           #Groundfish.Management.Area.Code == "EYKT" &
           Year >= 1990 &
           Age <= plus_group) %>% 
  ungroup()

# axisx <- tickr(agecompdat, year, 3)
# axisy <- tickr(agecompdat, age, 5)

ggplot(data = agecompdat,
       aes(x = Year, y = Age, size = count)) + #*FLAG* could swap size with proportion_scaled
  geom_point(shape = 21, fill = "black", colour = "black") +
  scale_size(range = c(0, 2)) +
  facet_wrap(~ Sex ) +
  labs(x = "Year", y = "Observed Age") +	#"\nYear", y = "Observed Age\n") +
  guides(size = "none")	+	#FALSE)# +
  #  theme_classic() +
  theme(text=element_text(size=12),
        axis.text.x = element_text(size=10,angle=45, hjust=1),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=15)) +	#, labels = axisx$labels) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=15)) 


#   scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
#   scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels)

ggsave("Figures/SSEO_bubble_port_agecomp_byyear.png", dpi=900, height=5, width=7, units="in")

# fishery
agecompdat <- agecomps %>% 
  filter(Sex %in% c("Female", "Male") &
           Source %in% c("LL fishery") &
           year >= 2002 &
           age <= plus_group) %>% 
  ungroup()

# axisx <- tickr(agecompdat, year, 3)
# axisy <- tickr(agecompdat, age, 5)

ggplot(data = agecompdat,
       aes(x = year, y = age, size = proportion)) + #*FLAG* could swap size with proportion_scaled
  geom_point(shape = 21, colour = "black", fill = "black") +
  scale_size(range = c(0, 4)) +
  facet_wrap(~ Sex) +
  labs(x = "\nYear", y = "Observed age\n") +
  guides(size = FALSE) #+
# scale_x_continuous(breaks = axisx$breaks, labels = axisx$labels) +
# scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels)

ggsave("figures/bubble_fishery_agecomp_byyear.png", 
       dpi=300, height=5, width=7.5, units="in")

##############################################################################################
# Length compositions ----

# Pers. comm. K. Fenske 2018-01-05: NMFS uses length bins 41, 43, 45 ... 99.
# These bins represent the center of the bin, so a 43 bin represents fish
# 42-43.9 cm. They omit fish smaller than 40 and fish larger than 100 cm are
# lumped into the 100 bin. I've maintained these conventions for easy
# comparison:
##Phil's try!

##JITTER OPTION		
lendat %>% 
  #filter(Groundfish.Management.Area.Code == "EYKT") %>% 
  filter(GFMU == "SSEO") %>% 
  #  filter(n > 150) %>%
  #  mutate(Source = derivedFactor("Survey" = Source == "LL survey",
  #                                "Fishery" = Source == "LL fishery",
  #                                .ordered = TRUE)) %>% 
  ggplot(aes(length, Year, group = Year, fill = Year)) + 
  geom_density_ridges(aes(point_fill = Year, point_color = Year),
                      alpha = 0.3, scale=1.5, jittered_points = TRUE,
                      rel_min_height = 0.01) +
  #  geom_vline(xintercept = 400, linetype = 4) +
  xlim(300, 900) + 
  xlab("Length (mm)") + 
  ylab(NULL) +
  #ylim(1985,2021) +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=10,angle=45, hjust=1),
        axis.text.y = element_text(size=12)) +
  #facet_wrap(~ Groundfish.Management.Area.Code + Sex) + 
  #facet_grid(Groundfish.Management.Area.Code ~ Sex)+	#+
  facet_wrap(~ GFMU + Sex) + 
  facet_grid(GFMU ~ Sex)+	
  scale_x_continuous(limits=c(300,900),breaks = c(seq(from=300, to=900, by=50))) +	#, labels = axisx$labels) +
  scale_y_continuous(limits=c(1984,2024),breaks = c(seq(from=1984, to=2024, by=2)))#scales::pretty_breaks(n=15)) 

ggsave("figures/EYKT_ridgejitter_lengthcomp_byyear.png", dpi=900, height=8, width=5, units="in")

##SAMPLE SIZE LABEL OPTION 
lendat %>% 
  #filter(Groundfish.Management.Area.Code == "SSEO") %>% 
  filter(GFMU == "SSEO") %>% 
  ggplot(aes(length, Year, group = Year, fill = Year)) + 
  geom_density_ridges(aes(point_fill = Year, point_color = Year),
                      alpha = 0.3, scale=2, #jittered_points = TRUE,
                      rel_min_height = 0.01) +
  xlim(300, 800) + 
  xlab("Length (mm)") + 
  ylab(NULL) +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=10,angle=45, hjust=1),
        axis.text.y = element_text(size=12)) +
  #facet_wrap(~ Groundfish.Management.Area.Code + Sex) + 
  #facet_grid(Groundfish.Management.Area.Code ~ Sex)+	#+
  facet_wrap(~ GFMU + Sex) + 
  facet_grid(GFMU ~ Sex)+	#+
  geom_label(data = lendat%>% 
             #  filter(Groundfish.Management.Area.Code == "EYKT") %>% group_by(Year), 
             filter(GFMU == "SSEO") %>% group_by(Year), 
             aes(y = Year, x = 775, 
                 label = n), colour="black", fill="white", nudge_y=0.5, size=3)+
  scale_x_continuous(limits=c(300,800),breaks = c(seq(from=300, to=800, by=50))) +	#, labels = axisx$labels) +
  #  scale_y_discrete(breaks = scales::pretty_breaks(n=15))  
  scale_y_continuous(breaks = c(seq(from=1984, to=2024, by=2)) )

ggsave("figures/SSEO_ridgelabeled_lengthcomp_byyear.png", dpi=900, height=8, width=5, units="in")

ggsave(paste0("figures/lengthcomp_ggridges_", YEAR, ".png"), 
       dpi=300, height=8, width=10, units="in")

##Bubble Age plots
lendat.bin$bin.no<-as.numeric(lendat$length_bin)
ggplot(data = lendat.bin %>% 
       #  filter(Groundfish.Management.Area.Code == "SSEO"),
       filter(GFMU == "SSEO"),
       aes(x = Year, y = bin.no, size = bin.count)) + #*FLAG* could swap size with proportion_scaled
  geom_point(shape = 21, fill = "black", colour = "black") +
  scale_size(range = c(0, 5)) +
  facet_wrap(~ Sex ) +
  labs(x = "\nYear", y = "Length Bin\n") +
  guides(size = "none")	+	#FALSE)# +
  #  theme_classic() +
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=12,angle=45, hjust=1),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=15)) +	#, labels = axisx$labels) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=15)) 

#raw lengths...str(lendat)
ggplot(data = lcomps %>% #filter(Groundfish.Management.Area.Code == "SSEO"),
       filter(GFMU == "SSEO"),
       aes(x = Year, y = length, size = proportion)) + #*FLAG* could swap size with proportion_scaled
  geom_point(shape = 21, fill = "black", colour = "black") +
  scale_size(range = c(0, 5)) +
  facet_wrap(~ Sex ) +
  labs(x = "\nYear", y = "Length\n") +
  guides(size = "none")	+	#FALSE)# +
  #  theme_classic() +
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=12,angle=45, hjust=1),
        panel.background = element_rect(fill = "white", colour = "grey50"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=15)) +	#, labels = axisx$labels) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=15))


##Janes code:
bind_rows(srv_bio %>% 
            filter(year >= 1997 &
                     Sex %in% c("Female", "Male") &
                     !is.na(length)) %>% 
            select(year, Sex, length) %>% 
            mutate(Source = "LL survey"),
          fsh_bio %>% 
            filter(year >= 2002 &
                     Sex %in% c("Female", "Male") &
                     !is.na(length)) %>% 
            select(year, Sex, length) %>% 
            mutate(Source = "LL fishery")#,
          # potsrv_bio %>% 
          #   filter(Sex %in% c("Female", "Male") &
          #            !is.na(length)) %>% 
          #   select(year, Sex, length) %>% 
          #   mutate(Source = "Pot survey")
) %>% 
  filter(!c(length < 40)) %>% 
  mutate(length2 = ifelse(length < 41, 41,
                          ifelse(length > 99, 99, length)),
         length_bin = cut(length2, breaks = seq(39.9, 99.9, 2),
                          labels = paste(seq(41, 99, 2)))) %>% 
  select(-length2) -> lendat


###done with Janes... now modifying Janes ... 
lendat %>% 
  # Length comps by Source, year, and Sex 
  count(Sex, Year, length_bin) %>%							#count(Source, Sex, year, length_bin) %>%
  group_by(Sex, Year) %>% 								#group_by(Source, Sex, year) %>%
  mutate(proportion = round( n / sum(n), 4)) %>% 
  bind_rows(lendat %>% # Sexes combined
              count(Year, length_bin) %>%						#count(Source, year, length_bin) %>%
              group_by(Year) %>% 							#group_by(Source, year) %>%
              mutate(proportion = round( n / sum(n), 4),
                     Sex.c = "Sex combined")) -> lencomps			#Sex = "Sex combined")) -> lencomps
view(lencomps)
# complete() was behaving weirdly. Expand to grid to include all length combos
expand.grid(year = unique(lencomps$Year), 
            #Source = unique(lencomps$Source),
            Sex = unique(lencomps$Sex),
            length_bin = sort(unique(lendat$length_bin)))  %>% 
  data.frame()  %>% 
  full_join(lencomps) %>%
  fill_by_value(n, proportion, value = 0) %>% 
  mutate(#length_bin = factor(length_bin),
    proportion = round(proportion, 4)) #%>%
# Keep only relevant years for each Source
#  filter(c(Source == "LL fishery" & year >= 2002) |
#           c(Source == "LL survey" & year >= 1997) #|
# c(Source == "Pot survey" & year %in% pot_yrs$year)
#        ) 
-> lencomps

# Check that they sum to 1
lencomps %>% 
  group_by(Sex, Year) %>% 						#group_by(Source, Sex, year) %>%
  summarise(sum(proportion)) #%>% View()

write_csv(lencomps, paste0("output/lengthcomps_", YEAR, ".csv"))

lendat %>% 
  # Mean length comp for comparison
  count(Source, Sex, length_bin) %>%
  group_by(Source, Sex) %>% 
  mutate(proportion = round( n / sum(n), 4)) %>% 
  arrange(Source, Sex, length_bin) %>% 
  # Fill in the blanks with 0's
  complete(Source, length_bin,
           fill = list(n = 0, proportion = 0)) %>% 
  bind_rows(lendat %>% # Sexes combined
              count(Source, length_bin) %>%
              group_by(Source) %>% 
              mutate(proportion = round( n / sum(n), 4),
                     Sex = "Sex combined")) -> mu_lencomps

mu_lencomps %>% 
  group_by(Source, Sex) %>% 
  summarise(sum(proportion))

s_lencomps <- lencomps %>% 
  filter(Source == "LL survey")

f_lencomps <- lencomps %>% 
  filter(Source == "LL fishery")

# ggridge plots

lendat %>% 
  filter(Source != "Pot survey") %>% 
  mutate(Source = derivedFactor("Survey" = Source == "LL survey",
                                "Fishery" = Source == "LL fishery",
                                .ordered = TRUE)) %>% 
  ggplot(aes(length, year, group = year, fill = year)) + 
  geom_density_ridges(aes(point_fill = year, point_color = year),
                      alpha = 0.3) +
  geom_vline(xintercept = 63, linetype = 4) +
  xlim(40, 90) + 
  xlab("\nLength (cm)") + 
  ylab(NULL) +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  facet_wrap(~ Source)

ggsave(paste0("figures/lengthcomp_ggridges_", YEAR, ".png"), 
       dpi=300, height=8, width=10, units="in")

# ggride plot for len dat by sex (for TMB inputs)

lendat %>% 
  filter(! c(Source %in% c("Pot survey", "LL fishery"))) %>% 
  mutate(Source = derivedFactor("Survey" = Source == "LL survey")) %>% 
  ggplot(aes(length, year, group = year, fill = year)) + 
  geom_density_ridges(aes(point_fill = year, point_color = year), alpha = 0.3) +
  #geom_vline(xintercept = 61, linetype = 4) + # L50
  xlim(40, 90) + 
  labs(x = "\nLength (cm)", y = "Year\n") +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  facet_wrap(~ Sex) +
  ggtitle("Survey")

ggsave(paste0("figures/tmb/lencomp_srv_",YEAR,".png"), 
       dpi=300, height=8, width=10, units="in")

# fishery
lendat %>% 
  filter(! c(Source %in% c("Pot survey", "LL survey"))) %>% 
  mutate(Source = derivedFactor("Fishery" = Source == "LL fishery")) %>% 
  ggplot(aes(length, year, group = year, fill = year)) + 
  geom_density_ridges(aes(point_fill = year, point_color = year), alpha = 0.3) +
  # geom_vline(xintercept = 61, linetype = 4) + # L50
  xlim(40, 90) + 
  labs(x = "\nLength (cm)", y = "Year\n") +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  facet_wrap(~ Sex) +
  ggtitle("Fishery")

ggsave(paste0("figures/tmb/lencomp_fsh_",YEAR,".png"), 
       dpi=300, height=8, width=10, units="in")

# All years smoothed by source
ggplot() +
  geom_point(data = lencomps %>% 
               filter(Sex == "Sex combined"),
             aes(x = length_bin, y = proportion, 
                 colour = Source),
             size = 1, alpha = 0.2) +
  # stat_smooth(size = 1.5, se = FALSE) +
  geom_line(data = mu_lencomps %>% 
              filter(Sex == "Sex combined"),
            aes(x = length_bin, y = proportion, colour = Source, 
                group = Source, linetype = Source), size = 1) +
  scale_x_discrete(breaks = seq(41, 99, 6),
                   labels = seq(41, 99, 6)) +
  scale_colour_grey() +
  xlab('\nFork length (cm)') +
  ylab('Proportion\n') +
  theme(legend.position = c(0.8, 0.8))

ggsave("figures/lengthcomp_bydatasource.png", 
       dpi=300, height=4.5, width=5, units="in")

# length comp figs, requested by AJ Lindley 2018-09-07
lencomps %>% 
  group_by(year, Source, Sex) %>% 
  dplyr::summarize(N = sum(n),
                   label = paste0("n = ", prettyNum(N, big.mark = ","))) %>% 
  ungroup() %>% 
  mutate(length_bin = "91", proportion = 0.18) -> labels 

# For survey
ggplot(data = lencomps %>% 
         # Last 10 years of data
         filter(year >= YEAR - 10 & 
                  Sex != "Sex combined" &
                  Source == "LL survey"), 
       aes(x = length_bin, y = proportion)) + 
  geom_bar(stat = "identity", colour = "lightgrey", fill = "lightgrey", width = 0.8) +
  geom_line(data = lencomps %>% 
              # Compare all past years to this year
              filter(year == YEAR & 
                       Sex != "Sex combined" &
                       Source == "LL survey") %>% 
              select(-year),
            aes(x = length_bin, y = proportion, group = 1),
            colour = "black") +
  geom_text(data = labels %>% 
              filter(year >= YEAR - 10 & 
                       Sex != "Sex combined" &
                       Source == "LL survey"),
            aes(x = length_bin, y = proportion, label = label),
            size = 3, family = "Times") +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = round(seq(0, 0.2, 0.1), 2),
                     labels =  round(seq(0, 0.2, 0.1), 2)) +
  scale_x_discrete(breaks = seq(41, 99, 6),
                   labels = seq(41, 99, 6)) +
  facet_grid(year ~ Sex) +
  labs(x = "\nFork length (cm)", y = "Proportion-at-length (longline survey)\n") +
  theme(strip.placement = "outside") 

ggsave(paste0("figures/llsrv_lencomps_", YEAR-10, "_", YEAR, ".png"), 
       dpi=300, height=8, width=6.5, units="in")

# For fishery
ggplot(data = lencomps %>% 
         # Last 10 years of data
         filter(year >= YEAR - 10 & 
                  Sex != "Sex combined" &
                  Source == "LL fishery"), 
       aes(x = length_bin, y = proportion)) + 
  geom_bar(stat = "identity", colour = "lightgrey", fill = "lightgrey", width = 0.8) +
  geom_line(data = lencomps %>% 
              # Compare all past years to this year
              filter(year == YEAR & 
                       Sex != "Sex combined" &
                       Source == "LL fishery") %>% 
              select(-year),
            aes(x = length_bin, y = proportion, group = 1),
            colour = "black") +
  geom_text(data = labels %>% 
              filter(year >= YEAR - 10 & 
                       Sex != "Sex combined" &
                       Source == "LL fishery"),
            aes(x = length_bin, y = proportion, label = label),
            size = 3, family = "Times") +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = round(seq(0, 0.2, 0.1), 2),
                     labels =  round(seq(0, 0.2, 0.1), 2)) +
  scale_x_discrete(breaks = seq(41, 99, 6),
                   labels = seq(41, 99, 6)) +
  facet_grid(year ~ Sex) +
  labs(x = "\nFork length (cm)", y = "Proportion-at-length (longline fishery)\n") +
  theme(strip.placement = "outside") 

ggsave(paste0("figures/llfsh_lencomps_", YEAR-10, "_", YEAR, ".png"), 
       dpi=300, height=8, width=6.5, units="in")

# Summary stats output for length comps (requested by AJ Linsley 20180907)
lendat %>% 
  filter(Source %in% c("LL survey", "LL fishery")) %>% 
  group_by(Source, Sex, year) %>% 
  dplyr::summarize(mean = mean(length),
                   min = min(length),
                   max = max(length)) %>% 
  mutate(variable = "Fork length") -> lensum

# axis <- tickr(lensum, year, 5)

lensum %>% 
  ggplot(aes(x = year, y = mean, colour = Source)) +
  geom_point() +
  geom_line() +
  scale_colour_grey(guide = FALSE) +
  facet_wrap(~ Sex) +
  # scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  labs(x = NULL, y = "Mean\nfork\nlength\n(cm)") +
  theme(axis.title.y = element_text(angle=0)) -> l

#quos() uses stand eval in dplyr, eval cols with nonstand eval using !!!
cols <- quos(Source, year, Sex, age) 

bind_rows(
  fsh_bio %>% mutate(Source = "LL fishery") %>% select(!!!cols), 
  srv_bio %>% mutate(Source = "LL survey") %>% select(!!!cols)) %>% 
  filter(year >= 1997 & Sex %in% c('Female', 'Male') & !is.na(age)) %>% 
  group_by(Source, Sex, year) %>% 
  dplyr::summarize(mean = mean(age),
                   min = min(age),
                   max = max(age)) %>% 
  mutate(variable = "Age") -> agesum

agesum %>% 
  mutate(age = round(mean, 0)) -> agesum1

# axisy <- tickr(agesum1, age, 3)

agesum %>% 
  ggplot(aes(x = year, y = mean, colour = Source)) +
  geom_point() +
  geom_line() +
  scale_colour_grey() +
  facet_wrap(~ Sex) +
  # scale_x_continuous(breaks = axis$breaks, labels = axis$labels) +
  # scale_y_continuous(breaks = axisy$breaks, labels = axisy$labels) +
  labs(x = NULL, y = "Mean\nage\n(yrs)") +
  theme(legend.position = "bottom",
        axis.title.y = element_text(angle=0)) -> a

cowplot::plot_grid(l, a, axis = "lrtb", align = "hv", ncol = 1) -> compare_comp_sums
compare_comp_sums
ggsave("figures/compare_comp_summaries.png",
       plot = compare_comp_sums,
       # dpi=300, height=5.5, width=6.5, units="in")
       dpi=300, height=7, width=7, units="in")

bind_rows(agesum, lensum) %>% 
  write_csv("Output/comps_summary.csv")

# survey length comps by stat area, requested by A Olson 2021-02-18
srv_bio %>% 
  filter(year >= 1997 &
           Sex %in% c("Female", "Male") &
           !is.na(length)) %>% 
  select(year, Stat, length) %>% 
  filter(!c(length < 40)) %>% 
  mutate(length2 = ifelse(length < 41, 41,
                          ifelse(length > 99, 99, length)),
         length_bin = cut(length2, breaks = seq(39.9, 99.9, 2),
                          labels = paste(seq(41, 99, 2)))) %>% 
  select(-length2) %>% 
  mutate(Stat = derivedFactor("345731" = Stat == "345731",
                              "345701" = Stat == "345701",
                              "345631" = Stat == "345631",
                              "345603" = Stat == "345603",
                              .ordered = TRUE)) -> statlen

statlen %>% 
  # Length comps by Source, year, and Sex 
  count(Stat, year, length_bin) %>%
  group_by(Stat, year) %>% 
  mutate(proportion = round( n / sum(n), 4)) -> statcomps

# complete() was behaving weirdly. Expand to grid to include all length combos
expand.grid(year = unique(statcomps$year), 
            Stat = unique(statcomps$Stat),
            length_bin = sort(unique(statlen$length_bin)))  %>% 
  data.frame()  %>% 
  full_join(statcomps) %>%
  fill_by_value(n, proportion, value = 0) %>% 
  mutate(proportion = round(proportion, 4)) -> statcomps

# Check that they sum to 1
statcomps %>% 
  group_by(Stat, year) %>% 
  summarise(sum(proportion)) 

statcomps %>% 
  group_by(year, Stat) %>% 
  dplyr::summarize(N = sum(n),
                   label = paste0("n = ", prettyNum(N, big.mark = ","))) %>% 
  ungroup() %>% 
  mutate(length_bin = "91", proportion = 0.18) -> labels 

# For survey
ggplot(data = statcomps %>% 
         # Last 10 years of data
         filter(year >= YEAR - 10), 
       aes(x = length_bin, y = proportion)) + 
  geom_bar(stat = "identity", colour = "lightgrey", fill = "lightgrey", width = 0.8) +
  geom_line(data = statcomps %>% 
              # Compare all past years to this year
              filter(year == YEAR) %>% 
              select(-year),
            aes(x = length_bin, y = proportion, group = 1),
            colour = "black") +
  geom_text(data = labels %>% 
              filter(year >= YEAR - 10),
            aes(x = length_bin, y = proportion, label = label),
            size = 3, family = "Times") +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = round(seq(0, 0.2, 0.1), 2),
                     labels =  round(seq(0, 0.2, 0.1), 2)) +
  scale_x_discrete(breaks = seq(41, 99, 6),
                   labels = seq(41, 99, 6)) +
  facet_grid(year ~ Stat) +
  labs(x = "\nFork length (cm)", y = "Proportion-at-length (longline survey)\n") +
  theme(strip.placement = "outside") 

# ggridge plots

statlen %>% 
  filter(year >= YEAR - 5) %>% 
  ggplot(aes(length, year, group = year, fill = year)) + 
  geom_density_ridges(aes(point_fill = year, point_color = year),
                      alpha = 0.3) +
  xlim(40, 90) + 
  xlab("\nLength (cm)") + 
  ylab(NULL) +
  # scale_y_reverse() +
  theme(legend.position = "none") + 
  facet_wrap(~ Stat, ncol = 1)

statlen %>%
  filter(year >= YEAR - 10 & between(length, 40, 99)) %>%
  ggplot(aes(x = length, y = factor(year), group = interaction(year, Stat),
             fill = Stat)) + #, col = Stat)) +
  geom_density_ridges(alpha = 0.5, size = 0.5, col = "white") +
  # scale_fill_manual(values = c("grey30", "#00BFC4", "#da2055", "#daa520" )) +
  # scale_colour_manual(values = c("grey30", "#00BFC4", "#daa520", "#da2055")) +
  labs(x = "Length (cm)", y = NULL, fill = NULL, col = NULL) + #, title = "Dusky rockfish") +
  theme_light() +
  scale_x_continuous(limits = c(40, 99)) +
  theme(legend.position = "top")

ggsave(paste0("figures/lengthcomp_statarea_", YEAR, ".png"), 
       dpi=300, height=8, width=10, units="in")
