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
source("r_helper/Port_bio_function_2026.R")
}

#Biological data pulled on 6/23/25 - LSC

###############################################################################
# Run Port Bio Function ----
YEAR = 2026
Port<-port.bio(YEAR) 

###############################################################################
# DATA EXPLORATION ----

unique(Port$Year)

#How many samples are missing the GFMU?
unknowns <- Port %>%
  filter(is.na(GFMU))

knowns <- Port %>%
  filter(!is.na(GFMU))

# Unique GFMU values in known data
unique(knowns$GFMU)

# Counts
nrow(unknowns)
nrow(knowns)

# Unique values within unknowns
unique(unknowns$Groundfish.Management.Area.Code)
unique(unknowns$Groundfish.Stat.Area.Group)
unique(unknowns$Groundfish.Stat.Area)

#Let's only look at longline samples from the DSR management areas in SEO with known GFMU
Port <- Port %>% 
  filter(GFMU %in% c("NSEO", "SSEO", "CSEO", "EYKT", "SEO"),
         !Project %in% c("Commercial Jig Trip","Atypical Longline Sample","Atypical Jig Sample"))

unique(Port$Project)

# How many samples do we have per sample type?
Port.rand <- Port %>%
  filter(Sample.Type == "Random")

unique(Port.rand$Year) #1981,1984-1985, 1987-2005, 2008-2026
unique(Port.rand$Project) 
nrow(Port.rand) #61,076

Port.rand %>% 
  group_by(Year) %>% 
  summarise(n=n())

Port.direct <- Port %>%
  filter(Project == "Commercial Longline Trip")

unique(Port.direct$Year) #1981, 1984-1985, 1987-2005, 2008-2019, 2025-2026
nrow(Port.direct) #49,500

Port.hal <- Port %>%
  filter(Project == "Commercial Halibut Longline")
unique(Port.hal$Year) #2003, 2008-2012, 2014-2026
nrow(Port.hal) #12,014


agecomps<-Port[!is.na(Port$Sex) & !is.na(Port$Age),]
nrow(agecomps) #37,042
unique(Port$Sample.Type)
agecomps.rand<-Port.rand[!is.na(Port.rand$Sex) & !is.na(Port.rand$Age),]
nrow(agecomps.rand) #36,725

agecomps<-agecomps.rand
unique(agecomps$Groundfish.Management.Area.Code)
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
gmus <- unique(Port$GFMU)

unique(agecomps$Sample.Type)
unique(lendat$Sample.Type)


for (i in gmus){  
  # i<-gmus[5]
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
      filter(GFMU == i & Year > 1983) %>% #change to 1980 for CSEO
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
      filter(GFMU == i & Year > 1983) %>% #change to 1980 for CSEO
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
      scale_y_continuous(breaks = c(seq(from=1983, to=YEAR, by=2)) )
    
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

ggsave(paste0("Figures/Bio_plots_",YEAR,"/","SEAK","_ridgelabeled_lengthcomp_byyear.png",sep=""),
       dpi=900, height=8, width=5, units="in")


