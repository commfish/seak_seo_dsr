#############################################################################
## Bubble plots for 2021 DSR SAFE report
## Sept 2021
## Phil Joy
#############################################################################

library(tidyverse)
library(ggridges)
library(ggplot2)
source("Code/Port_bio_function.R")

#setwd("D:/Groundfish Biometrics/Yelloweye/SAFE reports")

{ Port1<-read.csv("Data/SEO_YE_port_sampling_bio_data_1980-1989.csv")
  Port2<-read.csv("Data/SEO_YE_port_sampling_bio_data_1990-1999.csv")
  Port3<-read.csv("Data/SEO_YE_port_sampling_bio_data_2000-2009.csv")
  Port4<-read.csv("Data/SEO_YE_port_sampling_bio_data_2010-2019.csv")
  Port5<-read.csv("Data/SEO_YE_port_sampling_bio_data_2020-2022.csv")
  
  Port<-rbind(Port1,Port2,Port3, Port4, Port5); nrow(Port)
  #str(Port)
  Port$Year<-as.integer(Port[,1])
  #str(Port)
    #EYAK = EYKT
  unique(Port$Groundfish.Stat.Area)
  unique(Port$Groundfish.Stat.Area.Group)
  unique(Port$Groundfish.Management.Area.Code)
  
  Port$Groundfish.Management.Area.Code[Port$Groundfish.Management.Area.Code == "EYAK"]<-"EYKT"
  
  statareas<-read.csv("Data/g_stat_area.csv")
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

str(Port)

Port<-port.bio(2022)
Port.rand<-Port[Port$Sample.Type=="Random",]

agecomps<-Port[!is.na(Port$Sex) & !is.na(Port$Age),]
nrow(agecomps)
agecomps.rand<-Port.rand[!is.na(Port.rand$Sex) & !is.na(Port.rand$Age),]
nrow(agecomps.rand)

agecomps<-agecomps.rand
#agecomps<-agecomps[agecomps$Project == "Halibut Longline" | agecomps$Project == "Commercial Longline Trip",]
unique(agecomps$Groundfish.Management.Area.Code)
unique(agecomps$GFMU)
agecomps<-agecomps[agecomps$GFMU == "NSEO" |
			agecomps$GFMU == "SSEO" |
			agecomps$GFMU == "CSEO" |
			agecomps$GFMU == "EYKT",]
nrow(agecomps)
hist(agecomps$Age)

ggplot(agecomps, aes(x=Age)) +
  geom_histogram(binwidth=10) +
  geom_vline(xintercept=50,col="blue") +
  geom_vline(xintercept=75,col="blue", linetype=2) +
  facet_wrap(~Year)

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

hist(lcomps$Length.Millimeters)
min(lcomps$Length.Millimeters)
max(lcomps$Length.Millimeters)

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
  mutate(bin.n=n()) %>%
  #group_by(length_bin,Year,Groundfish.Management.Area.Code, Sex) %>%
  group_by(length_bin,Year,GFMU, Sex) %>%
  mutate(bin.count=n()) %>%
  mutate(bin.proportion = bin.count/bin.n) 


tst<-lendat %>% filter(GFMU == "CSEO")
unique(tst$length)
unique(tst$Year)
unique(tst$GFMU)

###########################################################################
plus_group<-100
###############################################################################
## Autoplotting... skip this and go to INDIVIDUAL PLOTS to make sure your
## shit works... 

unique(Port$GFMU)
gmus<-unique(Port$GFMU)[-5]
YEAR<-2022

for (i in gmus){  #i<-gmus[1]
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

    ggsave(paste0("Figures/",YEAR,"Safe/",i,"_bubble_agecomp_byyear.png",sep=""),
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
      scale_y_continuous(limits=c(1984,YEAR),breaks = c(seq(from=1984, to=(YEAR+2), by=2)))#scales::pretty_breaks(n=15)) 
    
    ggsave(paste0("Figures/",YEAR,"Safe/",i,"_ridgejitter_lengthcomp_byyear.png",sep=""),
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
    
    ggsave(paste0("Figures/",YEAR,"Safe/",i,"_ridgelabeled_lengthcomp_byyear.png",sep=""),
           dpi=900, height=8, width=5, units="in")
}

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
  scale_y_continuous(limits=c(1984,2022),breaks = c(seq(from=1984, to=2024, by=2)))#scales::pretty_breaks(n=15)) 

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
  scale_y_continuous(breaks = c(seq(from=1984, to=2022, by=2)) )

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
  write_csv("output/comps_summary.csv")

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
