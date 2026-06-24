## Function for loading port sampling bio data
## I rewrote Phil's function - I think this still achieves what he was trying to do
## I included inside waters so this function can be easily used for the inside assessment

port.bio<-function(YEAR=2026){
  
  Port1 <- read.csv("Data_processing/Data/SEAK_YE_port_sampling_bio_data_1981-2010.csv")
  
  unique(Port1$Year)
  
  #Pull new data from Region 1/GroundFish/User Reports/Yelloweye Reports for Phil/port sampling bio data
  #Filters applied: species -- 145, year >=2011
  Port2 <- read.csv(paste0("Data_processing/Data/SEAK_YE_port_sampling_bio_data_2011-", YEAR, ".csv"),na.strings = c("", " ", "NA"))
  
  Port<-rbind(Port1,Port2)
  Port$Year<-as.integer(Port[,1])
  
  statareas<-read.csv("Data_processing/Data/g_stat_area.csv")
  
  SSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEO"])
  CSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "CSEO"])
  NSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEO"])
  EYc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "EYKT"])
  
  SSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEI"])
  NSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEI"])
  
  Port <- Port %>%
    mutate(GFMU = case_when(
      Groundfish.Management.Area.Code != "" ~ Groundfish.Management.Area.Code,
      Groundfish.Stat.Area %in% SSc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "SSEO" ~ "SSEO",
      Groundfish.Stat.Area %in% CSc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "CSEO" ~ "CSEO",
      Groundfish.Stat.Area %in% NSc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "NSEO" ~ "NSEO",
      Groundfish.Stat.Area %in% EYc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "EYKT" ~ "EYKT",
      Groundfish.Stat.Area %in% SSIc | substr(Groundfish.Stat.Area.Group, 1, 4) == "SSEI" ~ "SSEI",
      Groundfish.Stat.Area %in% NSIc | substr(Groundfish.Stat.Area.Group, 1, 4) == "NSEI" ~ "NSEI",TRUE ~ NA_character_)) %>% mutate(Sex = 
                                                                                                                                       case_when(Sex.Code == 1 ~ "Male",
                                                                                                                                                 Sex.Code == 2 ~ "Female")) %>% 
    mutate(GFMU = case_when(GFMU == "SSEOC"~"SSEO", TRUE~GFMU))
  
  Port$Sex<-as.factor(Port$Sex)	
  
  return(Port)
}
