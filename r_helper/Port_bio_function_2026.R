## Function for loading port sampling bio data
## I rewrote Phil's function - I think this still achieves what he was trying to do
## I included inside waters and all projects so this function can be easily used for all projects
## Laura Coleman Last Edited: 7/21/26

port.bio<-function(YEAR=2026){
  
  #Read in biodata from 1981 - 2010
  Port1 <- read.csv("Data_processing/Data/SEAK_YE_port_sampling_bio_data_1981-2010.csv")
  
  #Pull new data from Region 1/GroundFish/User Reports/Yelloweye Reports for Phil/port sampling bio data
  Port2 <- read.csv(paste0("Data_processing/Data/SEAK_YE_port_sampling_bio_data_2011-", YEAR, ".csv"),na.strings = c("", " ", "NA"))
  
  #combine Files
  Port<-rbind(Port1,Port2) %>% 
    filter(!Groundfish.Management.Area.Code=="WYAK")
  
  #Data Checks -----------------------------------------------------------------
  unique(Port$Groundfish.Management.Area.Code)
  
  Port$Year<-as.integer(Port[,1])
  
  #Data Cleaning ---------------------------------------------------------------
  
  #Read in spreadsheet with groundfish stat areas
  statareas<-read.csv("Data_processing/Data/g_stat_area.csv")
  
  #Build lists of stat area codes belonging to each groundfish management area
  #outside waters
  SSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEO"])
  CSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "CSEO"])
  NSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEO"])
  EYc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "EYKT"])
  
  #inside waters
  SSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEI"])
  NSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEI"])
  
  #Create GFMU column using the existing groundfish management area if present or will use
  #the stat area or stat area group
  Port <- Port %>%
    mutate(GFMU = case_when(
      
      Groundfish.Management.Area.Code != "" ~ Groundfish.Management.Area.Code,
      
      Groundfish.Stat.Area %in% SSc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "SSEO" ~ "SSEO",
      Groundfish.Stat.Area %in% CSc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "CSEO" ~ "CSEO",
      Groundfish.Stat.Area %in% NSc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "NSEO" ~ "NSEO",
      Groundfish.Stat.Area %in% EYc  | substr(Groundfish.Stat.Area.Group, 1, 4) == "EYKT" ~ "EYKT",
      Groundfish.Stat.Area %in% SSIc | substr(Groundfish.Stat.Area.Group, 1, 4) == "SSEI" ~ "SSEI",
      Groundfish.Stat.Area %in% NSIc | substr(Groundfish.Stat.Area.Group, 1, 4) == "NSEI" ~ "NSEI",TRUE ~ NA_character_)) %>% 
    
    #Create sex column to translate numeric code to words
    mutate(Sex = case_when(Sex.Code == 1 ~ "Male", Sex.Code == 2 ~ "Female")) 
  
  Port$Sex<-as.factor(Port$Sex)	
  
  return(Port)
}

