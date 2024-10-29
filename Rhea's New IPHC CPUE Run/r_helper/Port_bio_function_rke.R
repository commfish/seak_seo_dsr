## Load port sampling bio data

port.bio<-function(YEAR=2024){

Port1<-read.csv("data/SEO_YE_port_sampling_bio_data_1980-1989.csv")
Port2<-read.csv("data/SEO_YE_port_sampling_bio_data_1990-1999.csv")
Port3<-read.csv("data/SEO_YE_port_sampling_bio_data_2000-2009.csv")
Port4<-read.csv("data/SEO_YE_port_sampling_bio_data_2010-2019.csv")
Port5<-read.csv("data/SEO_YE_port_sampling_bio_data_2020-2024.csv")

Port<-rbind(Port1,Port2,Port3, Port4, Port5)
#str(Port)
Port$Year<-as.integer(Port[,1])
#str(Port)
#EYAK = EYKT
Port$Groundfish.Management.Area.Code[Port$Groundfish.Management.Area.Code == "EYAK"]<-"EYKT"
unique(Port$Groundfish.Management.Area.Code)
statareas<-read.csv("data/g_stat_area.csv")
#str(statareas)
unique(statareas$G_MANAGEMENT_AREA_CODE)

SSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEO"])
CSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "CSEO"])
NSc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEO"])
EYc<-unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "EYKT"])

SSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "SSEI"])
NSIc <- unique(statareas$G_STAT_AREA[statareas$G_MANAGEMENT_AREA_CODE == "NSEI"])

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



Port %>% mutate(Sex = 
                        case_when(Sex.Code == 1 ~ "Male",
                                  Sex.Code == 2 ~ "Female") -> Port
)
Port$Sex<-as.factor(Port$Sex)	

return(Port)
}

str(Port)
with(Port, table(Year, GFMU))

