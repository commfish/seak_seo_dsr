################################################################################
## Spatial changes in the halibut harvest over time
## Goal is to identify shift associated with the change from derby to IFQ in 1995
## Phil Joy, May 2022
################################################################################

Hal<-read.csv("Data/Halibut_landings_tns_1982_2021.csv", header=T, skip=3)
str(Hal)

#Load halibut harvests by stat area data
Hal %>% mutate(Year = Year1, Reg.Area=as.factor(IPHC.Regulatory.Area),
               Stat.Area=as.factor(IPHC..Statistical.Area),
               Net.wt=as.numeric(Pacific.halibut.Net.wt..t.),
               cpue=Net.wt/Unique.vessels) %>%
  select(Year,Reg.Area,Stat.Area,Net.wt,Unique.vessels,cpue) ->Hal

head(Hal)  

#load IPHC survey data to get stat.areas and SE district conversion 
{
  BUT2C<-read.csv("Data/IPHC Set and Pacific halibut data 2C.csv", header=T)
  YE2C<-read.csv("Data/Non-Pacific halibut data YE 2C.csv", header=T)
  Sur2C<- BUT2C %>% full_join(YE2C, by="Stlkey")
  Sur2C<- Sur2C %>%
    mutate(SEdist = ifelse(IPHC.Stat.Area %in% c(182:184,171,173,174,161:163),"NSEI", 
                           ifelse(IPHC.Stat.Area %in% c(142:144,152,153),"SSEI",
                                  ifelse(IPHC.Stat.Area %in% c(140,141,151,150,160,170,181),
                                         ifelse(-MidLon.fished>137,"EYKT",
                                                ifelse(MidLat.fished>=57.5,"NSEO",
                                                       ifelse(MidLat.fished<57.5 & MidLat.fished>=56,"CSEO",
                                                              ifelse(MidLat.fished<56,"SSEO",NA)))),NA))))
  
  ## Lets bring in 3A and pull out surveys in the Yakutat area...
  BUT3A<-read.csv("Data/IPHC Set and Pacific halibut data 3A.csv", header=T)
  YE3A<-read.csv("Data/Non-Pacific halibut data YE 3A.csv", header=T)
  Sur3A<- BUT3A %>% full_join(YE3A, by="Stlkey")
  Sur3A<-Sur3A %>%
    mutate(SEdist = ifelse(-MidLon.fished <= 137 & MidLat.fished >= 57.5,"NSEO",
                           ifelse(-MidLon.fished > 137 & -MidLon.fished <139 & MidLat.fished > 54.5,  #Lon should be less than 138 but included here to increase sample size of the index... 
                                  "EYKT","whocares")))
  Sur3A<-Sur3A[Sur3A$SEdist != "whocares",]
  
  ## connect the two data sets
  Survey<-rbind(Sur3A,Sur2C)
  #make one more column to define inside versus outside for dealing with older data sets...
  Survey<-Survey %>% 
    mutate(In.Out = ifelse(SEdist %in% c("NSEO","CSEO","SSEO","EYKT"),"SEO",
                           ifelse(SEdist %in% c("NSEI","SSEI"),"SEI",NA))) %>%
    rename(YE.obs = Number.Observed,
           AvgDepth.fm = AvgDepth..fm.,
           Year = Year.x)%>%
    mutate(depth_bin = cut(AvgDepth.fm, breaks = seq(0,400,50),
                           labels = paste (seq(50,400,50))),
           YE.obs = ifelse(is.na(YE.obs),0,YE.obs),
           HooksObserved  = ifelse(is.na(HooksObserved ),140,HooksObserved ))
  
  Survey$HooksRetrieved<-as.numeric(Survey$HooksRetrieved)
  Survey$YE.exp<-Survey$HooksRetrieved/Survey$HooksObserved
  
  Survey<-Survey[order(Survey$Year),]
}

#check to see which stat areas occur in more than one SE district
str(Survey)
check<-data.frame()
sas<-unique(Survey$IPHC.Stat.Area)
i<-1
for (s in sas){  #s<-sas[2]
  dat<-Survey[Survey$IPHC.Stat.Area == s,]
  
  if (length(unique(dat$SEdist)) == 1) {
    check[i,"stat.area"]<-s
    check[i,"district"]<-as.list(unique(dat$SEdist))
    i<-i+1
  } else {
    for (j in 1:length(unique(dat$SEdist))){
      check[i,"stat.area"]<-s
      check[i,"district"]<-unique(dat$SEdist)[j]
      i<-i+1
    }
  }
}
check; length(sas)

## Get southeast halibut harvests 
#calculate deviation from the mean harvest for each stat area in each year
Hal %>% filter(Stat.Area %in% sas) %>% 
  group_by(Year) %>% 
  mutate(catch_u = mean(Net.wt),
         cpue_u = mean(cpue),
         ttl.vessels = sum(Unique.vessels)) %>% 
  group_by(Stat.Area) %>% 
  mutate(catch_dev = Net.wt-catch_u,
         cpue_dev = cpue-cpue_u,
         prop.vessels = Unique.vessels/ttl.vessels)-> SE.Hal

head(SE.Hal[SE.Hal$Year==2021,], 30)
head(SE.Hal[SE.Hal$Year==2020,], 30)

#SE.Hal %>% group_by(Stat.Area) %>%
ggplot(SE.Hal,aes(Year,catch_dev, colour=Stat.Area)) + 
  geom_line(size=1)

#catch in area relative to average catch of all areas for that year
ggplot(SE.Hal,aes(Year,catch_dev)) + 
  geom_line(size=1) +
  geom_vline(xintercept=1995,color="blue")+
  geom_hline(yintercept=0,color="darkgrey") +
  facet_wrap(~Stat.Area)

#proportion of all vessels fishing that area
ggplot(SE.Hal,aes(Year,prop.vessels)) + 
  geom_line(size=1) +
  geom_vline(xintercept=1995,color="blue")+
  facet_wrap(~Stat.Area)

#cpue in area relative to mean cpue across all areas for that year
ggplot(SE.Hal,aes(Year,cpue_dev)) + 
  geom_line(size=1) +
  geom_vline(xintercept=1995,color="blue")+
  geom_hline(yintercept=0,color="darkgrey") +
  facet_wrap(~Stat.Area)

#calculate average deviation for IFQ and DERBY for easier metrics... 
