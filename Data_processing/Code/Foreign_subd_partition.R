#################################################################################
## Historical Foreign Fleet YE removals in SEO
## Data from Donnie Arthuer SFD June 2022
## Resolution at the SEO level
## For SS-SPM need to partition into subdistrict level using historical catch 
## proportions in historical times
## Modelled in SS-SPM with very large cv's (over 1)
#################################################################################

# 1) Historical catch proportions

Removals.subd<-read.csv("Data/SEsubdistrict_YE_removals.csv")  #YE reconstruction from Rhea
str(Removals.subd); unique(Removals.subd$Year)

#2) pre-1980 removals from SEO as a whole... 
# Removals.pre<-read.csv("Data/XXX.csv") #Waiting on Donnie; 
# generate fake data here for model development...
Tots<-Removals.subd %>% filter(Subdistrict != "NSEI" &
                                 Subdistrict != "SSEI" &
                                 Subdistrict != "SEO") %>%
  group_by(Year) %>% 
  mutate(SEO.harv = sum(Remove.mt)) %>%
  ungroup() %>%
  mutate(subd.prop = Remove.mt/SEO.harv) 


ggplot(data=Tots,aes(Year,subd.prop, color=as.factor(Subdistrict))) + geom_line() 
  
Tots %>% group_by(Subdistrict) %>%
  summarize(mean.prop = mean(subd.prop, na.rm=T),
         sd.prop = sd(subd.prop, na.rm=T)) -> props

Fish.cpue<-read.csv("Data/YE_dirfishery_nom_CPUE.csv")
min(Fish.cpue$Year); max(Fish.cpue$Year)
N.fish<-length(unique(Fish.cpue$Year))
N.fish.closed<-YEAR-max(Fish.cpue$Year)

#2) Load Donnie's foreign fleet stuff
Foreign<-read.csv("Data/Foreign_YERF_SEAK.csv", skip=1)
str(Foreign)

Years<-unique(Foreign$Year)
subs<-unique(Tots$Subdistrict)

Fsub<-data.frame()
i<-1

for (y in Years) {
  for (s in subs) {
    Fsub[i,"Year"]<-y
    Fsub[i,"Area"]<-s
    Fsub[i,"mean.prop"]<-props$mean.prop[props$Subdistrict == s]
    Fsub[i,"sd.prop"]<-props$sd.prop[props$Subdistrict == s]
    Fsub[i,"sub.remove.mt"]<-Fsub[i,"mean.prop"]*Foreign$Estimate[Foreign$Year == y]
    i<-i+1
  }
}

Fsub<-Fsub[Fsub$Year < 1982,]

write.csv(Fsub,file="Data/Foreign_removals_1960-1982.csv")



















 
 