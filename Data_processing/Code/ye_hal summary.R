# Summarizing halibut and yelloweye rockfish fish ticket data for Phil Joy
# Includes: CFEC Gross Earnings subject area in OceanAK 
# https://oceanak.adfg.alaska.gov/analytics/saw.dll?Answers&path=%2Fshared%2FCommercial%20Fisheries%2FRegion%20I%2FGroundFish%2FUser%20Reports%2FYelloweye%20Reports%20for%20Phil%2FCFEC%20Gross%20Earnings#resultsTab197eb6751ad
# Author:  Rhea Ehresmann (rhea.ehresmann@alaska.gov) 
# Last modified: July 8,2025 by LSC

# set up ----
source('Data_processing/Code/helper.r') 

###  set plotting theme to use TNR  ###
#font_import() #remove # to run this but only do this one time - it takes a while
loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_bw(base_size=16,base_family='Times New Roman')
          +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


##########################################################################################
### IMPORT DATA ###
##########################################################################################


#Load stat areas
SAcode<- read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/SE_Stat_areas.csv")
unique(SAcode$State.MU)
SSEIcode<-SAcode$Stat.Area[SAcode$State.MU == "SSEI"]
NSEIcode<-SAcode$Stat.Area[SAcode$State.MU == "NSEI"]
SSEOcode<-SAcode$Stat.Area[SAcode$State.MU == "SSEO"]
CSEOcode<-SAcode$Stat.Area[SAcode$State.MU == "CSEO"]
NSEOcode<-SAcode$Stat.Area[SAcode$State.MU == "NSEO"]
EYKTcode<-SAcode$Stat.Area[SAcode$State.MU == "EYKT"]
Code2C<-SAcode$Stat.Area[SAcode$IFQ.area == "2C"]
Code3A<-SAcode$Stat.Area[SAcode$IFQ.area == "3A"]


# Data from OceanAK query in CFEC Gross earnings for all YE and Halibut data 

D1<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/1975-1981.csv", fileEncoding = 'UTF-8-BOM') #1975-1981
D2<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/1982-1987.csv", fileEncoding = 'UTF-8-BOM') 
D3<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/1988-1992.csv", fileEncoding = 'UTF-8-BOM') #1988-1992 
D4<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/1993-1995.csv", fileEncoding = 'UTF-8-BOM') #1993-1995 
D5<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/1996-1997.csv", fileEncoding = 'UTF-8-BOM') #1996-1998
D6<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/1998-2000.csv", fileEncoding = 'UTF-8-BOM') %>% #1999-2000
  filter(Year.Landed!=1998)#1998 is included in D5, so this was creating issues
D7<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2001-2002.csv", fileEncoding = 'UTF-8-BOM') #2001-2002
D8<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2003-2004.csv", fileEncoding = 'UTF-8-BOM') #2003-2004
D9<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2005-2006.csv", fileEncoding = 'UTF-8-BOM') #2005-2006
D10<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2007-2008.csv", fileEncoding = 'UTF-8-BOM') #2007-2008
D11<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2009-2010.csv", fileEncoding = 'UTF-8-BOM') #2009-2010
D12<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2011-2012.csv", fileEncoding = 'UTF-8-BOM') #2011-2012
D13<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2013-2014.csv", fileEncoding = 'UTF-8-BOM') #2013-2014
D14<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2015-2016.csv", fileEncoding = 'UTF-8-BOM') #2015-2016
D15<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2017-2018.csv", fileEncoding = 'UTF-8-BOM') #2017-2018
D16<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2019-2020.csv", fileEncoding = 'UTF-8-BOM') #2019-2020
# D17<-read.csv("Data_processing/Data/Harvests/CFEC Gross Earnings Data/2021-2023.csv", fileEncoding = 'UTF-8-BOM') %>%  #2021-2023
# select("Year.Landed","Fish.Ticket.Number","Port.Name","Permit.Fishery","CFEC.Permit.Fishery",
#        "Species.Code","IPHC.Regulatory.Area","IPHC.Statistical.Area","ADFG.Management.Area.Code",
#        "Gear.Code","Gear.Description","Harvest.Code","Harvest.Description","Delivery.Code",
#        "Delivery.Description","Disposition.Code","Disposition.Description","Pounds","Whole.Pounds",
#        "CFEC.Whole.Pounds","Pre.Print.Ticket","Fish.Ticket.Number.1","Office.Name","IFQ.Halibut.Area",
#        "IFQ.Sable.Area","Groundfish.Mgt.Area.District","NMFS.Area","Stat.Area","Statistical.Area",
#        "Weight.Modifier","CFEC.Permit.Type","CFEC.Permit.Holder.Name","CFEC.Permit.Serial.Number",
#        "Port.Name.1","Date.Landed","Date.Fishing.Began","Date.Fishing.Ended","Vessel.ADFG.Number") %>% 
#   rename(Fish.TIcket.Number=Fish.Ticket.Number,
#          Fish.TIcket.Number.1=Fish.Ticket.Number.1)
# I realized that this output is not complete and there must be other filter I need to apply in OceanAK

# combine into one df
ye_hal<-rbind(D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,D16) 

# organize the df 
ye_hal <- ye_hal %>% 
  rename_all(tolower)


# filter out the fisheries and areas we definitely don't want, as best as we can broadly
# do not use iphc.regulatory.area -- WRONG AREA DATA
ye_hal2 <- ye_hal %>% 
  mutate(stat.area.new = ifelse(is.na(stat.area), statistical.area, stat.area), 
         stat.area.new2 = ifelse(stat.area.new %in% c(0, "000000", 999999, "NULL00", NA, ""), iphc.statistical.area, stat.area.new)) %>% 
  filter(!ifq.halibut.area %in% c("3B", "4A", "4B", "4C", "4D", "4E", "CL"), 
         !groundfish.mgt.area.district %in% c("WYAK",  "LOWCI", "PWS" , "PWSI", "MSAPE", "NGFW",  "PWSF",  "AISD",  "NGSW",  "NORTE",
                                              "AFOGK", "ESIDE", "SWEST", "MSAPW", "IBS", "LCHIG", "PWSE",  "PWSW", "TBBY",  "BSEA",  
                                              "CISW",  "MDLAN", "WSIDE", "CIFW",  "DNUT",  "MITRO", "CHIGB", "SUTWK", "INTL", "SEAST",
                                              "KOD",   "ZNS",   "CHUK"), 
         !stat.area.new2 %in% 400000:900000, 
         #!iphc.statistical.area %in% 210:500, 
         !iphc.statistical.area %in% 400000:900000,
         !stat.area.new2 %in% c("24140",  "25966",  "28225",  "25284",  "25238",  "25263",  "25267",  "25279",  "25282",  "30310",  "36281", "1065", 
                       "2065", "3065", "18310", "18140", "22150", "24420", "24590", "25231", "25233", "25252", "25898", "25910", "25967", "27220", 
                       '012930', "51", "517", "52", "540", "610", "620", "630", "640", "19110",	"20102",	"20300",	"20305",	"22140",	"22521",
                       '24460',	'245500',	'25120',	'25220',	'25230',	'25235',	'25254',	'25262',	'25283',	'25311',	'25331',	'25333',	'25410',
                       '255601',	'25640',	"25783", "25810",	"25855",	"25891",	"25892",	"25897",	"25961",	"27380",	"29122",	"29131",	
                       "29132",	"29142",	"29152",	"29153", "30215",	"30216",	"30217", "30221",	"30222",	"30230",	"30231",	"30240", 
                       "32422",	"32423","32621",	"32670", "18300", "35631", "4065", "5065", "680", "0", "305902", "259010", "NA"), 
         # gear.description %in% c("Longline", "Other/unspecified/missing"), #LSC removed this filter so all data was included in "Halibut harvest reconstruction 2025"
         !harvest.code %in% c(43, 44, 42)) %>%  # remove test fisheries  
         #port.name != "Metlakatla") %>%  
         #!permit.fishery %in% c("", "9998", "0000")) %>% 
  mutate(mgmt.area = ifelse(stat.area.new2 %in% c(10100,  10111,  10115,  10120,  10121,  10123,  10125,  10127,  10128,  10129,
                                            10141,  10143,  10144,  10145,  10146,  10147,  10153,  10180,  10185,  10190, 
                                            10200,  10210,  10220,  10230,  10240,  10250,  10260,  10270,  10280,  10300,
                                            10311,  10321,  10325,  10330,  10340,  10350,  10360,  10370,  10380,  10390,
                                            10500,  10510,  10520,  10541,  10550,  10600,  10610,  10620,  10630,  10641,  10642,
                                            10700,  10710,  10720,  10730,  10740,  10800,  10810,  10820,  10830,  10840, 10175, 10860, 
                                            10850, 10195, 10140, 10365, 10110, 305431, 305501, 305502, 305503, 305531, 305532, 315401, 
                                            315431, 315432, 315501, 315502, 315503, 315504, 315531, 315532, 315600, 315630, 325401, 325431, 
                                            325433, 325501, 325502, 325503, 325504, 325531, 325532, 325533, 325600, 325601, 325602, 325603, 
                                            325604, 325621, 325631, 325632, 335506, 335533, 335534, 335535, 335601, 335602, 335603,
                                            335632, 335633, SSEIcode, 142, 143, 144, 152, 153, 145), "SSEI", 
                            ifelse(stat.area.new2 %in% c(10900,  10910,  10920,  10930,  10950,  10951,  10961,  11000,  11015,  
                                                   11031,  11150,  11200,  11211,  11212,  11213,  11215,  11216,  11217,  
                                                   11218,  11219,  11221,  11280,  11400,  11421,  11423,  11500,  11531, 11425,  
                                                   11534, 11100, 11114, 10942, 11214, 11250, 11427, 11245, 11024, 10962, 11017, 11016, 
                                                   11431, 11222, 11241, 10911, 11510, 11013, 11470, 11450, 11355, 10963, 11011, 11460, 
                                                   325700, 335631, 335634, 335701, 335702, 335703, 335704, 335705, 335731, 335732, 335733, 
                                                   335734, 335735, 335830, 345534, 345535,  345603, 345604, 345605, 
                                                   345606, 345630, 345631, 345632, 345701, 345702, 345703, 345704, 345705, 
                                                   345706, 345731, 345732, 345801, 345802, 345803, 345830, 355707, 355731, 355732, 355733, 171, 172, 173, 174,
                                                   355801, 355802, 355830, 355900, 365804, 365830, 365831, NSEIcode, 161, 162, 163, 182, 183, 184), "NSEI", 
                                   ifelse(stat.area.new2 %in% c(10400,  10410,  10420,  10430, 10435, 10440,  10450, 15000, 15200, 325402, 325434, 325432, 345601,
                                                           335433, 335504, 335505, 335532, 335531, 345532, 345536, 345537, 345500, 345430, 335431, 335432,
                                                           335401, 335503, 345533, 345401, 355401, 355430, 355500, 355530, 345533, 335501, 335502, 345531, 
                                                           365430, 365500, 365530, 365434, SSEOcode, 140, 141, 150, 151), "SSEO", 
                                          ifelse(stat.area.new2 %in% c(11300,  11310,  11311,  11312,  11321,  11322,  11331, 11341, 11332, 1540,   
                                                                 11343,  11345,  11351,  11357,  11361,  11362,  11365, 15400, 11301, 11335, "01540.", 
                                                                 11337, 11342, 11330, 11338, 11334, 345602, 345607, 345608, 1540, 01540, 01540., 
                                                                 355601, 355602, 355631, 355632, 355633, 355634, 355635, 355701, 355702, 355703, 355704, 
                                                                 355705, 355706, 365600, 365630, 365701, 365702, CSEOcode, 160, 170), "CSEO", 
                                                 ifelse(stat.area.new2 %in% c(11371,  11381,  11391,  11395, 11600, 11605,  11611, 11396, 11393, 
                                                                         11394, 11397, 11392, 365731, 365732, 365733, 365734, 365735, 365801,
                                                                         365802, 365803, 365805, NSEOcode, 181, 185), "NSEO", 
                                                        ifelse(stat.area.new2 %in% c(11612,  11614,  11620,  11625,  15700,  18100, 18300, 1570, 01570, 01570., "01570." ,
                                                                               18120, 18190,  18900, 18270, 15600, 18160, 18110, 375532, 
                                                                               11613, 375530, 375600, 375630, 375700, 375730, 375801, 375802, 375831, 375832,
                                                                               375901, 377800, 385600, 385700, 385730, 385800, 385831, 385901, 385902, 395630,
                                                                               395700, 395730, 395800, 395830, 395900, 395901, 395902, 395500, 
                                                                               EYKTcode, 190, 200), "EYKT",
                                                               ifelse(nmfs.area %in% 650, "SEO","WESTWARD"))))))), 
         round.pounds = ifelse(cfec.whole.pounds < pounds, pounds, cfec.whole.pounds)) %>%  
  filter(!mgmt.area == "WESTWARD")

write.csv(ye_hal2,"Data_processing/Data/Harvests/CFEC Gross Earnings Data/data_check_ye_hal_summary.csv")


# run halibut numbers for Phil - halibut and rockfish, just replace/remove the ! before species.code to run for halibut or rockfish
ye_hal2 %>%
  filter(species.code == 200) %>% 
  group_by(year.landed, mgmt.area, permit.fishery, gear.description) %>% 
  summarise(round.lbs = sum(round.pounds)) -> halibut
table(halibut$permit.fishery, halibut$gear.description)   

write.csv(halibut, "Data_processing/Data/Harvests/halibut_catch_data_cfec.csv")



# next want to group by adfg number or permit holder name and DOL to combine trips, but need to get other species in for that. 
# end results is one line per trip that shows halibut and/or other target species with yelloweye 
# how to ensure these are directed halibut trips with bycatch listed 

#!permit.fishery %in% c("L99B",  "M05B",  "M06B",  "M17B",  "M61B",  "S03A",  "S05B",  "S15B",  "M07B",  "P09B",  "S01A", #might use this later
# "K09A",  "P91B",  "S04D",  "G03H",  "M99B",  "D09B",  "G01J",  "G03J",  "M09B", "G34H",  "M26B",
#"M91B",  "P07B", "S03H",  "G34K",  "S03T",  "S04K",   "M04B",  "M060",  "M16B",  "G34N",  "D91B",  "T09K",  
#  "S01K",  "S03M", "P91E",  "S03E",  "G01T",  "H34H",  "T09H",  "R18B",  "M02B",  
# "G34U",  "T09M",  "H07K",  "S04E",  "S04T",  "Z79B",  "T91H",  "S01E",  "P09E",  "Y05A", "Y06A",  "Y61A",  "Y26A",  
# "G34A",  "T91Q",  "O06B",  "S01L", "S04M",  "M37B",  "D91J",  "I25B",  "I26B",  "I05B", 
# "I06B",  "Y25A",  "M25B",  "P07E", "F06B", "MO7B",  "P17J",  "P17A",  "P17E",  "E",  
#"W22B", "P91A",  "Z99Z",  "M6AB",  "M26G",  "M06G",  "M7GB",  "M7GG",  "M7HB",  "M7HG",  "M7IB",  "M91G",
#"M7FB",  "M7FG",  "M7IG",  "M6AG", "M09G" , "M6BG",  "M6BB",  "M05G",  "I61B",  "S01M",  "M25G",  "T09E", "T91E"),
