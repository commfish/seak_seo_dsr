#############################################################################
## Some script for plotting output from Discard deep dive.R scpue function output
##
library(ggplot2)
library(grid)

Disc<-read.csv("Output/SE.wcpue.ests_6.24.22.csv")

colnames(Disc)

Disc.io<-Disc[Disc$mngmt.divisions == "In.Out",]

Disc.subd<-Disc[Disc$mngmt.divisions =="SEdist",]
Disc.SEO<-Disc.subd[Disc.subd$mngmt.area == "CSEO" |
                      Disc.subd$mngmt.area == "NSEO" |
                      Disc.subd$mngmt.area == "SSEO" |
                      Disc.subd$mngmt.area == "EYKT",]
sub<-factor(c("EYKT","NSEO","CSEO","SSEO"))
Disc.SEO$mngmt.area<-factor(Disc.SEO$mngmt.area, levels=sub) 

Disc.SEI<-Disc.subd[Disc.subd$mngmt.area == "NSEI" |
                      Disc.subd$mngmt.area == "SSEI" ,]


cols<-c("raw WCPUE" = "red","weighted WCPUE"  ="blue")

ggplot(data=Disc.io, aes(x=Year)) +
  geom_line(aes(y= WCPUE32.mean, col= "red")) +
  geom_line(aes(y = wWCPUE32.mean, col = "blue")) +
  geom_ribbon(aes(ymin=WCPUE32.lo95ci, ymax=WCPUE32.hi95ci, fill="red"),
              alpha=0.2) +
  geom_ribbon(aes(ymin=wWCPUE32.lo95ci, ymax=wWCPUE32.hi95ci, fill="blue"),
              alpha=0.2) +
  facet_wrap(~mngmt.area) +
  ylab("WCPUE (kg yelloweye/kg halibut)") +
  scale_colour_manual(values = c("red", "blue"),
                      labels=c("raw", "weighted")) +
  scale_fill_manual(values = c("red", "blue"), breaks=c(),labels=c("raw", "weighted") ) +
  guides(col=guide_legend("WCPUE calculation:")) +
  theme(legend.position = "bottom")
  
ggsave(paste0("Figures/inside_outside.png"), dpi=300,  height=4.5, width=7.5, units="in")

ggplot(data=Disc.subd, aes(x=Year)) +
  geom_line(aes(y= WCPUE32.mean, col= "red")) +
  geom_line(aes(y = wWCPUE32.mean, col = "blue")) +
  geom_ribbon(aes(ymin=WCPUE32.lo95ci, ymax=WCPUE32.hi95ci, fill="red"),
              alpha=0.2) +
  geom_ribbon(aes(ymin=wWCPUE32.lo95ci, ymax=wWCPUE32.hi95ci, fill="blue"),
              alpha=0.2) +
  facet_wrap(~mngmt.area) +
  ylab("WCPUE (kg yelloweye/kg halibut)") +
  scale_colour_manual(values = c("red", "blue"),
                      labels=c("raw", "weighted")) +
  scale_fill_manual(values = c("red", "blue"), breaks=c(),labels=c("raw", "weighted") ) +
  guides(col=guide_legend("WCPUE calculation:")) +
  theme(legend.position = "bottom",
        #strip.text.x = element_text(
        #  color = c("red","blue")),
        strip.background=element_rect(
          fill=c("red","blue","green","red","blue","green")
        )
        )

ggsave(paste0("Figures/subdistricts.png"), dpi=300,  height=6, width=8, units="in")

ggplot(data=Disc.SEO, aes(x=Year)) +
  geom_line(aes(y= WCPUE32.mean, col= "red")) +
  geom_line(aes(y = wWCPUE32.mean, col = "blue")) +
  geom_ribbon(aes(ymin=WCPUE32.lo95ci, ymax=WCPUE32.hi95ci, fill="red"),
              alpha=0.2) +
  geom_ribbon(aes(ymin=wWCPUE32.lo95ci, ymax=wWCPUE32.hi95ci, fill="blue"),
              alpha=0.2) +
  facet_wrap(~mngmt.area) +
  ylab("WCPUE (kg yelloweye/kg halibut)") +
  scale_colour_manual(values = c("red", "blue"),
                      labels=c("raw", "weighted")) +
  scale_fill_manual(values = c("red", "blue"), breaks=c(),labels=c("raw", "weighted") ) +
  guides(col=guide_legend("WCPUE calculation:")) +
  theme(legend.position = "bottom")

ggsave(paste0("Figures/SEO.png"), dpi=300,  height=6, width=7, units="in")
