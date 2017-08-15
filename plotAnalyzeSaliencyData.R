rm(list=ls())
library(R.matlab)
library(ggplot2)
library(stringr)
library(boot)

# load in data from matlab file 
#d <- readMat('/home/taylorcp/Dropbox/Exp/bhismalab/bhismalab_year1/Reward_Images_imStats/ratiosPsy_and_Image.mat')
setwd("C:/Users/Anthony/Dropbox/Anthony/wamp/www/EyeTracking_PlosOne_2017") # you need to change this to whichever directory you're storing your file

visual_ratios <- read.csv("Visual_Ratios.csv")

# strip out valence

data <- data.frame(ratio = (c(visual_ratios$Valence,
                              visual_ratios$Arousal,
                              visual_ratios$RMS,
                              visual_ratios$LocalRMS,
                              visual_ratios$Salience)), 
                   measure = (rep(c("Rating - Valence","Rating - Arousal","RMS - global","RMS - local","Salience - Koch"),each = 40)))
data$newMeasure = str_wrap(data$measure, width=10)

bp <- ggplot(data, aes(x=newMeasure,y=ratio))+
  geom_boxplot()+
  theme(text = element_text(size=30))+
  ylab("Social / Non-social ratio")+
  xlab("")
bp

ggsave(bp, filename = "Figure_1.png", width=10, height=10)



# one-sample t-tests with null mu = 1, a social/non-social ratio of 1 is balanced
# all tests are "two.sided" by default
t.test(visual_ratios$Arousal,mu = 1)
t.test(visual_ratios$Valence,mu = 1)
t.test(visual_ratios$RMS,mu = 1)
t.test(visual_ratios$LocalRMS,mu = 1)
t.test(visual_ratios$Salience, mu = 1) # this one comes out p < 0.04, because of outliers

# bootstrap ci using percentile method on the medians
fc_median <- function(d, i) {
  d2 <- d[i,]
  return(median(d2))
}

boot_median_rms <- boot(data = d$ratioRMS, statistic = fc_median, R = 999)
boot_median_local_rms <- boot(data = d$ratioLocalRMS, statistic = fc_median, R = 999)
boot_median_salience <- boot(data = d$ratioSalience, statistic = fc_median, R = 999)

boot.ci(boot_median_rms,type="perc")
boot.ci(boot_median_local_rms,type="perc")
boot.ci(boot_median_salience,type="perc")
