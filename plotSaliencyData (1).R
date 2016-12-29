rm(list=ls())
library(R.matlab)
library(ggplot2)

# load in data from matlab
setwd("C:/Users/Anthony/Dropbox/Anthony/wamp/www/EyeTracking_PlosOne_2017") # you need to change this to whichever directory you're storing your file
d <- readMat('ratiosPsy_and_Image.mat')

# transpose image vaiables
d$ratioRMS <- t(d$ratioRMS)
d$ratioLocalRMS <- t(d$ratioLocalRMS)
d$ratioSalience <- t(d$ratioSalience)

# strip out valence
val <- d$ratioValence[,2]

data <- data.frame(ratio = (c(val,d$ratioArousal,d$ratioRMS,d$ratioLocalRMS,d$ratioSalience)), measure = (rep(c("rating - valence","rating - arousal","RMS - global","RMS - local","Salience - Koch"),each = 40)))
bp <- ggplot(data, aes(x=measure,y=ratio))+geom_boxplot()+ylab("social / non-social ratio")
#ggsave(bp, filename = "C:/Users/Anthony/Dropbox/Anthony/Bhismalab/Papers/FV_2010_Data/socnonsoc_ratio_measures.pdf", width=10, height=10)
ggsave(bp, filename = "C:/Users/Anthony/Dropbox/Anthony/wamp/www/EyeTracking_PlosOne_2017/Figure_1.png", width=10, height=10)
