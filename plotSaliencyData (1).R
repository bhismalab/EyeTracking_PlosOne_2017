rm(list=ls())

#install.packages("R.matlab") # this is commented out for users who already have the package.
library(R.matlab)
library(ggplot2)

# load in data from matlab

#setwd("R:/BhismaLab/anthony/Studies & Data sets/2010-2011/Freeviewing/Github/EyeTracking_PlosOne_2017") # you need to change this to whichever directory you're storing your file
setwd("C:/Users/Anthony/Dropbox/Anthony/wamp/www/EyeTracking_PlosOne_2017") # you need to change this to whichever directory you're storing your file
d <- readMat('ratiosPsy_and_Image.mat')

###outlier checking

d_vars = c("ratioArousal","ratioValence","ratioRMS","ratioLocalRMS","ratioSalience")

outlier_array = data.frame(mean=rep(0,5))

d$ratioValence = d$ratioValence[,2]

for(i in 1:5){
  outlier_array$mean[i] =  mean(d[[d_vars[i]]])
  outlier_array$sd[i] = sd(d[[d_vars[i]]])
  this_t.test = t.test(d[[d_vars[i]]],mu=1)
  outlier_array$t_value[i] = this_t.test$statistic
  outlier_array$p_value[i] = this_t.test$p.value
  outlier_array$metric[i] = d_vars[i]
  
  ### outlier detection
  
  clean_data = d[[d_vars[i]]][abs(d[[d_vars[i]]]-outlier_array$mean[i])<3*outlier_array$sd[i]]
  
  
  outlier_array$clean_mean[i] =  mean(clean_data)
  outlier_array$clean_sd[i] = sd(clean_data)
  this_t.test = t.test(clean_data,mu=1)
  outlier_array$clean_t_value[i] = this_t.test$statistic
  outlier_array$clean_p_value[i] = this_t.test$p.value
  outlier_array$clean_metric[i] = d_vars[i]
  outlier_array$pics[i]=length(d[[d_vars[i]]])
  outlier_array$cleanpics[i]=length(clean_data)
  
}

### end of outlier checking


# transpose image vaiables
d$ratioRMS <- t(d$ratioRMS)
d$ratioLocalRMS <- t(d$ratioLocalRMS)
d$ratioSalience <- t(d$ratioSalience)

# strip out valence
val <- d$ratioValence[,2]

data <- data.frame(ratio = (c(val,d$ratioArousal,d$ratioRMS,d$ratioLocalRMS,d$ratioSalience)), measure = (rep(c("Rating - Valence","Rating - Arousal","RMS - global","RMS - local","Salience - Koch"),each = 40)))


bp <- ggplot(data, aes(x=measure,y=ratio))+
  geom_boxplot()+
  theme(text = element_text(size=20))+
  ylab("Social / Non-social ratio")+
  xlab("")
bp

ggsave(bp, filename = "C:/Users/Anthony/Dropbox/Anthony/wamp/www/EyeTracking_PlosOne_2017/Figure_1.png", width=10, height=10)


