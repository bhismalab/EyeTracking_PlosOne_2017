steiger.test<-function(r12,r13,r23){
  rmsqrd<-(r12$estimate^2+r13$estimate^2)/2;
  f=(1-r23$estimate)/(2*(1-rmsqrd));
  h=(1-(f*rmsqrd))/(1-rmsqrd);
  z12=0.5*log((1+r12$estimate)/(1-r12$estimate));
  z13=0.5*log((1+r13$estimate)/(1-r13$estimate));
  n=r12$parameter+2
  Z=(z12 - z13) * sqrt(n-3)/(sqrt(2*(1-r23$estimate)*h));
  p=2*pnorm(-abs(Z));  
  results<-list(Z_value=Z,P_Value=p);  
  return(results)
}

cohen.d<-function(group1,group2,mu=F){
  
  if(typeof(mu)=="logical"){
    result=abs(mean(group1)-mean(group2))/sqrt((sd(group1)^2+sd(group2)^2)/2)
  }
  else   {  
    result = abs(mean(group1)-mu)/sd(group1)  
  }
  return(result)
}

# Cousineau Function
cousineau.SE<-function(data_frame_inputted){
  across.means<-rowMeans(data_frame_inputted);
  supreme.mean<-mean(across.means);
  adjust.means<-across.means-supreme.mean;
  corrected.values<-data_frame_inputted-adjust.means;
  corrected.means<-colMeans(corrected.values, dims = 1);  
  corrected.means.array<-data.frame(corrected.means)
  corrected.means.array<-data.frame(rep(corrected.means.array,dim(corrected.values)[1]))
  
  corrected.means.array <- as.data.frame(t(corrected.means.array))
  colnames(corrected.means.array) <- corrected.means.array[1, ]
  corr.dev.means<-corrected.values-corrected.means.array
  corr.adj.sd<-corr.dev.means*corr.dev.means;
  corrected.SEs<-(sqrt(colSums(corr.adj.sd/(dim(corr.dev.means)[1]-1))))/sqrt(dim(corr.dev.means)[1])                         
  results<-corrected.SEs;
  return(results)
}


library(ggplot2)
library(grid)
# data wrangling made easier


image.list=read.csv("image_list.csv")

eq.correlations=data.frame(unsc.r.val=rep(0,length(image.list)),
                           unsc.p.val=rep(0,length(image.list)),
                           scr.r.val=rep(0,length(image.list)),
                           scr.p.val=rep(0,length(image.list)),
                           steiger.p=rep(0,length(image.list)))
### code to remove any image

big.pb <- winProgressBar(title = "big progress bar", min = 0,
                     max = length(image.list), width = 300)


for(all_of_it in 1:length(image.list)){
  #print(image.list[all_of_it])

  setWinProgressBar(big.pb, all_of_it, title=paste( round(all_of_it/length(image.list)*100, 0),"% done"))



# setwd("C:/Users/Anthony/Dropbox/Anthony/wamp/www/EyeTracking_PlosOne_2017")

setwd("R:/BhismaLab/anthony/Studies & Data sets/2010-2011/Freeviewing/Github/EyeTracking_PlosOne_2017")


#GE analysis
ge.raw<-read.csv("GE_data.csv")
fv.data.raw<-read.csv("FV_data_raw.csv")

############################################
### this is only for creating the table of how the data looks when removing each image pairing
### remove image.list[all_of_it]
fv.data.raw$image1[fv.data.raw$image1==image.list[all_of_it]]<-NA
fv.data.raw<-na.omit(fv.data.raw)
############################################


ge.raw$pptrial<-paste(ge.raw$PP,ge.raw$trial.number)

#reducing the data to first saccades
ge.raw.trials<-ge.raw[!duplicated(ge.raw$pptrial),]

## converting the data to matrices to facilitate analysis

ge.raw.trials$CURRENT_SAC_AMPLITUDE=data.matrix(ge.raw.trials$CURRENT_SAC_AMPLITUDE)
ge.raw.trials$CURRENT_SAC_ANGLE=data.matrix(ge.raw.trials$CURRENT_SAC_ANGLE)
ge.raw.trials$CURRENT_SAC_START_TIME=data.matrix(ge.raw.trials$CURRENT_SAC_START_TIME)

typeof(ge.raw.trials$CURRENT_SAC_ANGLE)
typeof(ge.raw.trials$CURRENT_SAC_AMPLITUDE)



#Identifying Scrambled
ge.raw.trials$Scrambled[ge.raw.trials$trial.code>40]<-"Scrambled"
ge.raw.trials$Scrambled[ge.raw.trials$trial.code<41]<-"Unscrambled"

#new trial code that ignores scrambled
ge.raw.trials$red.trial.code<-ge.raw.trials$trial.code 
ge.raw.trials$red.trial.code[ge.raw.trials$red.trial.code>40]<-ge.raw.trials$red.trial.code[ge.raw.trials$red.trial.code>40]-40

#clear "." trials
ge.raw.trials$CURRENT_SAC_ANGLE[ge.raw.trials$CURRENT_SAC_ANGLE=="."]<-NA
ge.raw.trials$CURRENT_SAC_ANGLE[ge.raw.trials$CURRENT_SAC_ANGLE==""]<-NA
ge.raw.trials<-na.omit(ge.raw.trials)

time.exclusion.trial=0

#Identifying sociality of first saccade
for (i in 1:length(ge.raw.trials$PP)) {
  if (ge.raw.trials$red.trial.code[i]<21 && ge.raw.trials$CURRENT_SAC_ANGLE[i]>0 || 
        ge.raw.trials$red.trial.code[i]>20 && ge.raw.trials$CURRENT_SAC_ANGLE[i]<0){
    ge.raw.trials$sociality[i]<-"Social"
  } else {
    ge.raw.trials$sociality[i]<-"Nonsocial"
  }     
  
  
  #Deviaton angle (as opposed to absolute angle)
  trial.angle=as.numeric(ge.raw.trials$CURRENT_SAC_ANGLE[i])
  valid.saccade=0
  # If images presented on side 1
  if(ge.raw.trials$side[i]==1){
    #if social image presented in the Top half
    if(ge.raw.trials$up[i]==1){
      if(trial.angle > 90){ ##this used to be greater than 90
        ge.raw.trials$dev.angle[i]=180-trial.angle
        valid.saccade=1
      }

      if (trial.angle < -90){
        ge.raw.trials$dev.angle[i]=-180-trial.angle
        valid.saccade=1
      }    
    }
    ## if social images presented in Bottom half
    if(ge.raw.trials$up[i]==2){
      if(trial.angle > 90){
        ge.raw.trials$dev.angle[i]=trial.angle-180
        valid.saccade=1
      } 
      if (trial.angle < -90){
        ge.raw.trials$dev.angle[i]=180+trial.angle
        valid.saccade=1
      }
    }
  }
  # If images presented on the right
  
  if(ge.raw.trials$side[i]==2){
    # and social image presented in top-half  
   if(ge.raw.trials$up[i]==1){
    if(trial.angle < 90 && trial.angle > 0){
      ge.raw.trials$dev.angle[i]=trial.angle
      valid.saccade=1
      }
     if (trial.angle > -90 && trial.angle < 0){
        ge.raw.trials$dev.angle[i]=trial.angle
        valid.saccade=1
    }    
   }
    # and social images presented in bottom-half
  if(ge.raw.trials$up[i]==2){
    if(trial.angle < 90 && trial.angle > 0){
      ge.raw.trials$dev.angle[i]=-trial.angle
      valid.saccade=1
      } 
    if (trial.angle > -90 && trial.angle < 0){
      ge.raw.trials$dev.angle[i]=abs(trial.angle)
      valid.saccade=1
      }
    }
  }
  if (valid.saccade==0){
    ge.raw.trials$dev.angle[i]=NA
  }
  
  #Identifying trials as NA in which participants first saccade was outside of 70-500ms of trial onset and highlighting relevant data
  if(ge.raw.trials$CURRENT_SAC_START_TIME[i]<70){
    ge.raw.trials$time.exclusion.70.lower[i]=1
  } else {
    ge.raw.trials$time.exclusion.70.lower[i]=0
  }

  if(ge.raw.trials$CURRENT_SAC_START_TIME[i]>500){
    ge.raw.trials$time.exclusion.500.higher[i]=1
  } else {
    ge.raw.trials$time.exclusion.500.higher[i]=0
  }

}

ge.raw.trials<-na.omit(ge.raw.trials)

#summarising data for each participant
summary.ge<-data.frame(pps=unique(ge.raw$PP))


# loop through participants to count time excluded trials
for(i in 1:length(summary.ge$pps)){
  this.data<-subset(ge.raw.trials,PP==summary.ge$pps[i])
  summary.ge$time.exclusion.70.lower[i] = sum(this.data$time.exclusion.70.lower)
  summary.ge$time.exclusion.500.higher[i] = sum(this.data$time.exclusion.500.higher)
}


ge.raw.trials$pptrial[ge.raw.trials$CURRENT_SAC_START_TIME<70]<-NA
ge.raw.trials$pptrial[ge.raw.trials$CURRENT_SAC_START_TIME>500]<-NA
ge.raw.trials<-na.omit(ge.raw.trials)



for(i in 1:length(summary.ge$pps)){
  this.data<-subset(ge.raw.trials,PP==summary.ge$pps[i])
  #mean angle for unscrambled
  unscrambled.angle.data<-subset(this.data,select=c(dev.angle,Scrambled),Scrambled=="Unscrambled")
  summary.ge$angle.unsc[i]=mean(unscrambled.angle.data$dev.angle)
  #mean angle for scrambled
  scrambled.angle.data<-subset(this.data,select=c(dev.angle,Scrambled),Scrambled=="Scrambled")
  summary.ge$angle.scr[i]=mean(scrambled.angle.data$dev.angle)
  
  #mean latency for unscrambled social
  unscrambled.lat.soc.data<-subset(this.data,select=c(CURRENT_SAC_START_TIME,Scrambled,sociality),
                             Scrambled=="Unscrambled" & sociality=="Social")
  summary.ge$lat.soc.unscr[i]=mean(unscrambled.lat.soc.data$CURRENT_SAC_START_TIME)
  

  #mean latency for unscrambled nonsocial
  unscrambled.lat.nonsoc.data<-subset(this.data,select=c(CURRENT_SAC_START_TIME,Scrambled,sociality),
                             Scrambled=="Unscrambled" & sociality=="Nonsocial")
  summary.ge$lat.nonsoc.unscr[i]=mean(unscrambled.lat.nonsoc.data$CURRENT_SAC_START_TIME)

  #mean latency for scrambled social
  scrambled.lat.soc.data<-subset(this.data,select=c(CURRENT_SAC_START_TIME,Scrambled,sociality),
                                 Scrambled=="Scrambled" & sociality=="Social")
  summary.ge$lat.soc.scr[i]=mean(scrambled.lat.soc.data$CURRENT_SAC_START_TIME)
  
  
  #mean latency for scrambled nonsocial
  scrambled.lat.nonsoc.data<-subset(this.data,select=c(CURRENT_SAC_START_TIME,Scrambled,sociality),
                                    Scrambled=="Scrambled" & sociality=="Nonsocial")
  summary.ge$lat.nonsoc.scr[i]=mean(scrambled.lat.nonsoc.data$CURRENT_SAC_START_TIME)
  
  ### Count no. valid unscrambled and scrambled trials
  
  #unscrambled
  summary.ge$unsc.trials[i]=sum(this.data$Scrambled=="Unscrambled")
  #scrambled
  summary.ge$scr.trials[i]=sum(this.data$Scrambled=="Scrambled")

} 

summary.ge$pps[summary.ge$unsc.trials<30]<-NA
summary.ge$pps[summary.ge$scr.trials<30]<-NA

summary.ge<-na.omit(summary.ge)

ge.total.time.exclusions=sum(summary.ge$time.exclusion.500.higher+summary.ge$time.exclusion.70.lower)
ge.total.time.exclusions

#total trials for each condition (scrambled/unscrambled)
length(summary.ge$pps)*80         #total number of trials
ge.total.time.exclusions/(length(summary.ge$pps)*80) #total number of exclusions
ge.total.time.exclusions/(length(summary.ge$pps)) #average number of trials per participant


#GE summary is now ready


################# FV analysis ###############



fv.data.raw$image.code[fv.data.raw$trial_1>40]<-2
fv.data.raw$image.code[fv.data.raw$trial_1<41]<-1


#progress bar
pb <- winProgressBar(title = "progress bar", min = 0,
                     max = length(fv.data.raw$RECORDING_SESSION_LABEL), width = 300)

fv.data.raw$left.dur=NA
fv.data.raw$right.dur=NA
fv.data.raw$social.dur=NA
fv.data.raw$nonsoc.dur=NA

for(i in 2:length(fv.data.raw$RECORDING_SESSION_LABEL)){
#for(i in 2:1000){  
    setWinProgressBar(pb, i, title=paste( round(i/length(fv.data.raw$RECORDING_SESSION_LABEL)*100, 0),"% done"))

  valid.saccade=0
  valid.saccade.ant=0
  
  #checking whether this saccade is from the same or different trial to previous one
  if (fv.data.raw$RECORDING_SESSION_LABEL[i]==fv.data.raw$RECORDING_SESSION_LABEL[i-1] &&
      fv.data.raw$trial_1[i]==fv.data.raw$trial_1[i-1]){
    #i.e. same trial
    # Note: saccade start time indicates the time that the participant made a saccade AWAY from an image.
    #confirm that participant was looking within valid y-axis
    if (fv.data.raw$CURRENT_SAC_START_Y[i]<822 && fv.data.raw$CURRENT_SAC_START_Y[i]>378) {
      #identify whether the image was within x co-ordinates of the left image
      if (fv.data.raw$CURRENT_SAC_START_X[i]<712 && fv.data.raw$CURRENT_SAC_START_X[i]> 168){
        #adding gaze duration to left image per saccade
        fv.data.raw$left.dur[i]=fv.data.raw$CURRENT_SAC_START_TIME[i]-fv.data.raw$CURRENT_SAC_START_TIME[i-1]
        fv.data.raw$right.dur[i]=NA
        
        #identifying whether left is social or nonsocial
        if(fv.data.raw$side[i]==1){
          fv.data.raw$social.dur[i]=fv.data.raw$CURRENT_SAC_START_TIME[i]-fv.data.raw$CURRENT_SAC_START_TIME[i-1]
          fv.data.raw$nonsoc.dur[i]=NA
        }
        else {
          fv.data.raw$nonsoc.dur[i]=fv.data.raw$CURRENT_SAC_START_TIME[i]-fv.data.raw$CURRENT_SAC_START_TIME[i-1]
          fv.data.raw$social.dur[i]=NA
        }
        
        valid.saccade=1
      } 
      #identify whether the image was within x co-ordinates of right image  
      if (fv.data.raw$CURRENT_SAC_START_X[i]>888 && fv.data.raw$CURRENT_SAC_START_X[i] < 1432){
        #adding gaze duarion to the right image per saccade

        fv.data.raw$right.dur[i]=fv.data.raw$CURRENT_SAC_START_TIME[i]-fv.data.raw$CURRENT_SAC_START_TIME[i-1]
        fv.data.raw$left.dur[i]=NA
        
        #identifying whether right is social or nonsocial
        
        if(fv.data.raw$side[i]==2){
          fv.data.raw$social.dur[i]=fv.data.raw$CURRENT_SAC_START_TIME[i]-fv.data.raw$CURRENT_SAC_START_TIME[i-1]
          fv.data.raw$nonsoc.dur[i]=NA
        }
        else {
          fv.data.raw$nonsoc.dur[i]=fv.data.raw$CURRENT_SAC_START_TIME[i]-fv.data.raw$CURRENT_SAC_START_TIME[i-1]
          fv.data.raw$social.dur[i]=NA
        }
        valid.saccade=1
      }
    }

    if (valid.saccade==0){
      fv.data.raw$left.dur[i]=NA
      fv.data.raw$right.dur[i]=NA
      fv.data.raw$social.dur[i]=NA
      fv.data.raw$nonsoc.dur[i]=NA
    }
      
    
  } else {
    #i.e. start of new trial
    #no duration information for first saccade
    fv.data.raw$left.dur[i]=NA
    fv.data.raw$right.dur[i]=NA
    fv.data.raw$social.dur[i]=NA
    fv.data.raw$nonsoc.dur[i]=NA
  }
  
}
close(pb)

summary.fv<-data.frame(pps=unique(fv.data.raw$RECORDING_SESSION_LABEL))

for (i in 1:length(summary.fv$pps)){

  this.data<-subset(fv.data.raw,RECORDING_SESSION_LABEL==summary.fv$pps[i])
  
  # Number of First saccades to each stimulus type and number of valid trials (unscrambled and scrambled)
    unsc.soc.first=0
    unsc.nonsoc.first=0
    scr.soc.first=0
    scr.nonsoc.first=0
    unsc.trials=0
    scr.trials=0
      for(j in 1:80){
      trial.data=subset(this.data,trial_1==j)
      if (!length(trial.data$RECORDING_SESSION_LABEL)==0) {
        soc.responses=trial.data$social.dur
        nonsoc.responses=trial.data$nonsoc.dur
        soc.responses=na.omit(soc.responses)
        nonsoc.responses=na.omit(nonsoc.responses)
        if(length(soc.responses)>0 || length(nonsoc.responses)>0){
          if(trial.data$image.code[1]==1){
            unsc.trials=unsc.trials+1
          }
          if(trial.data$image.code[1]==2){
            scr.trials=scr.trials+1
          }
        }
      
        for (k in 1:length(trial.data$RECORDING_SESSION_LABEL)){
          #social
          if(!is.na(trial.data$social.dur[k])){
            #unscrambled
            if(trial.data$image.code[k]==1){
              unsc.soc.first=unsc.soc.first+1
            } else {
              #scrambled
              scr.soc.first=scr.soc.first+1
            }
            break
          }
          #nonsocial
          if(!is.na(trial.data$nonsoc.dur[k])){
            #unscrambled
            if(trial.data$image.code[k]==1){
              unsc.nonsoc.first=unsc.nonsoc.first+1
            } else {
              #scrambled
              scr.nonsoc.first=scr.nonsoc.first+1
            }
            break
          }
        }  
      }
      
      }
    #No. first fixations to each stimulus type
    summary.fv$first.unsc.soc[i]=unsc.soc.first
    summary.fv$first.unsc.nonsoc[i]=unsc.nonsoc.first
    summary.fv$first.scr.soc[i]=scr.soc.first
    summary.fv$first.scr.nonsoc[i]=scr.nonsoc.first
    
    #No. valid trials for unscrambled and scrambled
    summary.fv$unsc.trials[i]=unsc.trials
    summary.fv$scr.trials[i]=scr.trials
  
    #gaze duration for unscrambled social
    
    summary.fv$unsc.soc.dur[i]=sum(na.omit(this.data$social.dur[this.data$image.code==1]))
    summary.fv$unsc.soc.dur[i]
    #gaze duration for unscrambled nonsocial
    summary.fv$unsc.nonsoc.dur[i]=sum(na.omit(this.data$nonsoc.dur[this.data$image.code==1]))
    summary.fv$unsc.nonsoc.dur[i]
    
    #gaze duration for scrambled social
    summary.fv$scr.soc.dur[i]=sum(na.omit(this.data$social.dur[this.data$image.code==2]))
    summary.fv$scr.soc.dur[i]
    #gaze duration for scrambled nonsocial
    summary.fv$scr.nonsoc.dur[i]=sum(na.omit(this.data$nonsoc.dur[this.data$image.code==2]))
    summary.fv$scr.nonsoc.dur[i]
    
} #FV summary is now ready

#remove participant(s) with fewer than 75% of trials in which gaze data was captured in either condition
summary.fv$pps[summary.fv$unsc.trials<30]=NA
summary.fv$pps[summary.fv$scr.trials<30]=NA

summary.fv=na.omit(summary.fv)

#proportion data
summary.fv$unsc.prop=summary.fv$unsc.soc.dur/(summary.fv$unsc.soc.dur+summary.fv$unsc.nonsoc.dur)
summary.fv$scr.prop=summary.fv$scr.soc.dur/(summary.fv$scr.soc.dur+summary.fv$scr.nonsoc.dur)


################ EQ data and lining up with pp numbers #################

#EQ data has been processed as outlined in Baron-Cohen et al (2001)
EQ.data <- read.csv("EQ_Data.csv")  
summary.EQ<-subset(EQ.data,select = c(pps,EQ.Total))
  
### removing participants with less than 30 valid trials in unscrambled or scrambled condition ###


#Gender info lined up
gender.data<-read.csv("Gender_Data.csv")
summary.fv=merge(summary.fv,gender.data,by = "pps")
summary.ge=merge(gender.data,summary.ge,by = "pps")

# lineup data between all three

summary.ge.eq=merge(summary.EQ,summary.ge,by = "pps")
summary.fv.eq=merge(summary.EQ,summary.fv,by = "pps")

# Summary of gender numbers for each analysis

sum(summary.fv$Gender=="f")#females in freeview
sum(summary.ge$Gender=="f")#females in global effect
sum(summary.fv.eq$Gender=="f")#females in freeview EQ correlation
sum(summary.ge.eq$Gender=="f")#females in global effect EQ correlation

#saving .csv files of summaries ofGE and FV data with and without EQ data

write.csv(summary.ge, file="GlobalEffectSummary.csv")
write.csv(summary.fv, file="FreeviewSummary.csv")

write.csv(summary.ge.eq, file="GlobalEffectSummary_EQ.csv")
write.csv(summary.fv.eq, file="FreeviewSummary_EQ.csv")


######ANALYSIS######

#GE task

##Main effects

#latency analysis

mean(summary.ge$lat.soc.unscr)
sd(summary.ge$lat.soc.unscr)
mean(summary.ge$lat.nonsoc.unscr)
sd(summary.ge$lat.nonsoc.unscr)
mean(summary.ge$lat.soc.scr)
sd(summary.ge$lat.soc.scr)
mean(summary.ge$lat.nonsoc.scr)
sd(summary.ge$lat.nonsoc.scr)

#one sampled t-tests
t.test(summary.ge$angle.unsc, mu=0)
cohen.d(summary.ge$angle.unsc, mu=0)
t.test(summary.ge$angle.scr, mu=0)
cohen.d(summary.ge$angle.scr, mu=0)

#paired
t.test(summary.ge$angle.unsc,summary.ge$angle.scr, paired=T)
cohen.d(summary.ge$angle.unsc,summary.ge$angle.scr)


#correlations

#EQ GE Unscrambled
cor.test(summary.ge.eq$EQ.Total,summary.ge.eq$angle.unsc)
#EQ GE Scrambled
cor.test(summary.ge.eq$EQ.Total,summary.ge.eq$angle.scr)



#FV task
#Main effects
#one sample tests
t.test(summary.fv$unsc.prop, mu=.5)
cohen.d(summary.fv$unsc.prop,mu = .5)
t.test(summary.fv$scr.prop, mu=.5)
cohen.d(summary.fv$scr.prop,mu = .5)

#paired
t.test(summary.fv$unsc.prop,summary.fv$scr.prop, paired=T)
cohen.d(summary.fv$unsc.prop,summary.fv$scr.prop)

#correlations
#EQ FV Unscrambled
#prop
eq.fv.unsc.cor=cor.test(summary.fv.eq$EQ.Total,summary.fv.eq$unsc.prop)
eq.fv.unsc.cor
#EQ FV Scrambled
eq.fv.scr.cor=cor.test(summary.fv.eq$EQ.Total,summary.fv.eq$scr.prop)
eq.fv.scr.cor

#Steiger test
FV.Unsc.Scr.cor=cor.test(summary.fv.eq$unsc.prop,summary.fv.eq$scr.prop)
steiger.result=steiger.test(eq.fv.unsc.cor,eq.fv.scr.cor,FV.Unsc.Scr.cor)

eq.correlations$unsc.r.val[all_of_it]=eq.fv.unsc.cor$estimate
eq.correlations$unsc.p.val[all_of_it]=eq.fv.unsc.cor$p.value
eq.correlations$scr.r.val[all_of_it]=eq.fv.scr.cor$estimate
eq.correlations$scr.p.val[all_of_it]=eq.fv.scr.cor$p.value
eq.correlations$steiger.p[all_of_it]=steiger.result$P_Value

}

close(big.pb) 

write.csv(eq.correlations,"eq_fv_correlations.csv")

#### FIGURES ####

#GE data
#bar chart

ge.cousineau.se=subset(summary.ge, select = (c(angle.unsc,angle.scr)))


ge.dat=data.frame(Condition=c("Unscrambled","Scrambled"),
                  Angle=c(mean(summary.ge$angle.unsc),mean(summary.ge$angle.scr)),
                  SE=cousineau.SE(ge.cousineau.se))

ge.dat=transform(ge.dat, Condition=reorder(Condition, order (Condition,decreasing = TRUE)))

soc.angle.plot <-ggplot(ge.dat,aes(x=Condition, y=Angle, fill=Condition)) +
  coord_cartesian(ylim=c(-.2, 1.8)) +
  #  zlab("test") +
  geom_bar(position=position_dodge(),stat="identity", fill=c("azure 4","grey"))+    #colour="black"
  geom_errorbar(aes(ymin=Angle-SE,ymax=Angle+SE),width=.2,position=position_dodge(.9), size=1)+
  #labs(title = "Angle towards social images")+
  xlab("")+
  ylab("Average deviation to social images (Angle)")+
  guides(fill=FALSE) +
  theme_bw() +
  theme(plot.background = element_blank(),        
        text = element_text(size=20,face = "bold")        
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()) +annotate("text",x=2.5,y=1.7,label="(a)", size=10)+
        geom_segment(aes(x=1,y=1.55,xend=1,yend=1.6), size=1) +
        geom_segment(aes(x=1,y=1.6,xend=2,yend=1.6), size=1) +
        geom_segment(aes(x=2,y=1.6,xend=2,yend=1.55),size=1)+annotate("text",x=1.5,y=1.7,label="p<.001", size=8)+
        geom_segment(aes(x=1.5,y=1.6,xend=1.5,yend=1.65),size=1)
soc.angle.plot

ggsave(soc.angle.plot, filename = "GE_Un_v_Scrambled_bar.png", width=10, height=10)



#FV graph

#bar chart

fv.cousineau.se=subset(summary.fv, select = (c(unsc.prop,scr.prop)))


fv.dat=data.frame(Condition=c("Unscrambled","Scrambled"),
                  prop=c(mean(summary.fv$unsc.prop),mean(summary.fv$scr.prop)),
                  SE=cousineau.SE(fv.cousineau.se))

fv.dat=transform(fv.dat, Condition=reorder(Condition, order (Condition,decreasing = TRUE)))

soc.dur.plot <-ggplot(fv.dat,aes(x=Condition, y=prop, fill=Condition)) +
  coord_cartesian(ylim=c(.4, .63)) +
  #  zlab("test") +
  geom_bar(position=position_dodge(),stat="identity", fill=c("azure 4","grey"))+    #colour="black"
  geom_errorbar(aes(ymin=prop-SE,ymax=prop+SE),width=.2,position=position_dodge(.9), size=1)+
  #labs(title = "Angle towards social images")+
  geom_hline(aes(yintercept=0.5), color="black", linetype="dashed", size=1)+
  xlab("")+
  ylab("Proportion of gaze duration to social images")+
  guides(fill=FALSE) +
  theme_bw() +
  theme(plot.background = element_blank(),        
        text = element_text(size=20,face = "bold")        
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()) +annotate("text",x=2.5,y=.62,label="(b)", size=10)+
        geom_segment(aes(x=1,y=.58,xend=1,yend=.585), size=1) +
        geom_segment(aes(x=1,y=.585,xend=2,yend=.585), size=1) +
        geom_segment(aes(x=2,y=.585,xend=2,yend=.58),size=1)+annotate("text",x=1.5,y=.597, label="p<.001", size=8)+
        geom_segment(aes(x=1.5,y=.585,xend=1.5,yend=.59),size=1)
soc.dur.plot

ggsave(soc.dur.plot, filename = "FV_Un_v_Scrambled_bar.png", width=10, height=10)


#EQ FV correlation

EQ <- rep(summary.fv.eq$EQ.Total,2)
EQ_FV <- c(summary.fv.eq$unsc.prop, summary.fv.eq$scr.prop)
Condition <- c(rep("Unscrambled",length(summary.fv.eq$pps)), rep("Scrambled",length(summary.fv.eq$pps)))
# convert to data frame and omit missing values
datEQ_FV <- data.frame(EQ,EQ_FV,Condition)
datEQ_FV <- na.omit(datEQ_FV)

# Order scrambled and unscrambled in order of "Unscrambled" then "scrambled"
# adding column for desired order
datEQ_FV$order=c(rep(1,length(summary.fv.eq$pps)),rep(2,length(summary.fv.eq$pps)))
datEQ_FV=transform(datEQ_FV, Condition=reorder(Condition, order(order, order, decreasing = FALSE)))


# a nicer color map
cbbPalette <- c("#000000", "#A0A0A0", "#CC79A7", "#F0E442", "#56B4E9", "#000000", "#E69F00", "#009E73", "#0072B2", "#D55E00")

scatEQ_FV <- ggplot(datEQ_FV, aes(x=EQ,y=EQ_FV,cond=Condition,color=Condition))+
  geom_point(aes(shape=Condition),size = 4) + 
  stat_smooth(method="lm", se=TRUE, fill="gray")+xlab("EQ") + ylab("Proportion Looking Time (Social)") + 
  scale_color_manual(values=cbbPalette)+
  theme(plot.background = element_blank()
        ,text = element_text(size=20, face="bold"))


scatEQ_FV
ggsave(scatEQ_FV, filename = "EQ_v_FV.png", width=10, height=8)






############################ POST REVIEWER ANALYSIS ##################################

## Range of angles & SD scatter plot (person by person)


EQ_GEPlot <- rep(summary.ge.eq$EQ.Total,2)
EQ_GE <- c(summary.ge.eq$angle.unsc, summary.ge.eq$angle.scr)
GE_Condition <- c(rep("Unscrambled",length(summary.ge.eq$pps)), rep("Scrambled",length(summary.ge.eq$pps)))
datEQ_GE <- data.frame(EQ_GEPlot,EQ_GE,GE_Condition)
datEQ_GE <- na.omit(datEQ_GE)

# Order scrambled and unscrambled in order of "Unscrambled" then "scrambled"
# adding column for desired order
datEQ_GE$order=c(rep(1,length(summary.ge.eq$pps)),rep(2,length(summary.ge.eq$pps)))
datEQ_GE=transform(datEQ_GE, GE_Condition=reorder(GE_Condition, order(order, order, decreasing = FALSE)))

scatEQ_GE <- ggplot(datEQ_GE, aes(x=EQ_GEPlot,y=EQ_GE,cond=GE_Condition,color=GE_Condition))+
  geom_point(aes(shape=GE_Condition),size = 4,cond=GE_Condition) + 
  stat_smooth(method="lm", se=TRUE, fill="gray")+xlab("EQ") + ylab("Angle of first saccade (positive is towards social)") + 
  scale_color_manual(values=cbbPalette)+
  theme(plot.background = element_blank()
        ,text = element_text(size=20, face="bold"))


scatEQ_GE

ggsave(scatEQ_GE, filename = "EQ_v_GE.png", width=10, height=8)


## proportion of first fixations to social vs. nonsocial in Freeview task

summary.fv$first.unsc.prop=summary.fv$first.unsc.soc/(summary.fv$first.unsc.soc+summary.fv$first.unsc.nonsoc)
summary.fv$first.scr.prop=summary.fv$first.scr.soc/(summary.fv$first.scr.soc+summary.fv$first.scr.nonsoc)

summary.fv.eq$first.unsc.prop=summary.fv.eq$first.unsc.soc/(summary.fv.eq$first.unsc.soc+summary.fv.eq$first.unsc.nonsoc)
summary.fv.eq$first.scr.prop=summary.fv.eq$first.scr.soc/(summary.fv.eq$first.scr.soc+summary.fv.eq$first.scr.nonsoc)




#main effects
t.test(summary.fv$first.unsc.prop, mu=.5)
cohen.d(summary.fv$first.unsc.prop,mu = .5)
t.test(summary.fv$first.scr.prop, mu=.5)
cohen.d(summary.fv$first.scr.prop,mu = .5)
t.test(summary.fv$first.unsc.prop, summary.fv$first.scr.prop, paired=T)
cohen.d(summary.fv$first.unsc.prop, summary.fv$first.scr.prop)


#correlations
ff.eq.unsc.cor=cor.test(summary.fv.eq$EQ.Total,summary.fv.eq$first.unsc.prop)
ff.eq.scr.cor=cor.test(summary.fv.eq$EQ.Total,summary.fv.eq$first.scr.prop)
ff.unsc.scr.cor=cor.test(summary.fv.eq$first.unsc.prop,summary.fv.eq$first.scr.prop)
steiger.test(ff.eq.unsc.cor,ff.eq.scr.cor,ff.unsc.scr.cor)

cohen.d(summary.fv$unsc.prop,summary.fv$scr.prop)


### First fix plot

EQ_FFPlot <- rep(summary.fv.eq$EQ.Total,2)
EQ_FF <- c(summary.fv.eq$first.unsc.prop, summary.fv.eq$first.scr.prop)
FF_Condition <- c(rep("Unscrambled",length(summary.fv.eq$pps)), rep("Scrambled",length(summary.fv.eq$pps)))
datEQ_FF <- data.frame(EQ_FFPlot,EQ_FF,FF_Condition)
datEQ_FF <- na.omit(datEQ_FF)

# Order scrambled and unscrambled in order of "Unscrambled" then "scrambled"
# adding column for desired order
datEQ_FF$order=c(rep(1,length(summary.fv.eq$pps)),rep(2,length(summary.fv.eq$pps)))
datEQ_FF=transform(datEQ_FF, FF_Condition=reorder(FF_Condition, order(order, order, decreasing = FALSE)))

scatEQ_FF <- ggplot(datEQ_FF, aes(x=EQ_FFPlot,y=EQ_FF,cond=FF_Condition,color=FF_Condition))+
  geom_point(aes(shape=FF_Condition),size = 4,cond=FF_Condition) + 
  stat_smooth(method="lm", se=TRUE, fill="gray")+xlab("EQ") + ylab("proportion of first saccades to social images") + 
  scale_color_manual(values=cbbPalette)+
  theme(plot.background = element_blank()
        ,text = element_text(size=20, face="bold"))


scatEQ_FF

ggsave(scatEQ_FF, filename = "EQ_v_FF.png", width=10, height=8)



### mean latencies for social vs. nonsocial stimuli in GE task

mean(summary.ge$lat.soc.unscr)
mean(summary.ge$lat.nonsoc.unscr)
t.test(summary.ge$lat.soc.unscr,summary.ge$lat.nonsoc.unscr, paired=T)
mean(summary.ge$lat.soc.scr)
mean(summary.ge$lat.nonsoc.scr)
t.test(summary.ge$lat.soc.scr,summary.ge$lat.nonsoc.scr, paired=T)



### gender split analysis - same direction for both, but underpowered separately

male.summary.fv.eq=subset(summary.fv.eq,summary.fv.eq$Gender=="m")
female.summary.fv.eq=subset(summary.fv.eq,summary.fv.eq$Gender=="f")
cor.test(male.summary.fv.eq$EQ.Total,male.summary.fv.eq$unsc.prop)
cor.test(female.summary.fv.eq$EQ.Total,female.summary.fv.eq$unsc.prop)


### deletable - I think ####


### cronbach alpha calcs...?

# GE

ge.alpha.pps = unique(ge.raw.trials$PP)
ge.alpha.pairs = unique(ge.raw.trials$trial.code)
ge.alpha.matrix = matrix(nrow = length(ge.alpha.pps), ncol=length(ge.alpha.pairs))

for(i in 1:length(ge.alpha.pps)){
  for(j in 1:length(ge.alpha.pairs)){
    ge.alpha.matrix = ge.raw.trials$dev.angle
    [ge.raw.trials]
  }
}





