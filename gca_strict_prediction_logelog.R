# Gambi, Gorrie, Pickering, Rabagliati (2017). The development of linguistic prediction: Predictions of sound and meaning in 2-to-5 year olds.
# Analysis of fixation proportions from Determiner Onset to Noun Onset + 100 ms (Prediction Window) using Growth Curves
# Experiment 1
#################### Statistical Analyses reported in the paper ##################
###################### Plan of analyses ##########################################
### We are going to do the following
### Elog Growth Curve Analysis both by subjects and items
### We do this separately for Phon and Sem, but do between-condition comparisons
### We compare each age-group to the age group immediately younger, with adults as the reference level
### We include Vocabulary as a continuous predictor
### We then re-run all the analyses split by occurrence (App == 1 vs. 2) -- see separate .R file


# ELOG GROWTH CURVE BY SUBJECTS
# Note that we can use TARGETFIXBIN both in combined Sem Phon and in separate analyses
gcaT<-read.table("SemPhonDOplus1300_trial.txt", header=T)
library(dplyr)
#(Following Brock Ferguson example)
#http://brockferguson.com/tutorials/R-gca-workshop/tutorial4-growth-curve-analysis.html#orthogonal-polynomial-growth-curve-analysis
bySubj <- gcaT %>%
  group_by(Participant,Condition,DetP,DetS,DET,AgeGroup,Age,BPVS,App,time) %>% # aggregate within time slices
  summarise(TARGETFIXBIN.m = mean(TARGETFIXBIN),yT = sum(TARGETFIXBIN), N_T = length(TARGETFIXBIN)) %>%
  mutate(elogT = log( (yT + .5) / (N_T - yT + .5) )) %>%  # empirical logit
  #wts = 1/(y + .5) + 1/(N - y + .5), # optional weights
  #Arcsin = asin(sqrt(PropAnimal))) # arcsin-sqrt
  ungroup()   

## Sem vs. Phon
t<-poly(unique(bySubj$time),1)
time<-as.vector(unique(bySubj$time))
t<-cbind(t,time)
t<-as.data.frame(t)
bySubj<-bySubj[order(bySubj$time),]
bySubj$t1<-NA
##gcaT$t2<-NA
for (i in (1:nrow(bySubj))){
  bySubj$t1[i]<-t[t$time==bySubj$time[i],1] 
  ##gcaT$t2[i]<-t[t$time==gcaT$time[i],2] 
}
summary(bySubj)
## since these analyses compare Phon to Sem, use DET
bySubj$DC<-ifelse(bySubj$DET=="Incongruent",-.5,.5)
bySubj$DCC<-scale(bySubj$DC, T, F)
## cannot add add vocabulary as a continuous, centred covariate because this analysis includes adults
#gcaT$VCC<-scale(gcaT$BPVS, T, F)
## Include Age using backward difference coding (backward because we compare each level to previous one)
# 3yo vs. 2yo
bySubj$A32<-ifelse(bySubj$AgeGroup=="2yo",-3/4,1/4)
# 45yo vs. 3yo
bySubj$A453<-ifelse(bySubj$AgeGroup=="2yo",-1/2,ifelse(bySubj$AgeGroup=="3yo",-1/2,1/2))
# Adult vs. 45yo
bySubj$AA45<-ifelse(bySubj$AgeGroup=="Adult",3/4,-1/4)
# Condition
bySubj$CC<-ifelse(bySubj$Condition=="Sem",-.5,.5)
bySubj$CCC<-scale(bySubj$CC, T, F)

library(lme4)
# Full model with Age contrasts
mT<-lmer(elogT~1+t1*DCC*CCC*(A32+A453+AA45)+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj)
summary(mT)

### only 2yo
bySubj2<-bySubj[bySubj$AgeGroup=="2yo",]
## since these analyses compare Phon to Sem, use DET
bySubj2$DC<-ifelse(bySubj2$DET=="Incongruent",-.5,.5)
bySubj2$DCC<-scale(bySubj2$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubj2$VCC<-scale(bySubj2$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubj2$ACC<-scale(bySubj2$Age, T, F)
# Condition
bySubj2$CC<-ifelse(bySubj2$Condition=="Sem",-.5,.5)
bySubj2$CCC<-scale(bySubj2$CC, T, F)

mT2<-lmer(elogT~1+t1*DCC*CCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj2)
summary(mT2)
confint(mT2,method="Wald")
#loglik
mT2.i<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+t1:DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj2)
anova(mT2,mT2.i)#1.5657      1     0.2108
mT2.l<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj2)
anova(mT2,mT2.l)#2.5199      1     0.1124

### only 3yo
bySubj3<-bySubj[bySubj$AgeGroup=="3yo",]
## since these analyses compare Phon to Sem, use DET
bySubj3$DC<-ifelse(bySubj3$DET=="Incongruent",-.5,.5)
bySubj3$DCC<-scale(bySubj3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubj3$VCC<-scale(bySubj3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubj3$ACC<-scale(bySubj3$Age, T, F)
# Condition
bySubj3$CC<-ifelse(bySubj3$Condition=="Sem",-.5,.5)
bySubj3$CCC<-scale(bySubj3$CC, T, F)

mT3<-lmer(elogT~1+t1*DCC*CCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj3)
summary(mT3)
confint(mT3,method="Wald")

#loglik
mT3.i<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+t1:DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj3)
anova(mT3,mT3.i)#9.6584      1   0.001885 **
mT3.l<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj3)
anova(mT3,mT3.l)#14.189      1  0.0001653 ***

### only 45yo
bySubj45<-bySubj[bySubj$AgeGroup=="4-5yo",]
## since these analyses compare Phon to Sem, use DET
bySubj45$DC<-ifelse(bySubj45$DET=="Incongruent",-.5,.5)
bySubj45$DCC<-scale(bySubj45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubj45$VCC<-scale(bySubj45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubj45$ACC<-scale(bySubj45$Age, T, F)
# Condition
bySubj45$CC<-ifelse(bySubj45$Condition=="Sem",-.5,.5)
bySubj45$CCC<-scale(bySubj45$CC, T, F)

mT45<-lmer(elogT~1+t1*DCC*CCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj45)
summary(mT45)
confint(mT45, method="Wald")

#loglik
mT45.i<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+t1:DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj45)
anova(mT45,mT45.i)#10.572      1   0.001148 **
mT45.l<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubj45)
anova(mT45,mT45.l)#18.161      1   2.03e-05 ***

### only Adults
bySubjA<-bySubj[bySubj$AgeGroup=="Adult",]
## since these analyses compare Phon to Sem, use DET
bySubjA$DC<-ifelse(bySubjA$DET=="Incongruent",-.5,.5)
bySubjA$DCC<-scale(bySubjA$DC, T, F)
## cannot add vocabulary
#bySubj3$VCC<-scale(bySubj3$BPVS, T, F)
# Condition
bySubjA$CC<-ifelse(bySubjA$Condition=="Sem",-.5,.5)
bySubjA$CCC<-scale(bySubjA$CC, T, F)

mTA<-lmer(elogT~1+t1*DCC*CCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubjA)
summary(mTA)
confint(mTA,method="Wald")
#loglik
mTA.i<-lmer(elogT~1+t1*DCC+t1*CCC+t1:DCC:CCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubjA)
anova(mTA,mTA.i)#21.477      1   3.58e-06 ***
mTA.l<-lmer(elogT~1+t1*DCC+t1*CCC+DCC:CCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:CCC)+(1+t1||Participant:DCC:CCC), data=bySubjA)
anova(mTA,mTA.l)#18.655      1  1.566e-05 ***

## Sem
bySubjS<-bySubj[bySubj$Condition=="Sem",]
## since these analyses look at Sem only use DetS
bySubjS$DC<-ifelse(bySubjS$DetS=="One",-.5,.5)
bySubjS$DCC<-scale(bySubjS$DC, T, F)
## cannot add add vocabulary as a continuous, centred covariate because this analysis includes adults
#gcaT$VCC<-scale(gcaT$BPVS, T, F)
## Include Age using backward difference coding (backward because we compare each level to previous one)
# 3yo vs. 2yo
bySubjS$A32<-ifelse(bySubjS$AgeGroup=="2yo",-3/4,1/4)
# 45yo vs. 3yo
bySubjS$A453<-ifelse(bySubjS$AgeGroup=="2yo",-1/2,ifelse(bySubjS$AgeGroup=="3yo",-1/2,1/2))
# Adult vs. 45yo
bySubjS$AA45<-ifelse(bySubjS$AgeGroup=="Adult",3/4,-1/4)

# Full model with Age contrasts
mS<-lmer(elogT~1+t1*DCC*(A32+A453+AA45)+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS)
summary(mS)
confint(mS,method="Wald")

#loglik
mS.i.32<-lmer(elogT~1+t1*DCC*(A453+AA45)+A32+t1:A32+t1:DCC:A32+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS)
anova(mS.i.32,mS)#10.19      1   0.001412 **
mS.i.A45<-lmer(elogT~1+t1*DCC*(A453+A32)+AA45+t1:AA45+t1:DCC:AA45+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS)
anova(mS.i.A45,mS)# 31.404      1  2.096e-08 ***
mS.l.32<-lmer(elogT~1+t1*DCC*(A453+AA45)+A32+t1:A32+DCC:A32+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS)
anova(mS.l.32,mS)#1.8695      1     0.1715
mS.l.A45<-lmer(elogT~1+t1*DCC*(A453+A32)+AA45+t1:AA45+DCC:AA45+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS)
anova(mS.l.A45,mS)#22.182      1   2.48e-06 ***

### now analyse each age group separately and add vocabulary in children analyses

## only 2 yo
bySubjS2<-bySubjS[bySubjS$AgeGroup=="2yo",]
## since these analyses look at Sem only use DetS
bySubjS2$DC<-ifelse(bySubjS2$DetS=="One",-.5,.5)
bySubjS2$DCC<-scale(bySubjS2$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjS2$VCC<-scale(bySubjS2$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjS2$ACC<-scale(bySubjS2$Age, T, F)

mS2<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS2)
summary(mS2)
confint(mS2,method="Wald")

#loglik
mS2.l<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS2)
anova(mS2.l,mS2)# 10.516      1   0.001183 **
mS2.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+t1:DCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS2)
anova(mS2,mS2.i)# 2.0213      1     0.1551

## only 3 yo
bySubjS3<-bySubjS[bySubjS$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
bySubjS3$DC<-ifelse(bySubjS3$DetS=="One",-.5,.5)
bySubjS3$DCC<-scale(bySubjS3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjS3$VCC<-scale(bySubjS3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjS3$ACC<-scale(bySubjS3$Age, T, F)

mS3<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS3)
summary(mS3)
confint(mS3,method="Wald")
#loglik
mS3.l<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS3)
anova(mS3.l,mS3)# 28.567      1  9.053e-08 ***
mS3.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+t1:DCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS3)
anova(mS3,mS3.i)# 25.318      1  4.863e-07 ***
mS3.l.voc<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS3)
anova(mS3.l.voc,mS3)#6.3736      1    0.01158 *

## only 45 yo
bySubjS45<-bySubjS[bySubjS$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
bySubjS45$DC<-ifelse(bySubjS45$DetS=="One",-.5,.5)
bySubjS45$DCC<-scale(bySubjS45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjS45$ACC<-scale(bySubjS45$Age, T, F)

mS45<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS45)
summary(mS45)
confint(mS45,method="Wald")
#loglik
mS45.l<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS45)
anova(mS45.l,mS45)#34.941      1  3.399e-09 ***
mS45.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+t1:DCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS45)
anova(mS45,mS45.i)#32.822      1   1.01e-08 ***
mS45.l.voc<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjS45)
anova(mS45.l.voc,mS45)#5.4914      1    0.01911 *

## only Adults
bySubjSA<-bySubjS[bySubjS$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
bySubjSA$DC<-ifelse(bySubjSA$DetS=="One",-.5,.5)
bySubjSA$DCC<-scale(bySubjSA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)

mSA<-lmer(elogT~1+t1*DCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjSA)
summary(mSA)
confint(mSA, method="Wald")
# loglik
mSA.l<-lmer(elogT~1+t1+DCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjSA)
anova(mSA,mSA.l)#84.457      1  < 2.2e-16 ***
mSA.i<-lmer(elogT~1+t1+DCC:t1+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjSA)
anova(mSA,mSA.i)#68.947      1  < 2.2e-16 ***

## Phon
bySubjP<-bySubj[bySubj$Condition=="Phon",]
## since these analyses look at Phon only use DetP
bySubjP$DC<-ifelse(bySubjP$DetP=="A",-.5,.5)
bySubjP$DCC<-scale(bySubjP$DC, T, F)
## cannot add add vocabulary as a continuous, centred covariate because this analysis includes adults
#gcaT$VCC<-scale(gcaT$BPVS, T, F)
## Include Age using backward difference coding (backward because we compare each level to previous one)
# 3yo vs. 2yo
bySubjP$A32<-ifelse(bySubjP$AgeGroup=="2yo",-3/4,1/4)
# 45yo vs. 3yo
bySubjP$A453<-ifelse(bySubjP$AgeGroup=="2yo",-1/2,ifelse(bySubjP$AgeGroup=="3yo",-1/2,1/2))
# Adult vs. 45yo
bySubjP$AA45<-ifelse(bySubjP$AgeGroup=="Adult",3/4,-1/4)

mP<-lmer(elogT~1+t1*DCC*(A32+A453+AA45)+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP)
summary(mP)
confint(mP, method="Wald")
mP.iA45<-lmer(elogT~1+t1*DCC*(A32+A453)+t1:AA45+AA45+t1:DCC:AA45+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP)
anova(mP,mP.iA45)#11.593      1  0.0006621 ***
mP.lA45<-lmer(elogT~1+t1*DCC*(A32+A453)+t1:AA45+AA45+DCC:AA45+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP)
anova(mP,mP.lA45)#28.742      1  8.268e-08 ***

### now analyse each age group separately and add vocabulary in children analyses

## only 2 yo
bySubjP2<-bySubjP[bySubjP$AgeGroup=="2yo",]
## since these analyses look at Sem only use DetS
bySubjP2$DC<-ifelse(bySubjP2$DetP=="A",-.5,.5)
bySubjP2$DCC<-scale(bySubjP2$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjP2$VCC<-scale(bySubjP2$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjP2$ACC<-scale(bySubjP2$Age, T, F)

mP2<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP2)
summary(mP2)
confint(mP2,method="Wald")
# loglik
mP2.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC+DCC:VCC+DCC:ACC+DCC:VCC:ACC+ t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP2)
anova(mP2,mP2.i)#0.0577      1     0.8101
mP2.l<-lmer(elogT~1+t1*VCC*ACC+DCC+DCC:VCC+DCC:ACC+DCC:VCC:ACC+ t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP2)
anova(mP2,mP2.l)#1.0594      1     0.3033
mP2.i.Voc<-lmer(elogT~1+t1*VCC*ACC+t1:DCC+DCC+DCC:ACC+DCC:VCC:ACC+ t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP2)
anova(mP2,mP2.i.Voc)#12.842      1  0.0003389 ***

## only 3 yo
bySubjP3<-bySubjP[bySubjP$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
bySubjP3$DC<-ifelse(bySubjP3$DetP=="A",-.5,.5)
bySubjP3$DCC<-scale(bySubjP3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjP3$VCC<-scale(bySubjP3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjP3$ACC<-scale(bySubjP3$Age, T, F)

mP3<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP3)
summary(mP3)
confint(mP3,method="Wald")
# loglik
mP3.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC+DCC:VCC+DCC:ACC+DCC:VCC:ACC+ t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP3)
anova(mP3,mP3.i)#5.0775      1    0.02424 *
mP3.l<-lmer(elogT~1+t1*VCC*ACC+DCC+DCC:VCC+DCC:ACC+DCC:VCC:ACC+ t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP3)
anova(mP3,mP3.l)#0.6649      1     0.4148
mP3.l.voc<-lmer(elogT~1+t1*VCC*ACC+DCC+DCC:VCC+DCC:ACC+DCC:VCC:ACC+ t1:DCC:ACC+t1:DCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP3)
anova(mP3,mP3.l.voc)#3.0437      1    0.08105 .

## only 45 yo
bySubjP45<-bySubjP[bySubjP$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
bySubjP45$DC<-ifelse(bySubjP45$DetP=="A",-.5,.5)
bySubjP45$DCC<-scale(bySubjP45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjP45$VCC<-scale(bySubjP45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjP45$ACC<-scale(bySubjP45$Age, T, F)

# note that this model does not include DCC:ACC:VCC in the fixed effect structure - the DCC effect is not significant in a model that contains the three-way interaction
mP45<-lmer(elogT~1+t1*DCC*VCC+t1*DCC*ACC+ACC:VCC+t1:ACC:VCC+t1:DCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP45)
summary(mP45)
confint(mP45, method="Wald")
#loglik
mP45.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC+DCC:VCC+DCC:ACC+ t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP45)
anova(mP45,mP45.i)#5.4954      1    0.01907 *
mP45.l<-lmer(elogT~1+t1*VCC*ACC+DCC+DCC:VCC+DCC:ACC+ t1:DCC:ACC+t1:DCC:VCC+t1:DCC:VCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjP45)
anova(mP45,mP45.l)# 0.0447      1     0.8326

## only Adults
bySubjPA<-bySubjP[bySubjP$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
bySubjPA$DC<-ifelse(bySubjPA$DetP=="A",-.5,.5)
bySubjPA$DCC<-scale(bySubjPA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)

mPA<-lmer(elogT~1+t1*DCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjPA)
summary(mPA)
confint(mPA,method="Wald")
#loglik
mPA.i<-lmer(elogT~1+t1+t1:DCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjPA)
anova(mPA.i,mPA)#23.892      1  1.019e-06 ***
mPA.l<-lmer(elogT~1+t1+DCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=bySubjPA)
anova(mPA.l,mPA)#46.522      1  9.061e-12 ***

###### ELOG BY ITEMS
gcaT<-read.table("SemPhonDOplus1300_trial.txt", header=T)
library(dplyr)

byItem <- gcaT %>%
  group_by(ItemN,Condition,DetP,DetS,DET,AgeGroup,Age,BPVS,App,time) %>% # aggregate within time slices
  summarise(TARGETFIXBIN.m = mean(TARGETFIXBIN),yT = sum(TARGETFIXBIN), N_T = length(TARGETFIXBIN)) %>%
  mutate(elogT = log( (yT + .5) / (N_T - yT + .5) )) %>%  # empirical logit
  #wts = 1/(y + .5) + 1/(N - y + .5), # optional weights
  #Arcsin = asin(sqrt(PropAnimal))) # arcsin-sqrt
  ungroup()   

## Sem vs. Phon
t<-poly(unique(byItem$time),1)
time<-as.vector(unique(byItem$time))
t<-cbind(t,time)
t<-as.data.frame(t)
byItem<-byItem[order(byItem$time),]
byItem$t1<-NA
##gcaT$t2<-NA
for (i in (1:nrow(byItem))){
  byItem$t1[i]<-t[t$time==byItem$time[i],1] 
  ##gcaT$t2[i]<-t[t$time==gcaT$time[i],2] 
}
summary(byItem)
## since these analyses compare Phon to Sem, use DET
byItem$DC<-ifelse(byItem$DET=="Incongruent",-.5,.5)
byItem$DCC<-scale(byItem$DC, T, F)
## cannot add add vocabulary as a continuous, centred covariate because this analysis includes adults
#gcaT$VCC<-scale(gcaT$BPVS, T, F)
## Include Age using backward difference coding (backward because we compare each level to previous one)
# 3yo vs. 2yo
byItem$A32<-ifelse(byItem$AgeGroup=="2yo",-3/4,1/4)
# 45yo vs. 3yo
byItem$A453<-ifelse(byItem$AgeGroup=="2yo",-1/2,ifelse(byItem$AgeGroup=="3yo",-1/2,1/2))
# Adult vs. 45yo
byItem$AA45<-ifelse(byItem$AgeGroup=="Adult",3/4,-1/4)
# Condition
byItem$CC<-ifelse(byItem$Condition=="Sem",-.5,.5)
byItem$CCC<-scale(byItem$CC, T, F)

# Full model with Age contrasts
#mT.I<-lmer(elogT~1+t1*DCC*CCC*(A32+A453+AA45)+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem)# did not converge
mT.I<-lmer(elogT~1+t1*DCC*CCC+(A32+A453+AA45)+(A32+A453+AA45):DCC+(A32+A453+AA45):DCC:CCC+(A32+A453+AA45):t1:DCC+(A32+A453+AA45):t1:DCC:CCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem)
summary(mT.I)
### only 2yo
byItem2<-byItem[byItem$AgeGroup=="2yo",]
## since these analyses compare Phon to Sem, use DET
byItem2$DC<-ifelse(byItem2$DET=="Incongruent",-.5,.5)
byItem2$DCC<-scale(byItem2$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItem2$VCC<-scale(byItem2$BPVS, T, F)
## add vocabulary as a continuous, centred covariate 
byItem2$ACC<-scale(byItem2$Age, T, F)
# Condition
byItem2$CC<-ifelse(byItem2$Condition=="Sem",-.5,.5)
byItem2$CCC<-scale(byItem2$CC, T, F)

#mT2.I<-lmer(elogT~1+t1*DCC*CCC*VCC+(1+t1*VCC||ItemN)+(1+t1*VCC||ItemN:DCC)+(1+t1*VCC||ItemN:CCC)+(1+t1*VCC||ItemN:DCC:CCC), data=byItem2)# did not converge
#mT2.I<-lmer(elogT~1+t1*DCC*CCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:CCC)+(1+t1+VCC+ACC||ItemN:DCC:CCC), data=byItem2)# did not converge
# remove random slopes by Voc and Age from most complex term
mT2.I<-lmer(elogT~1+t1*DCC*CCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem2)
summary(mT2.I)
confint(mT2.I, method="Wald")
#loglik
mT2.I.i<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+t1:DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem2)#did not converge
anova(mT2.I,mT2.I.i)#0.9144      1     0.3389
mT2.I.l<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem2)#did not converge
anova(mT2.I,mT2.I.l)#4.1236      1    0.04229 *

### only 3yo
byItem3<-byItem[byItem$AgeGroup=="3yo",]
## since these analyses compare Phon to Sem, use DET
byItem3$DC<-ifelse(byItem3$DET=="Incongruent",-.5,.5)
byItem3$DCC<-scale(byItem3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItem3$VCC<-scale(byItem3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItem3$ACC<-scale(byItem3$Age, T, F)
# Condition
byItem3$CC<-ifelse(byItem3$Condition=="Sem",-.5,.5)
byItem3$CCC<-scale(byItem3$CC, T, F)

mT3.I<-lmer(elogT~1+t1*DCC*CCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:CCC)+(1+t1+VCC+ACC||ItemN:DCC:CCC), data=byItem3)
summary(mT3.I)
confint(mT3.I, method="Wald")
#loglik
mT3.I.i<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+t1:DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:CCC)+(1+t1+VCC+ACC||ItemN:DCC:CCC), data=byItem3)
anova(mT3.I.i,mT3.I)#12.198      1  0.0004785 ***
mT3.I.l<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:CCC)+(1+t1+VCC+ACC||ItemN:DCC:CCC), data=byItem3)
anova(mT3.I.l,mT3.I)#22.686      1  1.907e-06 ***

### only 45yo
byItem45<-byItem[byItem$AgeGroup=="4-5yo",]
## since these analyses compare Phon to Sem, use DET
byItem45$DC<-ifelse(byItem45$DET=="Incongruent",-.5,.5)
byItem45$DCC<-scale(byItem45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItem45$VCC<-scale(byItem45$BPVS, T, F)
## add vocabulary as a continuous, centred covariate 
byItem45$ACC<-scale(byItem45$Age, T, F)
# Condition
byItem45$CC<-ifelse(byItem45$Condition=="Sem",-.5,.5)
byItem45$CCC<-scale(byItem45$CC, T, F)

#mT45.I<-lmer(elogT~1+t1*DCC*CCC*VCC+(1+t1+VCC||ItemN)+(1+t1+VCC||ItemN:DCC)+(1+t1+VCC||ItemN:CCC)+(1+t1+VCC||ItemN:DCC:CCC), data=byItem45)#did not converge
#mT45.I<-lmer(elogT~1+t1*DCC*CCC+VCC+DCC:VCC+DCC:CCC:VCC+t1:DCC:VCC+t1:DCC:CCC:VCC+(1+t1+VCC||ItemN)+(1+t1+VCC||ItemN:DCC)+(1+t1+VCC||ItemN:CCC)+(1+t1+VCC||ItemN:DCC:CCC), data=byItem45)# did not converge
#mT45.I<-lmer(elogT~1+t1*DCC*CCC+VCC+DCC:CCC:VCC+t1:DCC:CCC:VCC+(1+t1+VCC||ItemN)+(1+t1+VCC||ItemN:DCC)+(1+t1+VCC||ItemN:CCC)+(1+t1+VCC||ItemN:DCC:CCC), data=byItem45)#did not converge
#mT45.I<-lmer(elogT~1+t1*DCC*CCC+VCC+DCC:CCC:VCC+t1:DCC:CCC:VCC+(1+t1||ItemN)+(1+t1+VCC||ItemN:DCC)+(1+t1+VCC||ItemN:CCC)+(1+t1+VCC||ItemN:DCC:CCC), data=byItem45)# did not converge
# removed Voc and Age from random structure to aid convergence
mT45.I<-lmer(elogT~1+t1*DCC*CCC*VCC*ACC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem45)
summary(mT45.I)
confint(mT45.I,method="Wald")
#loglik
mT45.I.i<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+t1:DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem45)
anova(mT45.I,mT45.I.i)#10.342      1     0.0013 **
mT45.I.l<-lmer(elogT~1+t1*DCC*VCC*ACC+t1*CCC*VCC*ACC+DCC:CCC+t1:DCC:CCC:ACC+t1:DCC:CCC:VCC+t1:DCC:CCC:ACC:VCC+DCC:CCC:ACC+DCC:CCC:VCC+DCC:CCC:ACC:VCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItem45)
anova(mT45.I,mT45.I.l)#25.921      1  3.557e-07 ***

### only Adults
byItemA<-byItem[byItem$AgeGroup=="Adult",]
## since these analyses compare Phon to Sem, use DET
byItemA$DC<-ifelse(byItemA$DET=="Incongruent",-.5,.5)
byItemA$DCC<-scale(byItemA$DC, T, F)
## cannot add vocabulary
#bySubj3$VCC<-scale(bySubj3$BPVS, T, F)
# Condition
byItemA$CC<-ifelse(byItemA$Condition=="Sem",-.5,.5)
byItemA$CCC<-scale(byItemA$CC, T, F)

mTA.I<-lmer(elogT~1+t1*DCC*CCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItemA)
summary(mTA.I)
confint(mTA.I, method="Wald")
#loglik
mTA.I.i<-lmer(elogT~1+t1*DCC+t1*CCC+t1:DCC:CCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItemA)
anova(mTA.I.i,mTA.I)#19.533      1   9.89e-06 ***
mTA.I.l<-lmer(elogT~1+t1*DCC+t1*CCC+DCC:CCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:CCC)+(1+t1||ItemN:DCC:CCC), data=byItemA)
anova(mTA.I.l,mTA.I)#10.571      1   0.001149 **

## Sem
byItemS<-byItem[byItem$Condition=="Sem",]
## since these analyses look at Sem only use DetS
byItemS$DC<-ifelse(byItemS$DetS=="One",-.5,.5)
byItemS$DCC<-scale(byItemS$DC, T, F)
## cannot add add vocabulary as a continuous, centred covariate because this analysis includes adults
#gcaT$VCC<-scale(gcaT$BPVS, T, F)
## Include Age using backward difference coding (backward because we compare each level to previous one)
# 3yo vs. 2yo
byItemS$A32<-ifelse(byItemS$AgeGroup=="2yo",-3/4,1/4)
# 45yo vs. 3yo
byItemS$A453<-ifelse(byItemS$AgeGroup=="2yo",-1/2,ifelse(byItemS$AgeGroup=="3yo",-1/2,1/2))
# Adult vs. 45yo
byItemS$AA45<-ifelse(byItemS$AgeGroup=="Adult",3/4,-1/4)

# Full model with Age contrasts
mS.I<-lmer(elogT~1+t1*DCC*(A32+A453+AA45)+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemS)
summary(mS.I)
confint(mS.I,method="Wald")
mS.I.i.32<-lmer(elogT~1+t1*DCC*(A453+AA45)+A32+t1:A32+t1:DCC:A32+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemS)
anova(mS.I,mS.I.i.32) #7.8306      1   0.005137 **
mS.I.i.A45<-lmer(elogT~1+t1*DCC*(A453+A32)+AA45+t1:AA45+t1:DCC:AA45+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemS)
anova(mS.I,mS.I.i.A45) #40.034      1  2.495e-10 ***
mS.I.l.32<-lmer(elogT~1+t1*DCC*(A453+AA45)+A32+t1:A32+DCC:A32+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemS)#did not converge!
anova(mS.I,mS.I.l.32) #12.615      1  0.0003827 *** (this did not return any warnings)
mS.I.l.A45<-lmer(elogT~1+t1*DCC*(A453+A32)+AA45+t1:AA45+DCC:AA45+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemS)
anova(mS.I,mS.I.l.A45) #131.68      1  < 2.2e-16 ***

### now analyse each age group separately and add vocabulary in children analyses

## only 2 yo
byItemS2<-byItemS[byItemS$AgeGroup=="2yo",]
## since these analyses look at Sem only use DetS
byItemS2$DC<-ifelse(byItemS2$DetS=="One",-.5,.5)
byItemS2$DCC<-scale(byItemS2$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemS2$VCC<-scale(byItemS2$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemS2$ACC<-scale(byItemS2$Age, T, F)

#mS2.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC*ACC||ItemN)+(1+t1*VCC*ACC||ItemN:DCC), data=byItemS2)# did not converge
# try removing interaction of Voc and age from all random terms
mS2.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS2)
summary(mS2.I)
confint(mS2.I,method="Wald")
#LogLink
mS2.I.l<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS2)
anova(mS2.I.l,mS2.I)# 11.075      1  0.0008752 ***
mS2.I.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:ACC:VCC+t1:DCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS2)
anova(mS2.I.i,mS2.I)# 1.0649      1     0.3021

## only 3 yo
byItemS3<-byItemS[byItemS$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
byItemS3$DC<-ifelse(byItemS3$DetS=="One",-.5,.5)
byItemS3$DCC<-scale(byItemS3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemS3$VCC<-scale(byItemS3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemS3$ACC<-scale(byItemS3$Age, T, F)

#mS3.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC*ACC||ItemN)+(1+t1*VCC*ACC||ItemN:DCC), data=byItemS3)# did not converge
# remove interaction of VCC and ACC from all random terms
mS3.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS3)
summary(mS3.I)
confint(mS3.I,method="Wald")
#LogLink
mS3.I.l<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS3)
anova(mS3.I.l,mS3.I)# 25.638      1  4.119e-07 ***
mS3.I.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:ACC:VCC+t1:DCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS3)
anova(mS3.I.i,mS3.I)# 25.538      1  4.337e-07 ***
mS3.I.l.voc<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS3)
anova(mS3.I.l.voc,mS3.I)#11.695      1  0.0006266 ***

## only 45 yo
byItemS45<-byItemS[byItemS$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
byItemS45$DC<-ifelse(byItemS45$DetS=="One",-.5,.5)
byItemS45$DCC<-scale(byItemS45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemS45$VCC<-scale(byItemS45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemS45$ACC<-scale(byItemS45$Age, T, F)

#mS45.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC*ACC||ItemN)+(1+t1*VCC*ACC||ItemN:DCC), data=byItemS45)# did not converge
# remove voc by age int from all random effects
mS45.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS45)
summary(mS45.I)
confint(mS45.I,method="Wald")
## loglik
mS45.I.l<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS45)
anova(mS45.I.l,mS45.I)#31.487      1  2.008e-08 ***
mS45.I.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC:ACC+t1:DCC:VCC+t1:DCC:ACC:VCC+t1:DCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS45)
anova(mS45.I.i,mS45.I)#26.833      1  2.218e-07 ***
mS45.I.l.voc<-lmer(elogT~1+t1*VCC*ACC+DCC*VCC*ACC+t1:DCC:ACC+t1:DCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemS45)
anova(mS45.I.l.voc,mS45.I)#6.0204      1    0.01414 *

## only Adults
byItemSA<-byItemS[byItemS$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
byItemSA$DC<-ifelse(byItemSA$DetS=="One",-.5,.5)
byItemSA$DCC<-scale(byItemSA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)

mSA.I<-lmer(elogT~1+t1*DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=byItemSA)
summary(mSA.I)
#logLik
mSA.I.l<-lmer(elogT~1+t1+DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=byItemSA)
anova(mSA.I,mSA.I.l)#54.731      1  1.382e-13 ***
mSA.I.i<-lmer(elogT~1+t1:DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=byItemSA)
anova(mSA.I,mSA.I.i)#52.628      2  3.733e-12 ***

## Phon
byItemP<-byItem[byItem$Condition=="Phon",]
## since these analyses look at Phon only use DetP
byItemP$DC<-ifelse(byItemP$DetP=="A",-.5,.5)
byItemP$DCC<-scale(byItemP$DC, T, F)
## cannot add add vocabulary as a continuous, centred covariate because this analysis includes adults
#gcaT$VCC<-scale(gcaT$BPVS, T, F)
## Include Age using backward difference coding (backward because we compare each level to previous one)
# 3yo vs. 2yo
byItemP$A32<-ifelse(byItemP$AgeGroup=="2yo",-3/4,1/4)
# 45yo vs. 3yo
byItemP$A453<-ifelse(byItemP$AgeGroup=="2yo",-1/2,ifelse(byItemP$AgeGroup=="3yo",-1/2,1/2))
# Adult vs. 45yo
byItemP$AA45<-ifelse(byItemP$AgeGroup=="Adult",3/4,-1/4)

mP.I<-lmer(elogT~1+t1*DCC*(A32+A453+AA45)+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemP)
summary(mP.I)
confint(mP.I,method="Wald")
mP.I.iA45<-lmer(elogT~1+t1*DCC*(A32+A453)+t1:AA45+t1:DCC:AA45+AA45+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemP)
anova(mP.I,mP.I.iA45)#14.606      1  0.0001325 ***
mP.I.lA45<-lmer(elogT~1+t1*DCC*(A32+A453)+t1:AA45+DCC:AA45+AA45+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC), data=byItemP)
anova(mP.I,mP.I.lA45)#137.77      1  < 2.2e-16 ***

### now analyse each age group separately and add vocabulary in children analyses

## only 2 yo
byItemP2<-byItemP[byItemP$AgeGroup=="2yo",]
## since these analyses look at Sem only use DetS
byItemP2$DC<-ifelse(byItemP2$DetP=="A",-.5,.5)
byItemP2$DCC<-scale(byItemP2$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemP2$VCC<-scale(byItemP2$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemP2$ACC<-scale(byItemP2$Age, T, F)

#mP2.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC*ACC||ItemN)+(1+t1*VCC*ACC||ItemN:DCC), data=byItemP2)# did not converge
# remove age by voc int from all random terms
mP2.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP2)
summary(mP2.I)
confint(mP2.I, method="Wald")
##loglik
mP2.I.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC+t1:DCC:VCC+t1:DCC:ACC+t1:DCC:ACC:VCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP2)
anova(mP2.I,mP2.I.i)#0.0716      1      0.789
mP2.I.l<-lmer(elogT~1+t1*VCC*ACC+DCC+t1:DCC:VCC+t1:DCC:ACC+t1:DCC:ACC:VCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP2)
anova(mP2.I,mP2.I.l)#0.2638      1     0.6075
mP2.I.i.Voc<-lmer(elogT~1+t1*VCC*ACC+t1:DCC+t1:DCC:VCC+t1:DCC:ACC+t1:DCC:ACC:VCC+DCC:ACC+DCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP2)
anova(mP2.I,mP2.I.i.Voc)#8.4896      1   0.003572 **

## only 3 yo
byItemP3<-byItemP[byItemP$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
byItemP3$DC<-ifelse(byItemP3$DetP=="A",-.5,.5)
byItemP3$DCC<-scale(byItemP3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemP3$VCC<-scale(byItemP3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemP3$ACC<-scale(byItemP3$Age, T, F)

#mP3.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC*ACC||ItemN)+(1+t1*VCC*ACC||ItemN:DCC), data=byItemP3)# did not converge
# remove voc by age int 
mP3.I<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP3)
summary(mP3.I)
confint(mP3.I,method="Wald")
##loglik
mP3.I.i<-lmer(elogT~1+t1*VCC*ACC+t1:DCC+t1:DCC:VCC+t1:DCC:ACC+t1:DCC:ACC:VCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP3)
anova(mP3.I,mP3.I.i)#11.492      1  0.0006991 ***
mP3.I.l<-lmer(elogT~1+t1*VCC*ACC+DCC+t1:DCC:VCC+t1:DCC:ACC+t1:DCC:ACC:VCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP3)
anova(mP3.I,mP3.I.l)#0.0824      1      0.774
mP3.I.l.voc<-lmer(elogT~1+t1*VCC*ACC+DCC+t1:DCC+t1:DCC:ACC+t1:DCC:ACC:VCC+DCC:ACC+DCC:VCC+DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP3)
anova(mP3.I,mP3.I.l.voc)#4.7533      1    0.02924 *

## only 45 yo
byItemP45<-byItemP[byItemP$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
byItemP45$DC<-ifelse(byItemP45$DetP=="A",-.5,.5)
byItemP45$DCC<-scale(byItemP45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemP45$VCC<-scale(byItemP45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemP45$ACC<-scale(byItemP45$Age, T, F)

# as in subject analyses, this model does not include DCC:ACC:VCC in the fixed effect structure
mP45.I2<-lmer(elogT~1+t1*DCC*VCC+t1*DCC*ACC+ACC:VCC+t1:ACC:VCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP45)
summary(mP45.I2)
confint(mP45.I2,method="Wald")
#LogLik
mP45.I2.i<-lmer(elogT~1+t1+t1:DCC+t1:VCC+DCC:VCC+t1:DCC:VCC+VCC+t1:ACC+DCC:ACC+t1:DCC:ACC+ACC+ACC:VCC+t1:ACC:VCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP45)
anova(mP45.I2,mP45.I2.i)#5.1875      1    0.02275 *
mP45.I2.l<-lmer(elogT~1+t1+DCC+t1:VCC+DCC:VCC+t1:DCC:VCC+VCC+t1:ACC+DCC:ACC+t1:DCC:ACC+ACC+ACC:VCC+t1:ACC:VCC+t1:DCC:ACC:VCC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=byItemP45)
anova(mP45.I2,mP45.I2.l)#0.0033      1     0.9539

## only Adults
byItemPA<-byItemP[byItemP$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
byItemPA$DC<-ifelse(byItemPA$DetP=="A",-.5,.5)
byItemPA$DCC<-scale(byItemPA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)

mPA.I<-lmer(elogT~1+t1*DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=byItemPA)
summary(mPA.I)
confint(mPA.I,method="Wald")
#loglik
mPA.I.i<-lmer(elogT~1+t1+t1:DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=byItemPA)
anova(mPA.I.i,mPA.I)#17.299      1  3.194e-05 ***
mPA.I.l<-lmer(elogT~1+t1+DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=byItemPA)
anova(mPA.I.l,mPA.I)#24.355      1   8.01e-07 ***