# Gambi, Gorrie, Pickering, Rabagliati (2017). The development of linguistic prediction: Predictions of sound and meaning in 2-to-5 year olds.
# Post-hoc analyses looking at the effect of picture occurrence
# These analyses are not reported in the paper
# Experiment 1

################### Statistical Analyses ##################
###################### Plan of analyses ###################
### We are going to do the following
### Elog Growth Curve Analysis both by subjects and items
### We do this separately for Phon and Sem
### We compare each age-group to the age group immediately younger, with adults as the reference level
### We include Vocabulary as a continuous predictor
### NOTE: All analyses split by occurrence (App == 1 vs. 2) below.

# ELOG BY SUBJECTS
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
### NOTE: we run no APP analyses comparing across Phon and Sem, as the purpose of these analyses is to ascertain whether
### the Det effect was present at the first occurrence in both Phon and Sem (rather than to judge if it was larger in Sem than Phon)

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
## App
summary(bySubjS$App)# numeric, 1 or 2
bySubjS$AppC<-ifelse(bySubjS$App==1,-.5,.5)
bySubjS$AppCC<-scale(bySubjS$AppC, T, F)

# Full model with Age contrasts (note that we had to expand the random structure as well as adding AppCC to fixed effects)
mS.App<-lmer(elogT~1+t1*DCC*AppCC*(A32+A453+AA45)+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjS)
summary(mS.App)

### now analyse each age group separately and add vocabulary in children analyses

## only 2 yo
bySubjS2<-bySubjS[bySubjS$AgeGroup=="2yo",]
## since these analyses look at Sem only use DetS
bySubjS2$DC<-ifelse(bySubjS2$DetS=="One",-.5,.5)
bySubjS2$DCC<-scale(bySubjS2$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjS2$VCC<-scale(bySubjS2$BPVS, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjS2$ACC<-scale(bySubjS2$Age, T, F)
## App
bySubjS2$AppC<-ifelse(bySubjS2$App==1,-.5,.5)
bySubjS2$AppCC<-scale(bySubjS2$AppC, T, F)

mS2.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjS2)
summary(mS2.App)

### only first occurrence
## App
mS2.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjS2,App==1))
summary(mS2.App.1)
#DCC             0.0565344  0.1526326   0.370
#t1:DCC          0.9060524  0.3947342   2.295
confint(mS2.App.1, method="Wald")
#DCC            -0.2426198908 0.355688755
#t1:DCC          0.1323877174 1.679717158

## only 3 yo
bySubjS3<-bySubjS[bySubjS$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
bySubjS3$DC<-ifelse(bySubjS3$DetS=="One",-.5,.5)
bySubjS3$DCC<-scale(bySubjS3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjS3$VCC<-scale(bySubjS3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjS3$ACC<-scale(bySubjS3$Age, T, F)
## App
bySubjS3$AppC<-ifelse(bySubjS3$App==1,-.5,.5)
bySubjS3$AppCC<-scale(bySubjS3$AppC, T, F)

# did not converge mS3.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjS3)
# remove VCC*ACC interactions from fixed part
mS3.App<-lmer(elogT~1+t1*DCC*AppCC*VCC+t1*DCC*AppCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjS3)
summary(mS3.App)

## only first occurrence
mS3.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjS3, App==1))
summary(mS3.App.1)
#DCC             0.4685949  0.1668463   2.809
#t1:DCC          2.1297460  0.4502651   4.730
confint(mS3.App.1,method="Wald")
#DCC             0.141582179 0.7956075765
#t1:DCC          1.247242721 3.0122493242

## only 45 yo
bySubjS45<-bySubjS[bySubjS$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
bySubjS45$DC<-ifelse(bySubjS45$DetS=="One",-.5,.5)
bySubjS45$DCC<-scale(bySubjS45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjS45$ACC<-scale(bySubjS45$Age, T, F)
## App
bySubjS45$AppC<-ifelse(bySubjS45$App==1,-.5,.5)
bySubjS45$AppCC<-scale(bySubjS45$AppC, T, F)

# did not converge mS45.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjS45)
# remove ACC*VCC interactions from fixed eggects
mS45.App<-lmer(elogT~1+t1*DCC*AppCC*VCC+t1*DCC*AppCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjS45)
summary(mS45.App)

#only first occurrence
mS45.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjS45, App==1))
summary(mS45.App.1)
#DCC             0.738690   0.119619   6.175
#t1:DCC          1.933714   0.441194   4.383
confint(mS45.App.1,method="Wald")
#DCC             0.504241163 0.973139541
#t1:DCC          1.068989981 2.798437044

## only Adults
bySubjSA<-bySubjS[bySubjS$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
bySubjSA$DC<-ifelse(bySubjSA$DetS=="One",-.5,.5)
bySubjSA$DCC<-scale(bySubjSA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)
## App
bySubjSA$AppC<-ifelse(bySubjSA$App==1,-.5,.5)
bySubjSA$AppCC<-scale(bySubjSA$AppC, T, F)

mSA.App<-lmer(elogT~1+t1*DCC*AppCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjSA)
summary(mSA.App)

#first occurrence only
mSA.App.1<-lmer(elogT~1+t1*DCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjSA, App==1))
summary(mSA.App.1)
#DCC          1.74454    0.17338  10.062
#t1:DCC       4.44027    0.44034  10.084
confint(mSA.App.1, method="Wald")
#DCC          1.40472469 2.0843517
#t1:DCC       3.57721460 5.3033307

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
## App
summary(bySubjP$App)# numeric, 1 or 2
bySubjP$AppC<-ifelse(bySubjP$App==1,-.5,.5)
bySubjP$AppCC<-scale(bySubjP$AppC, T, F)

# Full model with Age contrasts (note that we had to expand the random structure as well as adding AppCC to fixed effects)
mP.App<-lmer(elogT~1+t1*DCC*AppCC*(A32+A453+AA45)+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjP)
summary(mP.App)
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
## App
bySubjP2$AppC<-ifelse(bySubjP2$App==1,-.5,.5)
bySubjP2$AppCC<-scale(bySubjP2$AppC, T, F)

mP2.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjP2)
summary(mP2.App)

#only first occurrence
mP2.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjP2, App==1))
summary(mP2.App.1)
#DCC             0.1200711  0.1352987   0.888
#t1:DCC          0.2911351  0.4686292   0.621
#DCC:VCC         0.0306336  0.0125374   2.443
confint(mP2.App.1, method="Wald")
#DCC            -0.145109530 0.3852516559
#t1:DCC         -0.627361266 1.2096315177
#DCC:VCC         0.006060820 0.0552064506

#only second occurrence
mP2.App.2<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjP2, App==2))
summary(mP2.App.2)
#DCC            -0.132128   0.178768  -0.739
#t1:DCC          0.324366   0.339408   0.956
#0.032636   0.016105   2.026


## only 3 yo
bySubjP3<-bySubjP[bySubjP$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
bySubjP3$DC<-ifelse(bySubjP3$DetP=="A",-.5,.5)
bySubjP3$DCC<-scale(bySubjP3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjP3$VCC<-scale(bySubjP3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjP3$ACC<-scale(bySubjP3$Age, T, F)
## App
bySubjP3$AppC<-ifelse(bySubjP3$App==1,-.5,.5)
bySubjP3$AppCC<-scale(bySubjP3$AppC, T, F)

mP3.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjP3)
summary(mP3.App)

#only first occurrence
mP3.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjP3, App==1))
summary(mP3.App.1)
#DCC             0.172572   0.143498   1.203
#t1:DCC          0.147373   0.424745   0.347
confint(mP3.App.1, method="Wald")
#DCC            -0.108677803  0.453822403
#t1:DCC         -0.685112530  0.979858932

#only second occurrence
mP3.App.2<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjP3, App==2))
summary(mP3.App.2)
#DCC             0.305867   0.147625   2.072
#t1:DCC           0.297484   0.394330   0.754

## only 45 yo
bySubjP45<-bySubjP[bySubjP$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
bySubjP45$DC<-ifelse(bySubjP45$DetP=="A",-.5,.5)
bySubjP45$DCC<-scale(bySubjP45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
bySubjP45$VCC<-scale(bySubjP45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
bySubjP45$ACC<-scale(bySubjP45$Age, T, F)
## App
bySubjP45$AppC<-ifelse(bySubjP45$App==1,-.5,.5)
bySubjP45$AppCC<-scale(bySubjP45$AppC, T, F)

mP45.App<-lmer(elogT~1+t1*DCC*AppCC*VCC+1*DCC*AppCC*ACC+ACC:VCC+AppCC:ACC:VCC+t1:ACC:VCC+t1:DCC:ACC:VCC+t1:AppCC:ACC:VCC+t1:DCC:AppCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjP45)
summary(mP45.App)

# only first occurrence
mP45.App.1<-lmer(elogT~1+t1*DCC*VCC+t1*DCC*ACC+ACC:VCC+t1:ACC:VCC+t1:DCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjP45, App==1))
summary(mP45.App.1)
#DCC             0.219485   0.139912   1.569
#t1:DCC          0.034545   0.423366   0.082
confint(mP45.App.1,method="Wald")
#DCC            -0.0547375491 0.493706943
#t1:DCC         -0.7952368428 0.864327483

# only second occurrence
mP45.App.2<-lmer(elogT~1+t1*DCC*VCC+t1*DCC*ACC+ACC:VCC+t1:ACC:VCC+t1:DCC:ACC:VCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjP45, App==2))
summary(mP45.App.2)
#DCC             0.1991080  0.1294825   1.538
#t1:DCC          0.1940012  0.4757752   0.408


## only Adults
bySubjPA<-bySubjP[bySubjP$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
bySubjPA$DC<-ifelse(bySubjPA$DetP=="A",-.5,.5)
bySubjPA$DCC<-scale(bySubjPA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)
## App
bySubjPA$AppC<-ifelse(bySubjPA$App==1,-.5,.5)
bySubjPA$AppCC<-scale(bySubjPA$AppC, T, F)

mPA.App<-lmer(elogT~1+t1*DCC*AppCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:AppCC)+(1+t1||Participant:DCC:AppCC), data=bySubjPA)
summary(mPA.App)
##only first occurrence
mPA.App.1<-lmer(elogT~1+t1*DCC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(bySubjPA, App==1))
summary(mPA.App.1)
# DCC          0.66337    0.19565   3.391
# t1:DCC       2.30307    0.45911   5.016
confint(mPA.App.1, method="Wald")
# DCC          0.2799004 1.0468321
# t1:DCC       1.4032358 3.2028971

# ELOG BY ITEMS
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
### NOTE: we do not run these analyses (comparing across Sem and Phon) for post-hoc App checks.

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
#App
summary(byItemS$App)
byItemS$AppC<-ifelse(byItemS$App==1,-.5,.5)
byItemS$AppCC<-scale(byItemS$AppC, T, F)

# Full model with Age contrasts; added to random structure because age is within items
mS.I.App<-lmer(elogT~1+t1*DCC*AppCC*(A32+A453+AA45)+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC)+(1+t1+A32+A453+AA45||ItemN:AppCC)+(1+t1+A32+A453+AA45||ItemN:DCC:AppCC), data=byItemS)
summary(mS.I.App)
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
#App
byItemS2$AppC<-ifelse(byItemS2$App==1,-.5,.5)
byItemS2$AppCC<-scale(byItemS2$AppC, T, F)

# add Vocabulary and Age to random structure
# did not converge mS2.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1*VCC*ACC||ItemN)+(1+t1*VCC*ACC||ItemN:DCC)+(1+t1*VCC*ACC||ItemN:AppCC)+(1+t1*VCC*ACC||ItemN:DCC:AppCC), data=byItemS2)
# simplify random terms by removing interactions betwen ACC and VCC
# did not convergemS2.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC)+(1+t1*VCC+t1*ACC||ItemN:AppCC)+(1+t1*VCC+t1*ACC||ItemN:DCC:AppCC), data=byItemS2)
# simplify by removing all interactions from random structure
# did not converge mS2.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1+VCC+ACC||ItemN:DCC:AppCC), data=byItemS2)
# remove all random effects from most complex term
mS2.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1||ItemN:DCC:AppCC), data=byItemS2)
summary(mS2.I.App)

#only first occurrence
mS2.I.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=subset(byItemS2, App==1))
summary(mS2.I.App.1)
#DCC             0.0362101  0.1126498   0.321
#t1:DCC          0.6428526  0.2487494   2.584
confint(mS2.I.App.1,method="Wald")
#DCC            -0.184579415 0.2569995380
#t1:DCC          0.155312838 1.1303924505

## only 3 yo
byItemS3<-byItemS[byItemS$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
byItemS3$DC<-ifelse(byItemS3$DetS=="One",-.5,.5)
byItemS3$DCC<-scale(byItemS3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemS3$VCC<-scale(byItemS3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemS3$ACC<-scale(byItemS3$Age, T, F)
#App
byItemS3$AppC<-ifelse(byItemS3$App==1,-.5,.5)
byItemS3$AppCC<-scale(byItemS3$AppC, T, F)

# add Vocabulary and Age to random structure, but no int to aid convergence
mS3.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1+VCC+ACC||ItemN:DCC:AppCC), data=byItemS3)
summary(mS3.I.App)

#only first occurrence
# removed interaction of t1 and age from item random slopes to aid convergence
mS3.I.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1+ACC||ItemN:DCC), data=subset(byItemS3, App==1))
summary(mS3.I.App.1)
#DCC             0.2892480  0.0759588   3.808
#t1:DCC          1.3691327  0.2194664   6.238
confint(mS3.I.App.1, method="Wald")
#DCC             0.140371469  4.381246e-01
#t1:DCC          0.938986374  1.799279e+00

## only 45 yo
byItemS45<-byItemS[byItemS$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
byItemS45$DC<-ifelse(byItemS45$DetS=="One",-.5,.5)
byItemS45$DCC<-scale(byItemS45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemS45$VCC<-scale(byItemS45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemS45$ACC<-scale(byItemS45$Age, T, F)
#App
byItemS45$AppC<-ifelse(byItemS45$App==1,-.5,.5)
byItemS45$AppCC<-scale(byItemS45$AppC, T, F)
# first try with both voc and age in random structure but no interactions
# did not converge mS45.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1+VCC+ACC||ItemN:DCC:AppCC), data=byItemS45)
# remove both age and voc, starting from most complex random effect
# did not converge mS45.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1||ItemN:DCC:AppCC), data=byItemS45)
# so remove from all but most simple random term
mS45.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:AppCC)+(1+t1||ItemN:DCC:AppCC), data=byItemS45)
summary(mS45.I.App)
# reinstate Age and Voc effects on relevant random term to check this is not spurious
mS45.I.App2<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:AppCC)+(1+t1*ACC+VCC||ItemN:DCC:AppCC), data=byItemS45)
summary(mS45.I.App2)
#only first occurrence
mS45.I.App2.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=subset(byItemS45,App==1))
summary(mS45.I.App2.1)
#DCC             0.4626725  0.0510445   9.064
#t1:DCC          1.2817265  0.2597610   4.934
confint(mS45.I.App2.1, method="Wald")
#DCC             3.626271e-01  0.562717973
#t1:DCC          7.726043e-01  1.790848795
## only Adults
byItemSA<-byItemS[byItemS$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
byItemSA$DC<-ifelse(byItemSA$DetS=="One",-.5,.5)
byItemSA$DCC<-scale(byItemSA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)
#App
byItemSA$AppC<-ifelse(byItemSA$App==1,-.5,.5)
byItemSA$AppCC<-scale(byItemSA$AppC, T, F)

mSA.I.App<-lmer(elogT~1+t1*DCC*AppCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:AppCC)+(1+t1||ItemN:DCC:AppCC), data=byItemSA)
summary(mSA.I.App)

#only first occurrence
mSA.I.App.1<-lmer(elogT~1+t1*DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=subset(byItemSA, App==1))
summary(mSA.I.App.1)
# DCC           1.1475     0.1338   8.576
# t1:DCC        3.3000     0.2509  13.154
confint(mSA.I.App.1,method="Wald")
# DCC          0.88525259 1.40976623
# t1:DCC       2.80832127 3.79173554

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
#App
summary(byItemP$App)
byItemP$AppC<-ifelse(byItemP$App==1,-.5,.5)
byItemP$AppCC<-scale(byItemP$AppC, T, F)

# Full model with Age contrasts
mP.I.App<-lmer(elogT~1+t1*DCC*AppCC*(A32+A453+AA45)+(1+t1+A32+A453+AA45||ItemN)+(1+t1+A32+A453+AA45||ItemN:DCC)+(1+t1+A32+A453+AA45||ItemN:AppCC)+(1+t1+A32+A453+AA45||ItemN:DCC:AppCC), data=byItemP)
summary(mP.I.App)

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
#App
byItemP2$AppC<-ifelse(byItemP2$App==1,-.5,.5)
byItemP2$AppCC<-scale(byItemP2$AppC, T, F)

# add Vocabulary and age to random structure, but without interaction to increase convergence
mP2.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1+VCC+ACC||ItemN:DCC:AppCC), data=byItemP2)
summary(mP2.I.App)

mP2.I.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=subset(byItemP2, App==1))
summary(mP2.I.App.1)
#DCC            -0.0068240  0.1219654  -0.056
# DCC:VCC         0.0187840  0.0105216   1.785
# t1:DCC         -0.0450311  0.2882640  -0.156
# t1:DCC:VCC     -0.0178266  0.0185728  -0.960
confint(mP2.I.App.1, method="Wald")
#DCC            -0.245871807  0.2322237401
# DCC:VCC        -0.001837894  0.0394058433
# t1:DCC         -0.610018297  0.5199560009


## only 3 yo
byItemP3<-byItemP[byItemP$AgeGroup=="3yo",]
## since these analyses look at Sem only use DetS
byItemP3$DC<-ifelse(byItemP3$DetP=="A",-.5,.5)
byItemP3$DCC<-scale(byItemP3$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemP3$VCC<-scale(byItemP3$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemP3$ACC<-scale(byItemP3$Age, T, F)
#App
byItemP3$AppC<-ifelse(byItemP3$App==1,-.5,.5)
byItemP3$AppCC<-scale(byItemP3$AppC, T, F)

# add Vocabulary and age to random structure, but no int to aid convergence
mP3.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1+VCC+ACC||ItemN:DCC:AppCC), data=byItemP3)
summary(mP3.I.App)

# only first occurrence
mP3.I.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=subset(byItemP3, App==1))
summary(mP3.I.App.1)
#DCC             0.1520641  0.0640630   2.374
#t1:DCC          0.1211972  0.1780575   0.681
confint(mP3.I.App.1, method="Wald")
#DCC             0.0265028629  0.277625304
#t1:DCC         -0.2277891465  0.470183560


## only 45 yo
byItemP45<-byItemP[byItemP$AgeGroup=="4-5yo",]
## since these analyses look at Sem only use DetS
byItemP45$DC<-ifelse(byItemP45$DetP=="A",-.5,.5)
byItemP45$DCC<-scale(byItemP45$DC, T, F)
## add vocabulary as a continuous, centred covariate 
byItemP45$VCC<-scale(byItemP45$BPVS, T, F)
## add age in months as a continuous, centred covariate 
byItemP45$ACC<-scale(byItemP45$Age, T, F)
#App
byItemP45$AppC<-ifelse(byItemP45$App==1,-.5,.5)
byItemP45$AppCC<-scale(byItemP45$AppC, T, F)

# did not converge mP45.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1+VCC+ACC||ItemN:DCC:AppCC), data=byItemP45)
# remove voc and age from most complex random term
mP45.I.App<-lmer(elogT~1+t1*DCC*AppCC*VCC*ACC+(1+t1+VCC+ACC||ItemN)+(1+t1+VCC+ACC||ItemN:DCC)+(1+t1+VCC+ACC||ItemN:AppCC)+(1+t1||ItemN:DCC:AppCC), data=byItemP45)
summary(mP45.I.App)

# only first occurrence
mP45.I.App.1<-lmer(elogT~1+t1*DCC*VCC*ACC+(1+t1*VCC+t1*ACC||ItemN)+(1+t1*VCC+t1*ACC||ItemN:DCC), data=subset(byItemP45, App==1))
summary(mP45.I.App.1)
#DCC             0.1367107  0.1004821   1.361
#t1:DCC          0.0284112  0.2187222   0.130
confint(mP45.I.App.1,method="Wald")
#DCC            -0.0602306250 0.3336521014
#t1:DCC         -0.4002764280 0.4570988815

## only Adults
byItemPA<-byItemP[byItemP$AgeGroup=="Adult",]
## since these analyses look at Sem only use DetS
byItemPA$DC<-ifelse(byItemPA$DetP=="A",-.5,.5)
byItemPA$DCC<-scale(byItemPA$DC, T, F)
## cannot add vocabulary as a continuous, centred covariate 
#bySubjS45$VCC<-scale(bySubjS45$BPVS, T, F)
#App
byItemPA$AppC<-ifelse(byItemPA$App==1,-.5,.5)
byItemPA$AppCC<-scale(byItemPA$AppC, T, F)

mPA.I.App<-lmer(elogT~1+t1*DCC*AppCC+(1+t1||ItemN)+(1+t1||ItemN:DCC)+(1+t1||ItemN:AppCC)+(1+t1||ItemN:DCC:AppCC), data=byItemPA)
summary(mPA.I.App)

# only first occurrence
mPA.I.App.1<-lmer(elogT~1+t1*DCC+(1+t1||ItemN)+(1+t1||ItemN:DCC), data=subset(byItemPA, App==1))
summary(mPA.I.App.1)
# DCC          0.31869    0.15613   2.041
# t1:DCC       1.76009    0.38112   4.618
confint(mPA.I.App.1, method="Wald")
# DCC          0.01267067 0.6247023
# t1:DCC       1.01310093 2.5070758