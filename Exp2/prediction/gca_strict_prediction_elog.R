# Gambi, Gorrie, Pickering, Rabagliati (2017). The development of linguistic prediction: Predictions of sound and meaning in 2-to-5 year olds.
# Analysis of fixation proportions from Determiner Onset to Noun Onset + 100ms (Prediction window) using Growth Curves
# Experiment 2

#################### Start of statistical analyses ##################
# Effect of Onset, separately for New and Old Trials
# Effect of Trial Type (Old vs. New)

gcaVTP.elog.subj<-read.table("elog_Target_bysubj_GCA_DONO_strict.txt", header=T)
gcaVTP.elog.items<-read.table("elog_Target_byitem_GCA_DONO_strict.txt", header=T)

####################### Subjects ####################################
t<-poly(unique(gcaVTP.elog.subj$time),1)
time<-as.vector(unique(gcaVTP.elog.subj$time))
t<-cbind(t,time)
t<-as.data.frame(t)
gcaVTP.elog.subj<-gcaVTP.elog.subj[order(gcaVTP.elog.subj$time),]
gcaVTP.elog.subj$t1<-NA
#gcaP$t2<-NA
for (i in (1:nrow(gcaVTP.elog.subj))){
  gcaVTP.elog.subj$t1[i]<-t[t$time==gcaVTP.elog.subj$time[i],1] 
  #gcaP$t2[i]<-t[t$time==gcaP$time[i],2] 
}
gcaVTP.elog.subj$DC<-ifelse(gcaVTP.elog.subj$Onset=="Same",-.5,.5)
gcaVTP.elog.subj$DCC<-scale(gcaVTP.elog.subj$DC, T, F)
gcaVTP.elog.subj$OC<-ifelse(gcaVTP.elog.subj$OldNew=="Y",-.5,.5)
gcaVTP.elog.subj$OCC<-scale(gcaVTP.elog.subj$OC, T, F)
gcaVTP.elog.subj$AC<-ifelse(gcaVTP.elog.subj$Agecat=="Three",-.5,.5)
gcaVTP.elog.subj$ACC<-scale(gcaVTP.elog.subj$AC, T, F)
library(lme4)

mTP.elog.subj<-lmer(elogT~1+(t1)*DCC*OCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:OCC), data=gcaVTP.elog.subj)
summary(mTP.elog.subj)
confint(mTP.elog.subj, method="Wald")
#OCC            -0.55877473 -0.36109037
#t1:OCC         -0.41049563  0.06205345
#logLik
mTP.elog.subj.noOCC<-lmer(elogT~1+(t1)*DCC*ACC+OCC:DCC+OCC:DCC:t1+OCC:ACC+OCC:ACC:t1+OCC:DCC:ACC+OCC:DCC:ACC:t1+t1:OCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:OCC), data=gcaVTP.elog.subj)
anova(mTP.elog.subj.noOCC,mTP.elog.subj)#64.815      1  8.226e-16 ***
mTP.elog.subj.noOCC.lin<-lmer(elogT~1+(t1)*DCC*ACC+OCC:DCC+OCC:DCC:t1+OCC:ACC+OCC:ACC:t1+OCC:DCC:ACC+OCC:DCC:ACC:t1+OCC+(1+t1||Participant)+(1+t1||Participant:DCC)+(1+t1||Participant:OCC), data=gcaVTP.elog.subj)
anova(mTP.elog.subj.noOCC.lin,mTP.elog.subj)#2.1217      1     0.1452

## only new
mTP.elog.subj.N<-lmer(elogT~1+(t1)*DCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(gcaVTP.elog.subj, OldNew=="N"))
summary(mTP.elog.subj.N)
confint(mTP.elog.subj.N, method="Wald")
mTP.elog.subj.N.i<-lmer(elogT~1+(t1)*ACC+t1:DCC+DCC:ACC+t1:DCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(gcaVTP.elog.subj, OldNew=="N"))
anova(mTP.elog.subj.N,mTP.elog.subj.N.i)#4.6924      1     0.0303 *
mTP.elog.subj.N.l<-lmer(elogT~1+(t1)*ACC+DCC+DCC:ACC+t1:DCC:ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(gcaVTP.elog.subj, OldNew=="N"))
anova(mTP.elog.subj.N,mTP.elog.subj.N.l)#0.1321      1     0.7163

## only old
mTP.elog.subj.O<-lmer(elogT~1+(t1)*DCC*ACC+(1+t1||Participant)+(1+t1||Participant:DCC), data=subset(gcaVTP.elog.subj, OldNew=="Y"))
summary(mTP.elog.subj.O)
confint(mTP.elog.subj.O,method="Wald")

####################### Items ####################################
t<-poly(unique(gcaVTP.elog.items$time),1)
time<-as.vector(unique(gcaVTP.elog.items$time))
t<-cbind(t,time)
t<-as.data.frame(t)
gcaVTP.elog.items<-gcaVTP.elog.items[order(gcaVTP.elog.items$time),]
gcaVTP.elog.items$t1<-NA
#gcaP$t2<-NA
for (i in (1:nrow(gcaVTP.elog.items))){
  gcaVTP.elog.items$t1[i]<-t[t$time==gcaVTP.elog.items$time[i],1] 
  #gcaP$t2[i]<-t[t$time==gcaP$time[i],2] 
}
gcaVTP.elog.items$DC<-ifelse(gcaVTP.elog.items$Onset=="Same",-.5,.5)
gcaVTP.elog.items$DCC<-scale(gcaVTP.elog.items$DC, T, F)
gcaVTP.elog.items$OC<-ifelse(gcaVTP.elog.items$OldNew=="Y",-.5,.5)
gcaVTP.elog.items$OCC<-scale(gcaVTP.elog.items$OC, T, F)
gcaVTP.elog.items$AC<-ifelse(gcaVTP.elog.items$Agecat=="Three",-.5,.5)
gcaVTP.elog.items$ACC<-scale(gcaVTP.elog.items$AC, T, F)
library(lme4)

#include agecat
mTP.elog.items<-lmer(elogT~1+(t1)*DCC*OCC*ACC+(1+t1*ACC||ItemN)+(1+t1*ACC||ItemN:DCC)+(1+t1*ACC||ItemN:OCC), data=gcaVTP.elog.items)
summary(mTP.elog.items)
confint(mTP.elog.items, method="Wald")
#OCC            -0.62739022 -0.33638060
#loglik
mTP.elog.items.noOCC<-lmer(elogT~1+(t1)*DCC*ACC+t1:OCC+DCC:OCC+DCC:OCC:t1+ACC:OCC+ACC:OCC:t1+OCC:DCC:ACC+OCC:DCC:ACC:t1+(1+t1*ACC||ItemN)+(1+t1*ACC||ItemN:DCC)+(1+t1*ACC||ItemN:OCC), data=gcaVTP.elog.items)
anova(mTP.elog.items.noOCC,mTP.elog.items)#23.417      1  1.304e-06 ***

#only new
mTP.elog.items.N<-lmer(elogT~1+(t1)*DCC*ACC+(1+t1*ACC||ItemN)+(1+t1*ACC||ItemN:DCC), data=subset(gcaVTP.elog.items,OldNew=="N"))
summary(mTP.elog.items.N)
confint(mTP.elog.items.N, method="Wald")
#log lik
mTP.elog.items.N.i<-lmer(elogT~1+(t1)*ACC+t1:DCC+DCC:ACC+t1:DCC:ACC+(1+t1*ACC||ItemN)+(1+t1*ACC||ItemN:DCC), data=subset(gcaVTP.elog.items,OldNew=="N"))
anova(mTP.elog.items.N,mTP.elog.items.N.i)#0.6554      1     0.4182
mTP.elog.items.N.l<-lmer(elogT~1+(t1)*ACC+DCC+DCC:ACC+t1:DCC:ACC+(1+t1*ACC||ItemN)+(1+t1*ACC||ItemN:DCC), data=subset(gcaVTP.elog.items,OldNew=="N"))
anova(mTP.elog.items.N,mTP.elog.items.N.l)#2e-04      1     0.9876
#plot effect by items
byitems<-data.frame(ItemN=unique(gcaVTP.elog.items$ItemN),OnsetNew=NA)
Onsetbyitem<-fixef(mTP.elog.items.N)[3]+ranef(mTP.elog.items.N)$"ItemN:DCC"[rep(c(FALSE,TRUE),20),4]-ranef(mTP.elog.items.N)$"ItemN:DCC"[rep(c(TRUE,FALSE),20),4]
byitems$OnsetNew<-Onsetbyitem
plot(byitems$ItemN,byitems$OnsetNew)

#only old
mTP.elog.items.O<-lmer(elogT~1+(t1)*DCC*ACC+(1+t1*ACC||ItemN)+(1+t1*ACC||ItemN:DCC), data=subset(gcaVTP.elog.items,OldNew=="Y"))
summary(mTP.elog.items.O)
confint(mTP.elog.items.O,method="Wald")