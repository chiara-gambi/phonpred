# Gambi, Gorrie, Pickering, Rabagliati (2017). The development of linguistic prediction: Predictions of sound and meaning in 2-to-5 year olds.
# Time to first look to target on Distractor initial trials (from Determiner Onset)
# Experiment 2
fixsumbinT.D<-read.table("First_fix_Det.txt",header=T)
##means by Onset and OldNew and Agecat
summaryBy(FixStart.rel~Onset+OldNew+Agecat,data=fixsumbinT.D)
summaryBy(FixStart.rel~Onset+Agecat,data=fixsumbinT.D)
summaryBy(FixStart.rel~OldNew+Agecat,data=fixsumbinT.D)
summaryBy(FixStart.rel~Onset+OldNew,data=fixsumbinT.D)

# Onset OldNew    Agecat FixStart.rel.mean
# 1 Different      N Four-Five          957.3213
# 2 Different      N     Three         1081.8827
# 3 Different      Y Four-Five          926.8353
# 4 Different      Y     Three         1049.9537
# 5      Same      N Four-Five          994.8627
# 6      Same      N     Three         1098.6726
# 7      Same      Y Four-Five          912.2579
# 8      Same      Y     Three          921.7263


# plot and inspect for outliers and non-normality
library(lattice)
densityplot(~FixStart.rel, data=fixsumbinT.D) # no clear outliers, but very much a right tail, bimodal

#sqaure root transformation
densityplot(~sqrt(FixStart.rel), data=fixsumbinT.D) # no clear outliers, but bimodal
# log transformation
densityplot(~log(FixStart.rel), data=fixsumbinT.D) # no clear outliers, but very much a left tail, possibly bimodal
#inverse transformation
densityplot(~1/(FixStart.rel), data=fixsumbinT.D) # no clear outliers, but very much a right tail - at least does not look bimodal

#use inverse transformation below?

# code Onset/OldNew
fixsumbinT.D$DC<-ifelse(fixsumbinT.D$Onset=="Same",-.5,.5)
fixsumbinT.D$DCC<-scale(fixsumbinT.D$DC, T, F)
fixsumbinT.D$OC<-ifelse(fixsumbinT.D$OldNew=="Y",-.5,.5)
fixsumbinT.D$OCC<-scale(fixsumbinT.D$OC, T, F)
# code age
fixsumbinT.D$AC<-ifelse(fixsumbinT.D$Agecat=="Three",-.5,.5)
fixsumbinT.D$ACC<-scale(fixsumbinT.D$AC, T, F)
summary(fixsumbinT.D)
library(lme4)
m.ttff<-lmer(FixStart.rel~1+DCC*OCC*ACC+(1+DCC*OCC||Participant)+(1+DCC*OCC||ItemN), data=fixsumbinT.D,REML=F)
summary(m.ttff)

# try with inverse transformation
fixsumbinT.D$FixStart.rel.inv<-1/fixsumbinT.D$FixStart.rel
m.ttff.inv<-lmer(FixStart.rel.inv~1+DCC*OCC*ACC+(1+DCC*OCC||Participant)+(1+DCC*OCC||ItemN), data=fixsumbinT.D,REML=F)
summary(m.ttff.inv)

## only new
m.ttff.N<-lmer(FixStart.rel~1+DCC*ACC+(1+DCC||Participant)+(1+DCC||ItemN), data=subset(fixsumbinT.D,OldNew=="N"),REML=F)
summary(m.ttff.N)
confint(m.ttff.N, method="Wald")
## only old
m.ttff.O<-lmer(FixStart.rel~1+DCC*ACC+(1+DCC||Participant)+(1+DCC||ItemN), data=subset(fixsumbinT.D,OldNew=="Y"),REML=F)
summary(m.ttff.O)
confint(m.ttff.O, method="Wald")

# try with inverse transformation
## only new
m.ttff.inv.N<-lmer(FixStart.rel.inv~1+DCC*ACC+(1+DCC||Participant)+(1+DCC||ItemN), data=subset(fixsumbinT.D,OldNew=="N"),REML=F)
summary(m.ttff.inv.N)
## only old
m.ttff.inv.O<-lmer(FixStart.rel.inv~1+DCC*ACC+(1+DCC||Participant)+(1+DCC||ItemN), data=subset(fixsumbinT.D,OldNew=="Y"),REML=F)
summary(m.ttff.inv.O)
