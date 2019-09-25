library(ggplot2)
library(plyr)
library(gridExtra)
library(tableone)
library(magrittr)
library(dplyr)
library(ggseg)
library(tidyr)
library(psych)
library(lme4)
library(ppcor)
library(lemon)
library(data.table)
library(lmtest)
library(plm)
library(grid)
library(lattice)
##################################################################################################################
##################################################################################################################
#Functions
#source("W:/ances/julie/Data/StarkRequest/CommonFuncs.r")
source("C:/Users/julie.wisch/Documents/Stark/CommonFuncs.r")
source("C:/Users/julie.wisch/Documents/SexDiffs/SexDiffFuncs.r")
source("C:/Users/julie.wisch/Documents/SexDiffs/SexDiffDataCleaning_DR15.r")
##################################################################################################################
##################################################################################################################

df<-merge(df, Hippo.all[,c("ID", "VolRateofChange")], by = "ID", all = FALSE)


#############################################################################################################################
#Max Likelihood Model Selection to test for age-sex-amyloid differences
df.log<-df
#log transforming everything since it's so left heavy
for(i in c(3:96, 98:191)){
  df.log[,i]<-log(df[,i])
}



################################################       
df.log<-df.log[,-c(18, 18+95)]   
df.log<-df.log[,-c(55, 55+94)]
df.log<-df.log[,-c(62, 62+93)]
df.log<-df.log[, colSums(df.log != 0) > 0]
df.log$apoe4<-as.factor(df.log$apoe4)
result<-data.frame("Region" = character(), "Amyloid" = numeric(), "Sex" = numeric(), 
                   "APOE" = numeric(), "Interaction" = numeric(), "Sex_noAmyloid" = numeric(),
                   "APOE_noAmyloid" = numeric(), "Interaction_noAmyloid" = numeric(),
                   "Tau_Cortmean" = numeric(),
                   "Sex_Cortmean" = numeric(),
                   "APOE_Cortmean" = numeric(), "Interaction_Cortmean" = numeric(),
                   stringsAsFactors = FALSE)
k<-0

for(i in c(3:48)){
  k<-k+1
  result[k, 1]<-substring(names(df.log[i]), 15)
  
  Model.null<-lm(df.log[,i] ~ df.log[,i+92] + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
  result_hold<-coeftest(Model.null, vcov = vcovHC(Model.null))
  result[k, 2]<-result_hold[2, 4] #amyloid
  result[k, 3]<-result_hold[4, 4] #sex
  result[k, 4]<-result_hold[5, 4]#APOE
  result[k, 5]<-result_hold[7, 4] #interaction
  
  
  #No longer controlling for amyloid levels when doing the test...
  Model.null<-lm(df.log[,i] ~  Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
  result_hold<-coeftest(Model.null, vcov = vcovHC(Model.null))
  result[k, 6]<-result_hold[3, 4]
  result[k, 7]<-result_hold[4, 4]
  result[k, 8]<-result_hold[6, 4]
 
  
  
  #Correcting for global amyloid burden
  Model.null<-lm(df.log[,i] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
  result_hold<-coeftest(Model.null, vcov = vcovHC(Model.null))
  result[k, 9]<-result_hold[2, 4]
  result[k, 10]<-result_hold[4, 4]
  result[k, 11]<-result_hold[5, 4]
  result[k, 12]<-result_hold[7, 4]
}
rm(Model.new, Model.null, i, k)
result<-result[complete.cases(result$Amyloid),]
#Significant after bonferroni correction
result[result$Sex_Cortmean < 0.05/45,c(1, 3)]


Model.null<-lm(Tauopathy ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
coeftest(Model.null, vcov = vcovHC(Model.null)) #males have significantly lower tau than females for overall tauopathy p = 0.0249

Model.null<-lm(Tauopathy ~  Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
coeftest(Model.null, vcov = vcovHC(Model.null)) #sex difference persists even if you don't control for amyloid

#Sex differences exist if you don't correct for amyloid, if you correct for regional amyloid, or if you correct for cortical amyloid total
Model1<-lm(df.log[,17] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
Model2<-lm(df.log[,19] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
Model3<-lm(df.log[,25] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)#males more
Model4<-lm(df.log[,33] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
Model5<-lm(df.log[,34] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
Model6<-lm(df.log[,41] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
Model7<-lm(df.log[,44] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
Model8<-lm(Tauopathy ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + apoe4 + EDUC + apoe4:GENDER, data = df.log)
#caudal middle frontal, cuneus, insula, parahippocampus, pars opercularis, precuneus, superior temporal
library(stargazer)
stargazer(Model1, Model2, Model3, Model4, Model5, Model6, Model7, Model8, title = "Tau Differs Between Males and Females",
          align= TRUE, type = "html", out = "C:/Users/julie.wisch/Documents/SexDiffs/table_v5.htm",
          dep.var.labels = c("Cuneus", "Frontal Pole", "Lateral Occipital", "Pars Orbitalis",
                             "Pars Triangularis", "Rostral Middle Frontal", "Superior Temporal", "Global Summary"),
          omit.stat=c("LL","ser"),
          ci = FALSE, star.cutoffs = c(0.05/45, 0.01/45, 0.001/45))


coeftest(Model1, vcov = vcovHC(Model1)) 
coeftest(Model2, vcov = vcovHC(Model2)) 
coeftest(Model3, vcov = vcovHC(Model3)) 
coeftest(Model4, vcov = vcovHC(Model4)) 
coeftest(Model5, vcov = vcovHC(Model5)) 
coeftest(Model6, vcov = vcovHC(Model6)) 
coeftest(Model7, vcov = vcovHC(Model7)) 
coeftest(Model8, vcov = vcovHC(Model8)) 


results = data.frame(cbind(area=c("cuneus", "frontal pole", "lateral occipital",
                                  "pars triangularis", "pars orbitalis", "rostral middle frontal", "superior temporal"),
                           #em=result[ c(8:41), 4]),
                           em = c(replicate(7, 1))),
                     #hemi=replicate(7, "left"),
                     stringsAsFactors=F)

df.log<-df #Plotting the raw values rather than the log transformed values
pdf(file="C:/Users/julie.wisch/Documents/SexDiffs/RegionalFigures_simplified.pdf") 
results %>% 
  ggseg(mapping=aes(fill=as.factor(em)), position="stacked", colour="black",size=.2,
        show.legend = F) + xlab("") + ylab("") + ggtitle("Regions with Sex-Linked \nTau Accumulation Differences") +
        scale_fill_manual(values =c("#fdae61"), na.value="grey85")


LinePlot<-function(X, Y, XLABEL, YLABEL){
  plot<-ggplot(df.log, aes(x = X, y = Y, shape = GENDER, fill = GENDER, colour = GENDER)) + 
  geom_point() + geom_smooth(method = "lm") + scale_colour_manual(values = c("#d7191c", "#2b83ba")) + 
  scale_shape_manual(values = c(13, 3))+ scale_fill_manual(values = c("#d7191c", "#2b83ba")) + labs(fill = "Sex", shape = "Sex", colour = "Sex") +
  theme(panel.background = element_rect(fill = NA),legend.position = "bottom") +
  xlab(XLABEL) + ylab(YLABEL)
  return(plot)}

LinePlot(df.log$PIB_fSUVR_rsf_TOT_CTX_SUPERTMP, df.log$TAU_fSUVR_rsf_TOT_CTX_SUPERTMP,
         "Amyloid Accumulation\nSuperior Temporal", "Tau Accumulation\nSuperior Temporal")


p1<-LinePlot(df.log$PIB_fSUVR_rsf_TOT_CTX_CUNEUS, df.log$TAU_fSUVR_rsf_TOT_CTX_CUNEUS,
         "Amyloid Accumulation\nCuneus", "Tau Accumulation\nCuneus")
p2<-LinePlot(df.log$PIB_fSUVR_rsf_TOT_CTX_FRNPOLE, df.log$TAU_fSUVR_rsf_TOT_CTX_FRNPOLE,
         "Amyloid Accumulation\nFrontal Pole", "Tau Accumulation\nFrontal Pole")
p3<-LinePlot(df.log$PIB_fSUVR_rsf_TOT_CTX_LATOCC, df.log$TAU_fSUVR_rsf_TOT_CTX_LATOCC,
         "Amyloid Accumulation\nLateral Occipital", "Tau Accumulation\nLateral Occipital")
p4<-LinePlot(df.log$PIB_fSUVR_rsf_TOT_CTX_PARSORBLS, df.log$TAU_fSUVR_rsf_TOT_CTX_PARSORBLS,
         "Amyloid Accumulation\nPars Orbitalis", "Tau Accumulation\nPars Orbitalis")
grid.arrange(p1, p2, p3, p4, nrow = 4)

dev.off()
#################################################################################################################
#################################################################################################################

#Women have more tau than men in the following regions:
# "cuneus", "frontal pole", "lateral occipital",
# "pars triangularis", "pars orbitalis" (inf frontal gyrus), "rostral middle frontal", "superior temporal"

#https://www.sciencedirect.com/science/article/pii/S0197458018303130
#Stephanie doesn't find a difference in ANY of these regions when comparing A- and A+ participants

#https://www.sciencedirect.com/science/article/pii/S2352872916300732
#In 30 - 49 year old individuals, "Despite the overall low-average SUVRs,
#occipital and frontal regions tended to have higher mean SUVR compared to medial temporal regions."

#Braak stages
# I/II	Hippocampus
# III	Parahippocampal gyrus; fusiform gyrus; lingual gyrus; amygdala
# IV	Inferior temporal cortex; middle temporal cortex; temporal pole; thalamus; posterior cingulate; insula
# V	Frontal cortex; parietal cortex; occipital cortex; superior temporal cortex; precuneus; caudate nucleus; putamen
# VI	Precentral gyrus; postcentral gyrus; paracentral gyrus; cuneus

#So women have more tau than men in the later Braak Stages....is that why they really fall off the cliff?
#We only have 6 people with CDR > 0 and a paired amyloid and tau scan, so I don't think we can test this hypothesis

Meds<-read.csv("C:/Users/julie.wisch/Downloads/DR_mednames_190118.csv")


Meds$hold<-apply(Meds, 1, function(r) any(r %in% c("estradiol", "estrogen", "premarin", "raloxifene", "estrace", "estrogens", 
                                                   "prempro", "estrogel", "fempatch")))
Meds$Estrogen<-ifelse(Meds$hold == TRUE, 1, 0)


Meds$hold<-apply(Meds, 1, function(r) any(r %in% c("citalopram", "sertraline", "escitalopram", "fluoxetine",
                                                   "lexapro", "zoloft", "celexa", "prozac", "paxil", "paroxetine", "sumatriptan")))
Meds$SSRI<-ifelse(Meds$hold == TRUE, 1, 0)

Meds$testdate<-as.Date(Meds$testdate, format = "%d-%b-%y")

#Ever Used estrogen
EverEstrogen<-unique(Meds[Meds$Estrogen == 1,1])
EverSSRI<-unique(Meds[Meds$SSRI == 1,1])


for(i in c(3:96, 98:191)){
  df.log[,i]<-log(df.log[,i])
}

Meds$EverEstrogen<-as.factor(ifelse(Meds$ID %in% EverEstrogen, 1, 0))
Meds$EverSSRI<-as.factor(ifelse(Meds$ID %in% EverSSRI, 1, 0))

df.meds<-df.log[,c("ID", "PET_Date_tau", "TAU_fSUVR_rsf_TOT_CTX_CUNEUS", "TAU_fSUVR_rsf_TOT_CTX_FRNPOLE" ,
                   "TAU_fSUVR_rsf_TOT_CTX_LATOCC", "TAU_fSUVR_rsf_TOT_CTX_PARSORBLS", "TAU_fSUVR_rsf_TOT_CTX_PARSTRNGLS",
                   "TAU_fSUVR_rsf_TOT_CTX_ROSMIDFRN", "TAU_fSUVR_rsf_TOT_CTX_SUPERTMP", "Tauopathy",
                   "PIB_fSUVR_rsf_TOT_CTX_CUNEUS", "PIB_fSUVR_rsf_TOT_CTX_FRNPOLE" ,
                   "PIB_fSUVR_rsf_TOT_CTX_LATOCC", "PIB_fSUVR_rsf_TOT_CTX_PARSORBLS", "PIB_fSUVR_rsf_TOT_CTX_PARSTRNGLS",
                   "PIB_fSUVR_rsf_TOT_CTX_ROSMIDFRN", "PIB_fSUVR_rsf_TOT_CTX_SUPERTMP", "PIB_fSUVR_rsf_TOT_CORTMEAN",
                   "TimeBetween","BIRTH","GENDER", "EDUC", "race2","CDR","TESTDATE","Age","apoe4",
                   "PIBpos", "Taupos","AT","VolRateofChange" )]
df.meds<-MatchbyNearestDate(df.meds, Meds[,c("ID", "testdate", "Estrogen", "SSRI", "EverEstrogen", "EverSSRI")], "ID", "TESTDATE", "testdate")


df.meds_f<-df.meds[df.meds$GENDER == "female",]


#When you control for SSRI, Statins and ACE use, there's no significnat relationship EXCEPT in overall tauopathy (where there was not a sex difference)
#If you have eer used an SSRI, you have a slightly higher estimated tauopathy, even after controlling for amyloid
#Notable because that's what that epi study said, except they didn't have imaging

Model.null<-lm(TAU_fSUVR_rsf_TOT_CTX_CUNEUS ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null)) #significant with estrogen. estrogen treated is lower
Model.null<-lm(TAU_fSUVR_rsf_TOT_CTX_FRNPOLE ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null))
Model.null<-lm(TAU_fSUVR_rsf_TOT_CTX_LATOCC ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null))
Model.null<-lm(TAU_fSUVR_rsf_TOT_CTX_PARSORBLS ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null)) #Trend with estrogen. estrogen treated is lower
Model.null<-lm(TAU_fSUVR_rsf_TOT_CTX_PARSTRNGLS ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null))
Model.null<-lm(TAU_fSUVR_rsf_TOT_CTX_ROSMIDFRN ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null))
Model.null<-lm(TAU_fSUVR_rsf_TOT_CTX_SUPERTMP ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null)) #Trend with estrogen. Estrogen treated is lower
Model.null<-lm(Tauopathy ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null)) 

df.meds_f$EverEstrogenRev<-ifelse(df.meds_f$EverEstrogen == 1, 0, 1)


library(corrplot)
df.meds.corr<-df.meds
colnames(df.meds.corr)[3:10]<-c("Cuneus", "FrnPole", "LatOcc", "ParsOrbls", "ParsTrngls", "RosMidFrn", "SuperTmp", "Tauopahty")
p.mat <- cor.mtest(df.meds.corr[,3:10])$p
corrplot(cor(df.meds.corr[,3:10]), method = "number", type = "upper",
         p.mat = p.mat, sig.level = 0.01, insig = "blank")


MakeBarPlots<-function(Component, ComponentTitle, SubTitle){
  theme_set(theme_bw())  
  pcaPlot<-as.data.frame(Component)
  pcaPlot<-cbind(rownames(pcaPlot), pcaPlot[,1])
  colnames(pcaPlot)<-c("region", "contribution")
  pcaPlot<-as.data.frame(pcaPlot)
  pcaPlot$contribution<-as.numeric(as.character(pcaPlot$contribution))
  pcaPlot <- pcaPlot[order(-pcaPlot$contribution), ]  # sort
  # Diverging Barcharts
  p<-ggplot(pcaPlot, aes(x=reorder(region, contribution), y = contribution, label=contribution)) + 
    geom_bar(stat='identity', width=.5)  +
    #coord_cartesian(ylim = c(-0.4, 0.4)) +
    scale_y_continuous()+ ylab("Contribution")+xlab("")+
    labs(subtitle=SubTitle, 
         title= ComponentTitle) + 
    coord_flip() + 
    theme(legend.position = "none")
  return(p)}

df_pca <- prcomp(df.meds_f[,c(3:9)], scale = TRUE, rotate = "varimax")
#Tau sig explains 52%
loadings <- df_pca$rotation
sdev <- df_pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
var.cos2 <- var.coord^2
# Compute contributions of each brain region to the component
#result is between 0 and 1...should sum to 1
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
row.names(var.contrib)<-c("Cuneus", "Frontal Pole", "Lateral Occipital", "Pars Orbitalis", "Pars Triangularis",
                          "Rostral Middle Frontal", "Superior Temporal")
pcap1<-MakeBarPlots(subset(var.contrib[,1], var.contrib[,1] > 0), "Tau Summary", "52% of variance")
df.meds_f$TauSig<-as.numeric(as.matrix(df.meds_f[,c(3:9)])%*%-df_pca$rotation[,1])


Model.null<-lm(TauSig ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 EverEstrogen, data = df.meds_f)
coeftest(Model.null, vcov = vcovHC(Model.null))

ggplot(df.meds_f, aes(x = EverEstrogen, y = TauSig)) + geom_point() + geom_boxplot()

abs(mean(df.meds_f[df.meds_f$EverEstrogen == 1, "TauSig"]) - mean(df.meds_f[df.meds_f$EverEstrogen == 0, "TauSig"]))/
  sqrt((sd(df.meds_f[df.meds_f$EverEstrogen == 1, "TauSig"]) + sd(df.meds_f[df.meds_f$EverEstrogen == 0, "TauSig"]))/2)

#Effect Size = 0.206 which is just baaaaarely into medium sized effect territory

abs(mean(df.meds_f[df.meds_f$EverSSRI == 1, "TauSig"]) - mean(df.meds_f[df.meds_f$EverSSRI == 0, "TauSig"]))/
  sqrt((sd(df.meds_f[df.meds_f$EverSSRI== 1, "TauSig"]) + sd(df.meds_f[df.meds_f$EverSSRI == 0, "TauSig"]))/2)
#Effect Size = 0.249 which is  medium sized effect territory


#Estrogen and/or SSRI's may play a role. We need to collect hormones, get a fuller picture of post-menopausal life for these women.


library(dabestr)  



df.meds$TauSig<-as.numeric(as.matrix(df.meds[,c(3:9)])%*%-df_pca$rotation[,1])


 ADI<-read.csv("C:/Users/julie.wisch/Documents/ADRC/ADI_all.csv")
 
df.meds<-merge(df.meds, ADI[,c("id", "ADI_NATRANK")], by.x = "ID", by.y = "id", all = FALSE)
med.fit <- glm(ADI_NATRANK ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER +  apoe4 + EDUC, data = df.meds, family = gaussian(link = "identity"))
out.fit <- glm(TauSig ~ ADI_NATRANK*GENDER + PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC,
               data = df.meds, family = gaussian(link = "identity"))
med.out <- mediate(med.fit, out.fit, treat = "GENDER", mediator = "ADI_NATRANK",
                   robustSE = TRUE, sims = 1000)
summary(med.out)

#Same situation. ADI does not mediate the difference even though there's an ADI x Sex Difference


#Now need to repeat with all available tau scans, regardless of if there's a proximate amyloid scan
tau.all<-read.csv("C:/Users/julie.wisch/Documents/ADRC/HASD_ACS_DR15_TAU.csv")
tau.all<-tau.all[,c("Map", "PET_Date", "TAU_fSUVR_rsf_TOT_CTX_CUNEUS", "TAU_fSUVR_rsf_TOT_CTX_FRNPOLE" ,
                    "TAU_fSUVR_rsf_TOT_CTX_LATOCC", "TAU_fSUVR_rsf_TOT_CTX_PARSORBLS", "TAU_fSUVR_rsf_TOT_CTX_PARSTRNGLS",
                    "TAU_fSUVR_rsf_TOT_CTX_ROSMIDFRN", "TAU_fSUVR_rsf_TOT_CTX_SUPERTMP", "Tauopathy")]
colnames(tau.all)[1]<-"ID"

for(i in c(3:10)){
  tau.all[,i]<-log(tau.all[,i])
}


df.demog<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR_demographics_20190122.csv")
df.demog$BIRTH<-format(as.Date(df.demog$BIRTH, "%d-%b-%y"), "19%y-%m-%d")
df.demog$BIRTH<-as.Date(df.demog$BIRTH, format = "%Y-%m-%d")

tau.all<-merge(tau.all, df.demog, by = "ID", all.y = FALSE)
tau.all$PET_Date<-as.Date(tau.all$PET_Date, format = "%m/%d/%Y")

tau.all$Age<-as.numeric(tau.all$PET_Date - tau.all$BIRTH)/365

tau.all$apoe4<-as.factor(ifelse(tau.all$apoe == 24 | tau.all$apoe > 33, 1, 0))

tau.all$TauSig<-as.numeric(as.matrix(tau.all[,c(3:9)])%*%-df_pca$rotation[,1])

tau.all<-MatchbyNearestDate(tau.all, Meds[,c("ID", "testdate", "Estrogen", "SSRI", "EverEstrogen", "EverSSRI")], 
                            "ID", "PET_Date", "testdate")

tau.all<-tau.all[complete.cases(tau.all[,c("EverEstrogen", "EverSSRI", "TauSig")]),]

tau.all<-aggregate(tau.all, list(tau.all$ID), FUN=head, 1)
table(tau.all$GENDER, tau.all$EverEstrogen) #28 females
table(tau.all$GENDER, tau.all$EverSSRI) #34 females, #8 males

tau.all.f<-tau.all[tau.all$GENDER == "female",]



Model.null<-lm(TauSig ~  Age +  apoe4 + EDUC + 
                 EverEstrogen, data = tau.all.f)
coeftest(Model.null, vcov = vcovHC(Model.null))
#There's a TREND to suggest that people who have ever taken estrogen have lower tau levels p = -0.115
abs(mean(tau.all.f[tau.all.f$EverEstrogen == 1, "TauSig"]) - mean(tau.all.f[tau.all.f$EverEstrogen == 0, "TauSig"]))/
  sqrt((sd(tau.all.f[tau.all.f$EverEstrogen == 1, "TauSig"]) + sd(tau.all.f[tau.all.f$EverEstrogen == 0, "TauSig"]))/2)




#Significant total effect and significant direct effect, but no mediating effect
#THat being said, it's really hard to tell. Only 8 men. Need more.

tau.all$GROUP<-as.factor(paste(tau.all$GENDER, tau.all$EverEstrogen, sep = "-"))
levels(tau.all$GROUP)<-c("Female, No Estrogen", "Female, Estrogen Use", "Male")


Model.null<-lm(TauSig ~ GROUP + Age + apoe4 + EDUC, tau.all)
coeftest(Model.null, vcov = vcovHC(Model.null))

Model.null<-aov(formula = TauSig ~ GROUP + Age + apoe4 + EDUC, tau.all)
TukeyHSD(Model.null, "GROUP")

colnames(tau.all)[21]<-"TauSig"

source("C:/Users/julie.wisch/Documents/DabestrMine/dabestr-master/R/plot.R")
source("C:/Users/julie.wisch/Documents/DabestrMine/dabestr-master/R/dabestr.R")
source("C:/Users/julie.wisch/Documents/DabestrMine/dabestr-master/R/flat_violin.R")
source("C:/Users/julie.wisch/Documents/DabestrMine/dabestr-master/R/main.R")
source("C:/Users/julie.wisch/Documents/DabestrMine/dabestr-master/R/plot_helpers.R")


multi.group <- 
  tau.all %>%
  dabest(GROUP, TauSig, 
         idx = list(c("Female, No Estrogen",  "Female, Estrogen Use", "Male")),
         paired = FALSE
  )

multi.group
plot.dabest(multi.group)



myVars <- c("EDUC", "apoe4", "race2", "Age", "PIBpos", "Taupos", "TimeBetween", "ADI_NATRANK")
catVars <- c("apoe4",  "race2", "PIBpos", "Taupos")
CreateTableOne(vars = myVars, data = tau.all, factorVars = catVars, strata = c("GENDER", "EverEstrogen"))
df.log<-merge(df.log, df.meds[,c("ID", "EverEstrogen")], by = "ID", all.x = TRUE, all.y = FALSE)
CreateTableOne(vars = myVars, data = df.log, factorVars = catVars, strata = c("GENDER", "EverEstrogen"))
df.log$GROUP<-as.factor(paste(df.log$GENDER, df.log$EverEstrogen, sep = "-"))
df.log<-df.log[!(df.log$GROUP == "female-NA"),]
df.log$race2<-droplevels(df.log$race2)
res.aov <- aov(apoe4 ~ GROUP, data = df.log[!(df.log$GROUP == "female-1"),])
summary(res.aov)

chisq.test(table(df.log$race2, df.log$GROUP)) 



levels(df.log$GROUP)<-c("Female, No Estrogen", "Female, Estrogen Use", "NA", "Male")

#Now need to repeat with all available tau scans, regardless of if there's a proximate amyloid scan
PIB.all<-read.csv("C:/Users/julie.wisch/Documents/ADRC/HASD_ACS_DR14_PIB.csv")
PIB.all<-PIB.all[,c("ID", "PET_Date", "PUP_fSUVR_rsf_TOT_CTX_CUNEUS", "PUP_fSUVR_rsf_TOT_CTX_FRNPOLE" ,
                    "PUP_fSUVR_rsf_TOT_CTX_LATOCC", "PUP_fSUVR_rsf_TOT_CTX_PARSORBLS", "PUP_fSUVR_rsf_TOT_CTX_PARSTRNGLS",
                    "PUP_fSUVR_rsf_TOT_CTX_ROSMIDFRN", "PUP_fSUVR_rsf_TOT_CTX_SUPERTMP", "PUP_fSUVR_rsf_TOT_CORTMEAN")]


for(i in c(3:10)){
  PIB.all[,i]<-log(PIB.all[,i])
}



PIB.all<-merge(PIB.all, df.demog, by = "ID", all.y = FALSE)
PIB.all$PET_Date<-as.Date(PIB.all$PET_Date, format = "%m/%d/%Y")

PIB.all$Age<-as.numeric(PIB.all$PET_Date - PIB.all$BIRTH)/365

PIB.all$apoe4<-as.factor(ifelse(PIB.all$apoe == 24 | PIB.all$apoe > 33, 1, 0))

df_pca <- prcomp(PIB.all[,c(3:9)], scale = TRUE, rotate = "varimax")
#PIB sig explains 81%
loadings <- df_pca$rotation
sdev <- df_pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
var.cos2 <- var.coord^2
# Compute contributions of each brain region to the component
#result is between 0 and 1...should sum to 1
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
row.names(var.contrib)<-c("Cuneus", "Frontal Pole", "Lateral Occipital", "Pars Orbitalis", "Pars Triangularis",
                          "Rostral Middle Frontal", "Superior Temporal")

pcap2<-MakeBarPlots(subset(var.contrib[,1], var.contrib[,1] > 0), "PIB Summary", "81% of variance")
grid.arrange(pcap1, pcap2, nrow = 1)
#########################################################
PIB.all$PIBSig<-as.numeric(as.matrix(PIB.all[,c(3:9)])%*%-df_pca$rotation[,1])

PIB.all<-MatchbyNearestDate(PIB.all, Meds[,c("ID", "testdate", "Estrogen", "SSRI", "EverEstrogen", "EverSSRI")], 
                            "ID", "PET_Date", "testdate")

PIB.all<-PIB.all[complete.cases(PIB.all[,c("EverEstrogen",  "PIBSig")]),]

PIB.all<-aggregate(PIB.all, list(PIB.all$ID), FUN=head, 1)
table(PIB.all$GENDER, PIB.all$EverEstrogen) #28 females

PIB.all$GROUP<-as.factor(paste(PIB.all$GENDER, PIB.all$EverEstrogen))
PIB.all<-PIB.all[!(PIB.all$GROUP == "male 1"),]
levels(PIB.all$GROUP)<-c("Female, No Estrogen", "Female, Estrogen Use", "Male", "Dropped")

multi.group <- 
  PIB.all %>%
  dabest(GROUP, PUP_fSUVR_rsf_TOT_CORTMEAN, 
         idx = list(c("Female, No Estrogen",  "Female, Estrogen Use", "Male")),
         paired = FALSE
  )

multi.group
plot.dabest(multi.group)



Model.null<-lm(TauSig ~  Age +  apoe4 + EDUC + 
                 GROUP, data = tau.all)
coeftest(Model.null, vcov = vcovHC(Model.null))
#Female, estrogen -0.102, p = 0.094
#Male, -0.461, p < 2e-16


Model.null<-lm(PIBSig ~  Age +  apoe4 + EDUC + 
                 GROUP, data = PIB.all)
coeftest(Model.null, vcov = vcovHC(Model.null))
#Female, Estrogen 0.137, p = 0.355
#Male 0.130, p = 0.0798

#Males have more amyloid than females, but even after that, females have more tau
#This strengthens our point

#Comparing everything with cognition

df.psych<-read.csv("C:/Users/julie.wisch/Documents/Aschenbrenner/ances_psych_092618.csv")
df.clin<-read.csv("C:/Users/julie.wisch/Documents/Aschenbrenner/ances_clinical_011819.csv")


#srtfree, digsym, memunits (or lmdelay) and MMSE. 
df.psych<-df.psych[,c("ID", "psy_date", "srtfree", "digsym", "MEMUNITS", "lmdelay")]
df.clin<-df.clin[,c("ID", "cdr", "testdate", "MMSE")]

df.psych$ID<-as.factor(df.psych$ID)
df.psych$psy_date<-as.Date(df.psych$psy_date, format = "%d-%b-%y")

df.clin$ID<-as.factor(df.clin$ID)
df.clin$testdate<-as.Date(df.clin$testdate, format = "%d-%b-%y")

df<-MatchbyNearestDate(df.psych, df.clin, "ID", "psy_date", "testdate")

df$MMSE <- as.numeric(ifelse(!is.na(df$MMSE) & df$MMSE > 30, "NA", df$MMSE)) #Some MMSE values were > 30.  Doesn't make sense.  dropped them.

PIB.all<-MatchbyNearestDate(PIB.all, df, "ID", "PET_Date", "psy_date")
tau.all<-MatchbyNearestDate(tau.all, df, "ID", "PET_Date", "psy_date")


PIB.all[,c(29:31, 35)]<- apply(PIB.all[,c(29:31, 35)], 2, scale) 
tau.all[,c(29:31, 35)]<- apply(tau.all[,c(29:31, 35)], 2, scale) 


for(i in 1:length(PIB.all$ID)){
  PIB.all$nacount[i]<-sum(is.na(PIB.all[i, c(29:31, 35)]))}
for(i in 1:length(tau.all$ID)){
  tau.all$nacount[i]<-sum(is.na(tau.all[i, c(29:31, 35)]))}


GetPACC<-function(x){
  x[is.na(x)] <- 0
  PACC<-(x[,"srtfree"]+x[,"digsym"]+x[,"MEMUNITS"]+x[,"MMSE"])/(4-x[,"nacount"])
  return(PACC)}


for(i in 1:length(PIB.all$ID)){
  PIB.all$PACC[i]<-GetPACC(PIB.all[i,])}
for(i in 1:length(tau.all$ID)){
  tau.all$PACC[i]<-GetPACC(tau.all[i,])}



cor.test(tau.all$Tauopathy, tau.all$PACC)
cor.test(tau.all$TauSig, tau.all$PACC)

cor.test(PIB.all$PUP_fSUVR_rsf_TOT_CORTMEAN, PIB.all$PACC)
cor.test(PIB.all$PIBSig, PIB.all$PACC)


Model.null<-lm(PACC ~  PIBSig + Age +  apoe4 + EDUC + 
                 GROUP + PIBSig:GROUP, data = PIB.all)
coeftest(Model.null, vcov = vcovHC(Model.null))
#More amyloid means you do worse, regardless of group
#Also males have a lower baseline

Model.null<-lm(PACC ~  PUP_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                 GROUP + PUP_fSUVR_rsf_TOT_CORTMEAN:GROUP, data = PIB.all)
coeftest(Model.null, vcov = vcovHC(Model.null))
#More overall amyloid means you do worse, regardless of group
#Also males have a lower baseline

Model.null<-lm(PACC ~  TauSig + Age +  apoe4 + EDUC + 
                 GROUP + TauSig:GROUP, data = tau.all)
coeftest(Model.null, vcov = vcovHC(Model.null)) 
#More tau in these regions mean you do worse, regardless of group
#Also males have a lower baseline

Model.null<-lm(PACC ~  Tauopathy + Age +  apoe4 + EDUC + 
                 GROUP + Tauopathy:GROUP, data = tau.all)
coeftest(Model.null, vcov = vcovHC(Model.null))
#More tauopathy mean you do worse, regardless of group
#No sex difference in this model


PIB.all$GROUP<-droplevels(PIB.all$GROUP)
tau.all$GROUP<-droplevels(tau.all$GROUP)

PACC_pib<-ggplot(PIB.all, aes(x = PUP_fSUVR_rsf_TOT_CORTMEAN, y = PACC, colour = GROUP, shape = GROUP, group = GROUP)) + geom_point() + 
  geom_smooth(method = "lm", fullrange = TRUE) + theme(legend.position = "none",panel.grid.major = element_blank(), 
                                                         panel.grid.minor = element_blank(),
                                                         panel.background = element_blank(), 
                                                         axis.line = element_line(colour = "black"))+ 
  scale_color_manual(values = c("#d7191c", "#fdae61", "#2b83ba")) +
  scale_shape_manual(values = c(13, 1, 3)) + xlab("Global Amyloid Burden")


PACC_Tau<-ggplot(tau.all, aes(x = TauSig, y = PACC, colour = GROUP, shape = GROUP, group = GROUP)) + geom_point() + 
  geom_smooth(method = "lm", fullrange = TRUE) + theme(legend.position = "none",panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), 
                                                       axis.line = element_line(colour = "black"))+ 
  scale_color_manual(values = c("#d7191c", "#fdae61", "#2b83ba")) +
  scale_shape_manual(values = c(13, 1, 3)) + xlab("Regional Tau Accumulation Signature")


grid.arrange(PACC_pib, PACC_Tau, nrow = 1)




##############################


CreateTableOne(vars = myVars, data = PIB.all, factorVars = catVars, strata = c("GENDER", "EverEstrogen"))
CreateTableOne(vars = myVars, data = tau.all, factorVars = catVars, strata = c("GENDER", "EverEstrogen"))

Demogs<-data.frame(rbind(PIB.all[,c("ID", "EDUC", "apoe4", "race2", "Age", "GENDER", "EverEstrogen")],
                         tau.all[,c("ID", "EDUC", "apoe4", "race2", "Age", "GENDER", "EverEstrogen")]))
Demogs<-Demogs[!duplicated(Demogs$ID),]
Demogs$apoe4<-as.factor(Demogs$apoe4)
Demogs$race2<-as.factor(Demogs$race2)
Demogs$Group<-as.factor(paste(Demogs$GENDER, Demogs$EverEstrogen, sep = ""))
CreateTableOne(vars = myVars, data = Demogs, factorVars = catVars, strata = "Group")


