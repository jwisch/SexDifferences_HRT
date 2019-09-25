library(plyr)
Hippo<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR15/HASD_ACS_DR15_3TMR.csv")
colnames(Hippo)[5]<-"ID"
Hippo<-Hippo[!(Hippo$Scanner == "CCIR Vida"),]

#Hippo<-read.csv("C:/Users/julie.wisch/Documents/ADRC/HASD_ACS_DR14_3TMR.csv")
#Hippo<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/HASD_ACS_DR14_3TMR.csv")
#Hippo<-Hippo[,c(4:5, 21:22, 343, 55:56)]
Hippo<-Hippo[,c("ID", "MR_Date", "MR_LV_CBLL_CORTEX", "MR_RV_CBLL_CORTEX", "MR_TOTV_INTRACRANIAL",
                "MR_LV_HIPPOCAMPUS", "MR_RV_HIPPOCAMPUS")]

df.clin<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR_clinical_20190122.csv")
#df.clin<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/DR_clinical_20190122.csv")
df.clin$TESTDATE<-as.Date(df.clin$TESTDATE, format = "%Y-%m-%d")
df.clin<-df.clin[,2:4]

df.demog<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR_demographics_20190122.csv")
#df.demog<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/DR_demographics_20190122.csv")
df.demog$BIRTH<-format(as.Date(df.demog$BIRTH, "%d-%b-%y"), "19%y-%m-%d")

df.demog<-merge(df.demog, df.clin, by = "ID", all.y = FALSE)

Hippo$MR_Date<-as.Date(Hippo$MR_Date, format = "%m/%d/%Y")
Hippo<-MatchbyNearestDate(Hippo, df.demog, "ID", "MR_Date", "TESTDATE") #have multiple timepoints for a TON of these guys

Hippo$apoe4<-as.factor(ifelse(Hippo$apoe == 24 | Hippo$apoe > 33, 1, 0))
Hippo<-Hippo[complete.cases(Hippo$apoe4),]
Hippo<-Hippo[Hippo$CDR == 0,]#only cognitively normals
Hippo<-Hippo[!is.na(Hippo$MR_TOTV_INTRACRANIAL),]

Hippo_Ncat<-Hippo
Scale<-function(COLUMN, COLNAME){
  MEAN<-mean(Hippo_Ncat$MR_TOTV_INTRACRANIAL)
  model<-lm(COLUMN ~ MR_TOTV_INTRACRANIAL, data = Hippo_Ncat)
  Hippo_Ncat[,paste("Normalized", COLNAME, sep = "")]<-(COLUMN - (coef(model)[2] * (Hippo_Ncat$MR_TOTV_INTRACRANIAL - MEAN)))
  return(Hippo_Ncat)}

Hippo_Ncat<-Scale(Hippo_Ncat$MR_LV_HIPPOCAMPUS, "LHippo")
Hippo_Ncat<-Scale(Hippo_Ncat$MR_RV_HIPPOCAMPUS, "RHippo")
Hippo_Ncat$Normalized<-(Hippo_Ncat$NormalizedLHippo + Hippo_Ncat$NormalizedRHippo)/2
Hippo_Ncat$Ndegen<-ifelse(Hippo_Ncat$Normalized < (mean(Hippo_Ncat$Normalized) - 1.5*sd(Hippo_Ncat$Normalized)), 1, 0)


ScaleandNormalize<-function(COLUMN, COLNAME){
  MEAN<-mean(Hippo$MR_TOTV_INTRACRANIAL)
  model<-lm(COLUMN ~ MR_TOTV_INTRACRANIAL, data = Hippo)
  Hippo[,paste("Normalized", COLNAME, sep = "")]<-scale(COLUMN - (coef(model)[2] * (Hippo$MR_TOTV_INTRACRANIAL - MEAN)))
  return(Hippo)}

for(i in c(3:4, 6:7)){
  Hippo<-ScaleandNormalize(Hippo[,i],names(Hippo[i]))
}

Hippo$Vol<-(Hippo$NormalizedMR_LV_HIPPOCAMPUS + Hippo$NormalizedMR_RV_HIPPOCAMPUS)/2
Hippo$CerebVol<-(Hippo$NormalizedMR_LV_CBLL_CORTEX + Hippo$NormalizedMR_RV_CBLL_CORTEX)/2

Hippo$Age<-as.numeric(difftime(Hippo$MR_Date, Hippo$BIRTH, units = "weeks")/52)

Hippo<-Hippo[with(Hippo, order(ID, MR_Date)),]

#Get PACC
df.psych<-read.csv("C:/Users/julie.wisch/Documents/Aschenbrenner/ances_psych_092618.csv")
df.clin<-read.csv("C:/Users/julie.wisch/Documents/Aschenbrenner/ances_clinical_011819.csv")


#srtfree, digsym, memunits (or lmdelay) and MMSE. 
df.psych<-df.psych[,c("ID", "psy_date", "srtfree", "digsym", "MEMUNITS", "lmdelay")]
df.clin<-df.clin[,c("ID", "cdr", "testdate", "MMSE")]

df.psych$ID<-as.factor(df.psych$ID)
df.psych$psy_date<-as.Date(df.psych$psy_date, format = "%d-%b-%y")

df.clin$ID<-as.factor(df.clin$ID)
df.clin$testdate<-as.Date(df.clin$testdate, format = "%d-%b-%y")

df.psych<-MatchbyNearestDate(df.psych, df.clin, "ID", "psy_date", "testdate")
df.psych$MMSE<-as.numeric(ifelse(!is.na(df.psych$MMSE) & df.psych$MMSE > 30, "NA", df.psych$MMSE)) #Some MMSE values were > 30.  Doesn't make sense.  dropped them.

df.psych<-df.psych[,-5]


for(i in 1:length(df.psych$ID)){df.psych$nacount[i]<-sum(is.na(df.psych[i, c(3:5, 8)]))}

df.psych<-df.psych[df.psych$nacount < 2,]

df.psych[,c(3:6, 8)]<- apply(df.psych[,c(3:6, 8)], 2, scale) 




GetPACC<-function(x){
  x[is.na(x)] <- 0
  PACC<-(x[,"srtfree"]+x[,"digsym"]+x[,"lmdelay"]+x[,"MMSE"])/(4-x[,"nacount"])
  return(PACC)}


for(i in 1:length(df.psych$ID)){df.psych$PACC[i]<-GetPACC(df.psych[i,])}
df.demog<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR_demographics_20190122.csv")

df.psych<-merge(df.demog, df.psych, by = "ID", all.x = FALSE, all.y = TRUE)
rm(df2, df.clin, df.demog, cat1, cat2, i, p1, p2, p3, p4, p5, p6, Model1, Model2, Model7, Model8, hold)
df.psych$BIRTH<-as.Date(df.psych$BIRTH, format = "%d-%b-%y")
df.psych$testdate<-as.Date(df.psych$testdate, format = "%Y-%m-%d")
df.psych$BIRTH<-df.psych$BIRTH - 100*365.25
df.psych$Age<-as.numeric(df.psych$testdate - df.psych$BIRTH)/365
df.psych<-df.psych[df.psych$Age < 110,]

df.psych$apoe4<-as.factor(ifelse(df.psych$apoe == 24 | df.psych$apoe > 33, 1, 0))


###CSF

CSF<-read.csv("C:/Users/julie.wisch/Documents/ADRC/OI_Schindler_2018_01_01.csv")
#CSF<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/OI_Schindler_2018_01_01.csv")

CSF<-CSF[CSF$CDR == 0,] #only keeping cognitively normals
CSF<-CSF[,c(1:4, 6:11, 21:26)]
CSF$LP_date<-as.Date(CSF$LP_date, format = "%d-%b-%y")
CSF<-CSF[complete.cases(CSF$GENDER),]
CSF$GENDER<-revalue(as.factor(CSF$GENDER), c("2"="female", "1"="male"))
CSF<-CSF[with(CSF, order(MAP_ID, LP_date)),]
CSF$PIBposbyCSF<-ifelse(CSF$E_ptau/CSF$E_ab42 <  0.0198,  0, 1)



CSF<-CSF[complete.cases(CSF$age_at_lp),]
colnames(CSF)[5]<-"apoe4"
CSF$logAB42<-log(CSF$E_ab42)
CSF$logptau<-log(CSF$E_ptau)
CSF<-CSF[complete.cases(CSF$apoe4),]
CSF.once<-ddply(CSF,.(MAP_ID), tail, 1)

CSF$apoe4<-as.factor(CSF$apoe4)


#NFL
NFL<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR15/OI_Schindler_FullNFL_2019_06_10.csv")
NFL<-NFL[,c("MAP_ID", "AGE_AT_LP", "GENDER", "EDUC", "APOE_10", "Race", "CDR",
            "LP_date","NFL")]
NFL$logNFL<-log(NFL$NFL)
colnames(NFL)[1]<-"ID"
NFL$LP_date<-as.Date(NFL$LP_date, format= "%m/%d/%Y")
NFL<-NFL[NFL$CDR==0,]

NFL<-NFL[with(NFL, order(ID, LP_date)),]


NFL.recent<-ddply(NFL,.(ID), tail, 1) #Keeping just the most recent NFL value

colnames(NFL)[2]<-"Age"
colnames(NFL)[5]<-"apoe4"
NFL<-NFL[complete.cases(NFL$logNFL) & complete.cases(NFL$Age) & complete.cases(NFL$apoe4),]
Hippo<-Hippo[,c("ID", "MR_Date", "BIRTH", "GENDER", "EDUC", "apoe4", "race2", "Vol", "CerebVol", "Age", "TESTDATE", "CDR")]

Hippo.hold<-Hippo
Hippo.all<-Hippo
Hippo<-ddply(Hippo.all,.(ID), tail, 1) #Keeping just the most recent volume for hippo

Hippo.head<-ddply(Hippo.all,.(ID), head, 1)

Hipporesult<-merge(Hippo.head[,c("ID", "Age", "Vol", "CerebVol")], Hippo[,c("ID", "Age", "Vol", "CerebVol")], by = "ID")
Hipporesult$VolRateofChange<-(Hipporesult$Vol.y - Hipporesult$Vol.x)/(Hipporesult$Age.y - Hipporesult$Age.x)
Hipporesult$CerebVolRateofChange<-(Hipporesult$CerebVol.y - Hipporesult$CerebVol.x)/(Hipporesult$Age.y - Hipporesult$Age.x)
Hipporesult<-Hipporesult[Hipporesult$VolRateofChange < 10 | Hipporesult$VolRateofChange > -10,]
Hipporesult<-Hipporesult[complete.cases(Hipporesult),]
lm_coefs<-Hipporesult[,c("ID", "VolRateofChange")]
Hippo.all<-merge(Hippo.all, lm_coefs[,c("ID", "VolRateofChange")], by = "ID", all.x = FALSE, all.y = FALSE)

lm_coefs<-Hipporesult[,c("ID", "CerebVolRateofChange")]

Hippo.all<-merge(Hippo.all, lm_coefs[,c("ID", "CerebVolRateofChange")], by = "ID", all.x = FALSE, all.y = FALSE)
Hippo.all<-Hippo.all[,c("ID", "GENDER", "EDUC", "apoe4", "race2", "Age", "VolRateofChange", "CerebVolRateofChange")]
Hippo.all<-data.frame(Hippo.all[!duplicated(Hippo.all$ID), ])
Hippo.all$VolRateofChange<-as.numeric(Hippo.all$VolRateofChange)
Hippo.all$CerebVolRateofChange<-as.numeric(Hippo.all$CerebVolRateofChange)
rm(df.clin,  lm_coefs)

df.pib<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR15/HASD_ACS_DR15_PIB.csv")
df.tau<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR15/HASD_ACS_DR15_TAU.csv")
# df.pib<-read.csv("C:/Users/julie.wisch/Documents/ADRC/HASD_ACS_DR14_PIB.csv")
# df.tau<-read.csv("C:/Users/julie.wisch/Documents/SexDiffs/HASD_ACS_DR14_TAU.csv")

# df.pib<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/HASD_ACS_DR14_PIB.csv")
# df.tau<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/HASD_ACS_DR14_TAU.csv")
df.tau<-df.tau[!(df.tau$PUP_QC_Status == "Failed" | df.tau$PUP_QC_Status == "Quarantined"),]
df.tau<-df.tau[,c(5:6, 806:831, 867:879, 915:924, 938:941, 955:973, 1009:1021, 1057:1064, 25)]
df.tau$PET_Date<-as.Date(df.tau$PET_Date, format = "%m/%d/%Y")
colnames(df.tau)[1]<-"ID"

df.pib<-df.pib[!(df.pib$PUP_QC_Status == "Failed" | df.pib$PUP_QC_Status == "Quarantined"),]
df.pib<-df.pib[,c(5:6, 800:825, 861:873, 909:918, 932:935, 949:967, 1003:1015, 1051:1059)]
df.pib$PET_Date<-as.Date(df.pib$PET_Date, format = "%m/%d/%Y")
colnames(df.pib)[1]<-"ID"

#####Keeping only PIB positive people
#df.pib<-df.pib[df.pib$PUP_fSUVR_rsf_TOT_CORTMEAN > 1.42,]
#################
df.tau<-df.tau[!is.na(df.tau$Tauopathy),]
df.pib<-df.pib[!is.na(df.pib$PIB_fSUVR_rsf_TOT_CORTMEAN),]

getLongitudinal<-function(df){
  DT <- data.table(df)
  counts <- DT[, .(rowCount = .N), by = ID]
  counts<-counts[!is.na(counts$ID),]
  counts<-counts[counts$rowCount > 1,]
  df<-df[df[,"ID"] %in% counts$ID,]
  return(df)}

df.tau<-df.tau[with(df.tau, order(ID, PET_Date)),]
df.tau.long<-getLongitudinal(df.tau)
df.tau<-ddply(df.tau,.(ID), tail, 1)

df.pib<-df.pib[with(df.pib, order(ID, PET_Date)),]
df.pib.long<-getLongitudinal(df.pib)
df.pib<-ddply(df.pib,.(ID), tail, 1)

df.demog<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR_demographics_20190122.csv")
#df.demog<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/DR_demographics_20190122.csv")
df.demog$BIRTH<-format(as.Date(df.demog$BIRTH, "%d-%b-%y"), "19%y-%m-%d")

df.clin<-read.csv("C:/Users/julie.wisch/Documents/ADRC/DR_clinical_20190122.csv")
#df.clin<-read.csv("W:/ances/julie/Data/ADRC/Aggregated/DR_clinical_20190122.csv")
df.clin$TESTDATE<-as.Date(df.clin$TESTDATE, format = "%Y-%m-%d")

df.clin<-df.clin[,2:4]

df<-merge(df.demog, df.clin, by = "ID", all.y = FALSE)
df$BIRTH<-as.Date(df$BIRTH, format = "%Y-%m-%d")

colnames(df.pib)[2]<-"PET_Date_PIB"
colnames(df.tau)[2]<-"PET_Date_tau"
df.PET<-MatchbyNearestDate(df.tau, df.pib, "ID", "PET_Date_tau", "PET_Date_PIB")
df.PET$TimeBetween<-as.numeric((df.PET$PET_Date_PIB - df.PET$PET_Date_tau)/365)

df<-MatchbyNearestDate(df.PET, df, "ID", "PET_Date_tau", "TESTDATE")
df<-df[df$TimeBetween > -4,] #Making sure scans occurred within 4 years of each other




df$Age<-as.numeric(difftime(df$PET_Date_tau, df$BIRTH, units = "weeks")/52) #Age at tau scan
#df<-df[df$Age > 64.99,] #Restricting to only "older adults"

df<-df[df$CDR == 0,] #only cognitively normal

df$apoe4<-ifelse(df$apoe == 24 | df$apoe > 33, 1, 0)

df$PIBpos<-ifelse(df$PIB_fSUVR_rsf_TOT_CORTMEAN > 1.419999, 1, 0)

#Need to add back in after hearing back from Aylin

df$Taupos<-ifelse(df$Tauopathy > 1.22, 1, 0)
df$AT<-paste(df$PIBpos, df$Taupos, sep = "-")
#df<-df[!(df$AT == "0-1"),] #Dropping SNAP people

df.pib.long<-merge(df.demog[,c("ID", "BIRTH", "GENDER", "EDUC", "race2", "apoe")], df.pib.long, by = "ID", all.x = FALSE, all.y = TRUE)
df.pib.long$apoe4<-as.factor(ifelse(df.pib.long$apoe == 24 | df.pib.long$apoe > 33, 1, 0))
df.tau.long<-merge(df.demog[,c("ID", "BIRTH", "GENDER", "EDUC", "race2", "apoe")], df.tau.long, by = "ID", all.x = FALSE, all.y = TRUE)
df.tau.long$apoe4<-as.factor(ifelse(df.tau.long$apoe == 24 | df.tau.long$apoe > 33, 1, 0))
