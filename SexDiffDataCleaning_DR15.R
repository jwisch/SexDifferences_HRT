#Reading in PET Scans
FILEPATH<-"C:/Users/julie.wisch/Documents/ADRC/DR15/"

df.pib<-read.csv(paste(FILEPATH, "HASD_ACS_DR15_PIB.csv", sep = ""))
df.tau<-read.csv(paste(FILEPATH, "HASD_ACS_DR15_TAU.csv", sep = ""))
df.demog<-read.csv(paste(FILEPATH, "DR_demographics_20190122.csv", sep = ""))
df.clin<-read.csv(paste(FILEPATH, "DR_clinical_20190122.csv", sep = ""))



df.tau<-df.tau[!(df.tau$PUP_QC_Status == "Failed" | df.tau$PUP_QC_Status == "Quarantined"),]
df.tau<-df.tau[,c(5:6, 806:831, 867:879, 915:924, 938:941, 955:973, 1009:1021, 1057:1064, 25)]
df.tau$PET_Date<-as.Date(df.tau$PET_Date, format = "%m/%d/%Y")
colnames(df.tau)[1]<-"ID"

df.pib<-df.pib[!(df.pib$PUP_QC_Status == "Failed" | df.pib$PUP_QC_Status == "Quarantined"),]
df.pib<-df.pib[,c(5:6, 800:825, 861:873, 909:918, 932:935, 949:967, 1003:1015, 1051:1059)]
df.pib$PET_Date<-as.Date(df.pib$PET_Date, format = "%m/%d/%Y")
colnames(df.pib)[1]<-"ID"

#Dropping incomplete visits
df.tau<-df.tau[!is.na(df.tau$Tauopathy),]
df.pib<-df.pib[!is.na(df.pib$PIB_fSUVR_rsf_TOT_CORTMEAN),]



df.tau<-df.tau[with(df.tau, order(ID, PET_Date)),]


df.pib<-df.pib[with(df.pib, order(ID, PET_Date)),]


df.demog$BIRTH<-format(as.Date(df.demog$BIRTH, "%d-%b-%y"), "19%y-%m-%d")

df.clin$TESTDATE<-as.Date(df.clin$TESTDATE, format = "%Y-%m-%d")

df.clin<-df.clin[,2:4]

df<-merge(df.demog, df.clin, by = "ID", all.y = FALSE)
df$BIRTH<-as.Date(df$BIRTH, format = "%Y-%m-%d")

colnames(df.pib)[2]<-"PET_Date_PIB"
colnames(df.tau)[2]<-"PET_Date_tau"
df.PET<-MatchbyNearestDate(df.tau, df.pib, "ID", "PET_Date_tau", "PET_Date_PIB")
df.PET$TimeBetween<-as.numeric((df.PET$PET_Date_PIB - df.PET$PET_Date_tau)/365)

df<-MatchbyNearestDate(df.PET, df, "ID", "PET_Date_tau", "TESTDATE")
df<-df[order( df[,"ID"], df[,"TimeBetween"] ),]


df<-df[df$TimeBetween > -4,] #Making sure scans occurred within 4 years of each other



df$Age<-as.numeric(difftime(df$PET_Date_tau, df$BIRTH, units = "weeks")/52) #Age at tau scan

df<-df[df$CDR == 0,] #only cognitively normal

#Turning APOE into a dichotomous variable
df$apoe4<-ifelse(df$apoe == 24 | df$apoe > 33, 1, 0)

#Classifying people by amyloid and tau status
df$PIBpos<-ifelse(df$PIB_fSUVR_rsf_TOT_CORTMEAN > 1.419999, 1, 0)
df$Taupos<-ifelse(df$Tauopathy > 1.22, 1, 0)
