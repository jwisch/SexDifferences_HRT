
#Loading in functions used as well as cleaned dataframe for analysis
FILEPATH_OTHERCODE<-"C:/Users/julie.wisch/Documents/Transition/DraftedMS_SexDiffs/Code/"
source(paste(FILEPATH_OTHERCODE, "CommonFuncs.r", sep = ""))
source(paste(FILEPATH_OTHERCODE, "SexDiffFuncs.r", sep = ""))
source(paste(FILEPATH_OTHERCODE, "SexDiffDataCleaning_DR15_try2.r", sep = ""))

######################################################################################
#Loading in files not loaded in the datacleaning step
FILEPATH_DATA<-"C:/Users/julie.wisch/Documents/Transition/DraftedMS_SexDiffs/Data/"
Meds<-read.csv(paste(FILEPATH_DATA, "DR_mednames_190118.csv", sep = ""))
df.psych<-read.csv(paste(FILEPATH_DATA, "ances_psych_092618.csv", sep = ""))
df.clin<-read.csv(paste(FILEPATH_DATA, "ances_clinical_011819.csv", sep = ""))

########################################################################################
#Loading libraries in

library(eeptools)
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
library(RColorBrewer)


#################################################################################
#Testing for sex differences
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
df.log$apoe4<-as.factor(df.log$apoe4)
result<-data.frame("Region" = character(), "Amyloid" = numeric(), "Sex" = numeric(), 
                   "APOE" = numeric(), "Interaction" = numeric(), "Sex_noAmyloid" = numeric(),
                   "APOE_noAmyloid" = numeric(), "Interaction_noAmyloid" = numeric(),
                   "Tau_Cortmean" = numeric(),
                   "Sex_Cortmean" = numeric(),
                   "APOE_Cortmean" = numeric(), "Interaction_Cortmean" = numeric(),
                   stringsAsFactors = FALSE)
k<-0
#Looping through and testing for a difference in tau burden for a given level of amyloid for every region individually
#Using the function coeftest from the sandwich package to get heteroskedacity consistent estimators
#Lots of these models violate the heteroskedacity argument. this is one way around that
for(i in c(15:48)){
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
result$BonfCorrected<-p.adjust(result$Sex_Cortmean, method = "bonferroni", n = length(result$Sex_Cortmean))
result$BonfCorrected_NoA<-p.adjust(result$Sex_noAmyloid, method = "bonferroni", n = length(result$Sex_Cortmean))

#Significant after bonferroni correction
result[result$BonfCorrected < 0.05,c(1, 3, 6, 10, 13)] #After you control for amyloid
result[result$BonfCorrected_NoA < 0.05,c(1, 3, 6, 10, 14)] #If you don't control for amyloid
#All regions differ by sex regardless of whether or not you control for amyloid
#The one exception is the middle temporal which is only significant if you control for amyloid burden

#Visualizing the regions with significant sex differences
results = data.frame(cbind(area=c("cuneus", "frontal pole", "inferior parietal",
                                  "lateral occipital", "middle temporal",
                                  "pars triangularis", "pars orbitalis", "rostral middle frontal", "superior temporal",
                                  "temporal pole"),
                           #em=result[ c(8:41), 4]),
                           em = c(replicate(10, 1))),
                     #hemi=replicate(7, "left"),
                     stringsAsFactors=F)
brewer.pal(2, "Set1")

# to display that palette:
display.brewer.pal(2, "Set1")

results %>% 
  ggseg(mapping=aes(fill=as.factor(em)), position="stacked", colour="black",size=.2,
        show.legend = F) + xlab("") + ylab("") + ggtitle("Regions with Sex-Linked \nTau Accumulation Differences") +
  scale_fill_manual(values =c("#377EB8"), na.value="grey85")
########################################################################################

########################################################################################
#Now adding in HT use history

Meds$hold<-apply(Meds, 1, function(r) any(r %in% c("estradiol", "estrogen", "premarin", "raloxifene", "estrace", "estrogens", 
                                                   "prempro", "estrogel", "fempatch")))
Meds$Estrogen<-ifelse(Meds$hold == TRUE, 1, 0)
Meds$testdate<-as.Date(Meds$testdate, format = "%d-%b-%y")

#Ever Used estrogen
EverEstrogen<-unique(Meds[Meds$Estrogen == 1,1])

Meds$EverEstrogen<-as.factor(ifelse(Meds$ID %in% EverEstrogen, 1, 0))


df.meds<-df.log[,c("ID", "PET_Date_tau", "TAU_fSUVR_rsf_TOT_CTX_CUNEUS", "TAU_fSUVR_rsf_TOT_CTX_FRNPOLE" ,
                   "TAU_fSUVR_rsf_TOT_CTX_INFERPRTL" ,
                   "TAU_fSUVR_rsf_TOT_CTX_LATOCC", "TAU_fSUVR_rsf_TOT_CTX_MIDTMP",
                   "TAU_fSUVR_rsf_TOT_CTX_PARSORBLS", "TAU_fSUVR_rsf_TOT_CTX_PARSTRNGLS",
                   "TAU_fSUVR_rsf_TOT_CTX_ROSMIDFRN", "TAU_fSUVR_rsf_TOT_CTX_SUPERTMP", 
                   "TAU_fSUVR_rsf_TOT_CTX_TMPPOLE", "Tauopathy",
                   "PIB_fSUVR_rsf_TOT_CTX_CUNEUS", "PIB_fSUVR_rsf_TOT_CTX_FRNPOLE" ,
                   "PIB_fSUVR_rsf_TOT_CTX_INFERPRTL",
                   "PIB_fSUVR_rsf_TOT_CTX_LATOCC", "PIB_fSUVR_rsf_TOT_CTX_MIDTMP",
                   "PIB_fSUVR_rsf_TOT_CTX_PARSORBLS", "PIB_fSUVR_rsf_TOT_CTX_PARSTRNGLS",
                   "PIB_fSUVR_rsf_TOT_CTX_ROSMIDFRN", "PIB_fSUVR_rsf_TOT_CTX_SUPERTMP", 
                   "PIB_fSUVR_rsf_TOT_CTX_TMPPOLE", "PIB_fSUVR_rsf_TOT_CORTMEAN",
                   "TimeBetween","BIRTH","GENDER", "EDUC", "race2","CDR","TESTDATE","Age","apoe4",
                   "PIBpos", "Taupos")]
df.meds<-MatchbyNearestDate(df.meds, Meds[,c("ID", "testdate", "Estrogen", "EverEstrogen")], "ID", "TESTDATE", "testdate")
df.meds$Group<-paste(df.meds$GENDER, df.meds$EverEstrogen, sep = "-")
###########################################################

############################################################
#Testing for differences by HT use history

k<-0
plist<-data.frame("Region" = character(), "P" = numeric(), stringsAsFactors = FALSE)
for(i in 3:12){
  k<-k+1
  Model.null<-lm(df.meds[df.meds$GENDER == "female", i] ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age +  apoe4 + EDUC + 
                   EverEstrogen, data = df.meds[df.meds$GENDER == "female",])
  print(names(df.meds[i]))
  hold<-coeftest(Model.null, vcov = vcovHC(Model.null))
  plist[k, "Region"]<-as.character(names(df.meds[i]))
  plist[k, "P"]<-hold[6,4]}


plist$BonfCorrected<-p.adjust(plist$P, method = "bonferroni", n = length(plist$P))
plist
library(margins)
#I only tested for HT differences in the regions where we had already identified a sex difference
Model1<-lm(TAU_fSUVR_rsf_TOT_CTX_CUNEUS ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model2<-lm(TAU_fSUVR_rsf_TOT_CTX_FRNPOLE ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model3<-lm(TAU_fSUVR_rsf_TOT_CTX_INFERPRTL ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model4<-lm(TAU_fSUVR_rsf_TOT_CTX_LATOCC ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model5<-lm(TAU_fSUVR_rsf_TOT_CTX_MIDTMP ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model6<-lm(TAU_fSUVR_rsf_TOT_CTX_PARSORBLS ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model7<-lm(TAU_fSUVR_rsf_TOT_CTX_PARSTRNGLS ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model8<-lm(TAU_fSUVR_rsf_TOT_CTX_ROSMIDFRN ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model9<-lm(TAU_fSUVR_rsf_TOT_CTX_SUPERTMP ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model10<-lm(TAU_fSUVR_rsf_TOT_CTX_TMPPOLE ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model11<-lm(Tauopathy ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + apoe4 + EDUC + Group, data = df.meds)
Model12<-lm( PIB_fSUVR_rsf_TOT_CORTMEAN ~ Age + apoe4 + EDUC + Group, data = df.meds)

#Visualizing the results
library(jtools)
ForPlot<-data.frame("Region" = rep(NA, 10), "Coefficient" = rep(NA, 10), "StdErr" = rep(NA, 10))
ForPlot$Region<-c("Cuneus", "Frontal Pole", "Inferior Parietal", "Lateral Occipital", "Middle Temporal",
                  "Pars Orbitalis", "Pars Triangularis",
                  "Rostral Middle Frontal", "Superior Temporal", "Temporal Pole")
ForPlot$Coefficient[1]<-as.numeric(coef(Model1)[4])
ForPlot$StdErr[1]<-coef(summary(Model1))[,2][4]
ForPlot$Coefficient[2]<-as.numeric(coef(Model2)[4])
ForPlot$StdErr[2]<-coef(summary(Model2))[,2][4]
ForPlot$Coefficient[3]<-as.numeric(coef(Model3)[4])
ForPlot$StdErr[3]<-coef(summary(Model3))[,2][4]
ForPlot$Coefficient[4]<-as.numeric(coef(Model4)[4])
ForPlot$StdErr[4]<-coef(summary(Model4))[,2][4]
ForPlot$Coefficient[5]<-as.numeric(coef(Model5)[4])
ForPlot$StdErr[5]<-coef(summary(Model5))[,2][4]
ForPlot$Coefficient[6]<-as.numeric(coef(Model6)[4])
ForPlot$StdErr[6]<-coef(summary(Model6))[,2][4]
ForPlot$Coefficient[7]<-as.numeric(coef(Model7)[4])
ForPlot$StdErr[7]<-coef(summary(Model7))[,2][4]
ForPlot$Coefficient[8]<-as.numeric(coef(Model8)[4])
ForPlot$StdErr[8]<-coef(summary(Model8))[,2][4]
ForPlot$Coefficient[9]<-as.numeric(coef(Model9)[4])
ForPlot$StdErr[9]<-coef(summary(Model9))[,2][4]
ForPlot$Coefficient[10]<-as.numeric(coef(Model10)[4])
ForPlot$StdErr[10]<-coef(summary(Model10))[,2][4]

MarginalEffects<-rbind(data.frame(summary(margins(Model1, variables = "Group"))), c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model2, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model3, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model4, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model5, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model6, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model7, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model8, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model9, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA),
data.frame(summary(margins(Model10, variables = "Group"))),c(NA, NA, NA, NA, NA, NA, NA))
MarginalEffects$Region<-c(rep("Cuneus", 3), rep("Frontal Pole", 3), rep("Inferior Parietal", 3),
                          rep("Lateral Occipital", 3),
                          rep("Middle Temporal", 3),
                        rep("Pars Orbitalis", 3), rep("Pars Triangularis", 3), rep("Rostral Middle Frontal", 3),
                        rep("Superior Temporal", 3), rep("Temporal Pole", 3))
MarginalEffects$RowCount<-as.factor(seq(1:length(MarginalEffects$Region)))
MarginalEffects$factor<-as.factor(MarginalEffects$factor)
MarginalEffects$factor<-droplevels(MarginalEffects$factor)
levels(MarginalEffects$factor)<-c("Female, HT User", "Male")
#Plotting the estimated marginal effects of sex and HT usage on tau burden in the different regions
ggplot(MarginalEffects, aes(x = RowCount, y = AME, group = factor, color = factor)) + geom_point() + 
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  geom_errorbar(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom", axis.ticks = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey37") +
  #manually writing out the axes labels
  scale_x_discrete(labels=c("1" = "", "2" = "Cuneus", "3" = "", "4" = "Pole", "5" = "Frontal", "6" = "",
                            "7" = "Parietal", "8" = "Inferior", "9" = "",
                            "10" = "Occipital", "11" = "Lateral", "12" = "", 
                            "13" = "Temporal", "14" = "Middle", "15" = "",
                            "16" = "Orbitalis", "17" = "Pars", "18" = "",
                            "19" = "Triangularis", "20" = "Pars", "21" = "", 
                            "22" = "Frontal", "23" = "Middle", "24" = "Rostral",
                            "25" = "Temporal", "26" = "Superior", "27" = "", "28" = "Pole", "29" = "Temporal", "30" = "")) +
  xlab("") + ylab("Estimated Marginal Effect\nof Sex and HT History") + 
  labs(group = "", color = "") +
  coord_flip() + theme(axis.ticks = element_blank())
  


library(corrplot)
df.meds.corr<-df.meds
colnames(df.meds.corr)[3:13]<-c("Cuneus", "Frontal Pole", "Inferior Parietal", "Lateral Occipital", 
                                "Middle Temporal", "Pars Orbitalis", "Pars Triangularis", 
                                "Rostral Middle Frontal", 
                                "Superior Temporal", "Temporal Pole", "Tauopathy")
p.mat <- cor.mtest(df.meds.corr[,3:13])$p
corrplot(cor(df.meds.corr[,3:13]), method = "number", type = "upper",
         p.mat = p.mat, sig.level = 0.05, insig = "blank", tl.col="black", tl.cex = 0.75)

################################################################################
#Performing PCA to get a single summary tau value rather than the 10 regions to increase power for subsequent analysis
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
    geom_bar(stat='identity', width=.5, fill = "#377EB8")  +
    #coord_cartesian(ylim = c(-0.4, 0.4)) +
    scale_y_continuous()+ ylab("Contribution")+xlab("")+
    labs(title= ComponentTitle) + 
    coord_flip() + 
    theme(legend.position = "none")
  return(p)}

df_pca <- prcomp(df.meds[,c(3:12)], scale = TRUE, rotate = "varimax")
#Tau sig explains 55%
loadings <- df_pca$rotation
sdev <- df_pca$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
var.cos2 <- var.coord^2
# Compute contributions of each brain region to the component
#result is between 0 and 1...should sum to 1
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
row.names(var.contrib)<-c("Cuneus", "Frontal Pole", "Inferior Parietal", "Lateral Occipital", 
                          "Middle Temporal", "Pars Orbitalis", "Pars Triangularis", 
                          "Rostral Middle Frontal", 
                          "Superior Temporal", "Temporal Pole")
MakeBarPlots(subset(var.contrib[,1], var.contrib[,1] > 0), "Tau Summary Composition", "")
df.meds$TauSig<-as.numeric(as.matrix(df.meds[,c(3:12)])%*%-df_pca$rotation[,1])








#Testing for differences by sex/history of HT in amyloid accumulation
model1<-lm(PIB_fSUVR_rsf_TOT_CORTMEAN ~ Age + GENDER + EDUC + apoe4, data = df.meds)
model2<-lm(PIB_fSUVR_rsf_TOT_CORTMEAN ~ Age + Group + EDUC + apoe4, data = df.meds)
coeftest(model1, vcov = vcovHC(model1)) #No sex difference in amyloid if I just keep this small subset of everyone
coeftest(model2, vcov = vcovHC(model2)) #No sex difference in amyloid if I just keep this small subset of everyone

p1<-plot_summs(model1, model2, scale = TRUE, plot.distributions = FALSE, 
               coefs = c("Participant Age" = "Age",
                         "Years Education" = "EDUC",
                         "APOE4 Status" = "apoe4",
                         "Males compared to\nall females" = "GENDER",
                         "Males compared to\nfemales without HT" = "Groupmale-0", 
                         "Females with HT compared\nto females without HT" = "Groupfemale-1"), colors = "Set1") + 
  theme(legend.position = "none", axis.text.y = element_text(size = 8)) + ggtitle("Cortical Amyloid Burden")


#Testing for differences by sex/history of HT in tau accumulation
model1<-lm(TauSig ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + GENDER + EDUC + apoe4, data = df.meds)
model2<-lm(TauSig ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + Group + EDUC + apoe4, data = df.meds)
coeftest(model1, vcov = vcovHC(model1)) #significant sex difference
coeftest(model2, vcov = vcovHC(model2)) #trend for tau difference

#Getting Effect Sizes -> Cohen's f2
model3<-lm(TauSig ~ PIB_fSUVR_rsf_TOT_CORTMEAN + Age + Group + EDUC + apoe4, data = df.meds[df.meds$GENDER == "female",])
## store ANOVA from model
msum <- anova(model3)
## divide all Sums of Squares by the sum() of the Sums of Squares
msum[["Sum Sq"]]/sum(msum[["Sum Sq"]]) #This gives Cohen's f2
#Cohen's f2 > 0.02 is small, f2 ??? 0.15 is medium, and f2 ??? 0.35 is large



p2<-plot_summs(model1, model2, scale = TRUE, plot.distributions = FALSE, 
               coefs = c("Participant Age" = "Age",
                         "Years Education" = "EDUC",
                         "APOE4 Status" = "apoe4",
                         "Males compared to\nall females" = "GENDER",
                         "Males compared to\nfemales without HT" = "Groupmale-0", 
                         "Females with HT compared\nto females without HT" = "Groupfemale-1"), colors = "Set1") + 
  theme(legend.position = "none", axis.text.y = element_text(size = 8)) + ggtitle("Summary Tau Burden")

grid.arrange(p1, p2, nrow = 2)


########################################################################################3
#Now adding in the PACC



#srtfree, digsym, memunits (or lmdelay) and MMSE. 
df.psych<-df.psych[,c("ID", "psy_date", "srtfree", "digsym", "MEMUNITS", "lmdelay")]
df.clin<-df.clin[,c("ID", "cdr", "testdate", "MMSE")]

df.psych$ID<-as.factor(df.psych$ID)
df.psych$psy_date<-as.Date(df.psych$psy_date, format = "%d-%b-%y")

df.clin$ID<-as.factor(df.clin$ID)
df.clin$testdate<-as.Date(df.clin$testdate, format = "%d-%b-%y")

df<-MatchbyNearestDate(df.psych, df.clin, "ID", "psy_date", "testdate")

df$MMSE <- as.numeric(ifelse(!is.na(df$MMSE) & df$MMSE > 30, "NA", df$MMSE)) #Some MMSE values were > 30.  Doesn't make sense.  dropped them.

df.meds<-MatchbyNearestDate(df.meds, df, "ID", "PET_Date_tau", "psy_date")

for(i in 1:length(df.meds$ID)){
  df.meds$nacount[i]<-sum(is.na(df.meds[i, c(42:44, 48)]))}


GetPACC<-function(x){
  x[is.na(x)] <- 0
  PACC<-(x[,"srtfree"]+x[,"digsym"]+x[,"MEMUNITS"]+x[,"MMSE"])/(4-x[,"nacount"])
  return(PACC)}


for(i in 1:length(df.meds$ID)){
  df.meds$PACC[i]<-GetPACC(df.meds[i,])}

#Centering amyloid and tau burden numbers so that interactions are easier to interpret
df.meds$PIB_centered<-df.meds$PIB_fSUVR_rsf_TOT_CORTMEAN - mean(df.meds$PIB_fSUVR_rsf_TOT_CORTMEAN, na.rm = TRUE)
df.meds$Tau_centered<-df.meds$TauSig - mean(df.meds$TauSig, na.rm = TRUE)


model10<-lm(PACC ~  PIB_centered + Age +  apoe4 + EDUC + 
                 Group + PIB_centered:Group, data = df.meds)
coeftest(model10, vcov = vcovHC(model10))

#Getting Effect Sizes
## store ANOVA from model
msum <- anova(model10)
## divide all Sums of Squares by the sum() of the Sums of Squares
msum[["Sum Sq"]]/sum(msum[["Sum Sq"]]) #This gives Cohen's f2
#Cohen's f2 > 0.02 is small, f2 ??? 0.15 is medium, and f2 ??? 0.35 is large



model11<-lm(PACC ~  Tau_centered + Age +  apoe4 + EDUC + 
                 Group + Tau_centered:Group, data = df.meds)
coeftest(model11, vcov = vcovHC(model11))

msum <- anova(model11)
## divide all Sums of Squares by the sum() of the Sums of Squares
msum[["Sum Sq"]]/sum(msum[["Sum Sq"]])


coeftest(model10, vcov = vcovHC(model10)) #significant sex difference
coeftest(model11, vcov = vcovHC(model11)) #trend for tau difference


#Just pulling all the individual amyloid/tau/group/sex effect sizes. Edit these ad hoc to get what you need.
model10a<-lm(PACC ~  PIB_centered + Age +  apoe4 + EDUC + 
              Group + PIB_centered:Group, data = df.meds)
model10b<-lm(PACC ~   Age +  apoe4 + EDUC + 
              Group, data = df.meds)
(summary(model10a)$r.squared - summary(model10b)$r.squared)/(1 - summary(model10a)$r.squared)

model11a<-lm(PACC ~  PIB_centered + Age +  apoe4 + EDUC + 
               Group + PIB_centered:Group, data = df.meds[df.meds$GENDER == "female",])
model11b<-lm(PACC ~   PIB_centered + Age +  apoe4 + EDUC + Group, data = df.meds[df.meds$GENDER == "female",])
(summary(model10a)$r.squared - summary(model10b)$r.squared)/(1 - summary(model10a)$r.squared)


model11a<-lm(PACC ~  Tau_centered + Age +  apoe4 + EDUC + 
               Group + Tau_centered:Group, data = df.meds)
model11b<-lm(PACC ~   Age +  apoe4 + EDUC + 
               Group, data = df.meds)
(summary(model11a)$r.squared - summary(model11b)$r.squared)/(1 - summary(model11a)$r.squared)

model11a<-lm(PACC ~  Tau_centered + Age +  apoe4 + EDUC + 
               Group + Tau_centered:Group, data = df.meds[df.meds$GENDER == "female",])
model11b<-lm(PACC ~   Tau_centered + Age +  apoe4 + EDUC, data = df.meds[df.meds$GENDER == "female",])
(summary(model11a)$r.squared - summary(model11b)$r.squared)/(1 - summary(model11a)$r.squared)




p1<-plot_summs(model10, scale = TRUE, plot.distributions = FALSE, 
               coefs = c("Tau_centered" = "Centered Tau Signature",
                          "Participant Age" = "Age",
                         "Years Education" = "EDUC",
                         "APOE4 Status" = "apoe4",
                         "Males compared to\nall females" = "GENDER",
                         "Males compared to\nfemales without HT" = "Groupmale-0", 
                         "Females with HT compared\nto females without HT" = "Groupfemale-1",
                         "Tau Summary" = "Tau_centered",
                         "Cortical Amyloid" = "PIB_centered",
                          "Tau x Female with HT Interaction" ="Tau_centered:Groupfemale-1",
                          "Tau x Male Interaction" = "Tau_centered:Groupmale-0",
                          "Amyloid x Female with HT Interaction" ="PIB_centered:Groupfemale-1",
                          "Amyloid x Male Interaction" = "PIB_centered:Groupmale-0"), colors = "Set1") + 
  theme(legend.position = "none", axis.text.y = element_text(size = 8)) + ggtitle("PACC Performance")

p2<-plot_summs(model11, scale = TRUE, plot.distributions = FALSE, 
           coefs = c("Tau_centered" = "Centered Tau Signature",
                     "Participant Age" = "Age",
                     "Years Education" = "EDUC",
                     "APOE4 Status" = "apoe4",
                     "Males compared to\nall females" = "GENDER",
                     "Males compared to\nfemales without HT" = "Groupmale-0", 
                     "Females with HT compared\nto females without HT" = "Groupfemale-1",
                     "Tau Summary" = "Tau_centered",
                     "Cortical Amyloid" = "PIB_centered",
                     "Tau x Female with HT Interaction" ="Tau_centered:Groupfemale-1",
                     "Tau x Male Interaction" = "Tau_centered:Groupmale-0",
                     "Amyloid x Female with HT Interaction" ="PIB_centered:Groupfemale-1",
                     "Amyloid x Male Interaction" = "PIB_centered:Groupmale-0"), colors = "#377EB8") + 
  theme(legend.position = "none", axis.text.y = element_text(size = 8)) + ggtitle("PACC Performance")

grid.arrange(p1, p2, nrow = 2)



df.meds$Group<-as.factor(df.meds$Group)
levels(df.meds$Group)<-c("Female, no HT", "Female, HT Use", "Male")

p1<-ggplot(df.meds, aes(x = PIB_fSUVR_rsf_TOT_CORTMEAN, y = PACC, colour = Group, shape = Group, group = Group)) + geom_point() + 
  geom_smooth(method = "lm", fullrange = TRUE) + theme(legend.position = "bottom",panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), 
                                                       axis.line = element_line(colour = "black"))+ 
  scale_color_manual(values = c("#d7191c", "#fdae61", "#2b83ba")) +
  scale_shape_manual(values = c(13, 1, 3)) + xlab("Cortical Amyloid Burden")

p2<-ggplot(df.meds, aes(x = TauSig, y = PACC, colour = Group, shape = Group, group = Group)) + geom_point() + 
  geom_smooth(method = "lm", fullrange = TRUE) + theme(legend.position = "bottom",panel.grid.major = element_blank(), 
                                                       panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), 
                                                       axis.line = element_line(colour = "black"))+ 
  scale_color_manual(values = c("#d7191c", "#fdae61", "#2b83ba")) +
  scale_shape_manual(values = c(13, 1, 3)) + xlab("Summary Tau Burden")
grid.arrange(p1, p2, nrow = 2)

df.demogs<-df.meds

df.demogs[,3:13]<-exp(df.demogs[,3:13])
df.demogs$PIBpos<-ifelse(df.demogs$PIB_fSUVR_rsf_TOT_CORTMEAN > 1.42, 1, 0)

#Demographics table for all participants

myVars <- c("Age", "EDUC", "apoe4", "race2",  "TimeBetween", "PIB_fSUVR_rsf_TOT_CORTMEAN", "PIBpos", 
            "TAU_fSUVR_rsf_TOT_CTX_CUNEUS", "TAU_fSUVR_rsf_TOT_CTX_FRNPOLE", "TAU_fSUVR_rsf_TOT_CTX_INFERPRTL",
            "TAU_fSUVR_rsf_TOT_CTX_LATOCC", "TAU_fSUVR_rsf_TOT_CTX_MIDTMP", "TAU_fSUVR_rsf_TOT_CTX_PARSORBLS",
            "TAU_fSUVR_rsf_TOT_CTX_PARSTRNGLS","TAU_fSUVR_rsf_TOT_CTX_ROSMIDFRN",  "TAU_fSUVR_rsf_TOT_CTX_SUPERTMP",
            "TAU_fSUVR_rsf_TOT_CTX_TMPPOLE","Tauopathy")
catVars <- c("apoe4",  "race2", "PIBpos")
CreateTableOne(vars = myVars, data = df.demogs, factorVars = catVars, strata = c("Group"))

