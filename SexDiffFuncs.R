AssignNames<-function(DF, NAMES){
  rownames(DF)<-NAMES
  colnames(DF)<-NAMES
  return(DF)
}

SetUpBarCharts<-function(df){
  AmyloidTau_p<-matrix(, nrow = 88, ncol = 88)
  AmyloidTau_r<-matrix(, nrow = 88, ncol = 88)
  
  for(i in 1:88){
    for(j in 1:88){
      AmyloidTau_p[i, j]<-cor.test(df[,18+i], df[,137+j], method = "pearson")$p.value
      AmyloidTau_r[i, j]<-cor.test(df[,18+i], df[,137+j], method = "pearson")$estimate
    }}
  
  NAMES<-names(df[138:225])
  NAMES<-substring(NAMES, 15)
  NAMES = substr(NAMES,1,nchar(NAMES)-2)
  
  AmyloidTau_p<-AssignNames(AmyloidTau_p, NAMES)
  AmyloidTau_r<-AssignNames(AmyloidTau_r, NAMES)
  
  
  
  df.new<-data.frame("region" = character(), "Comparison1" = character(), 
                     "corr1" = numeric(), "p1" = numeric(), stringsAsFactors = FALSE)
  for(i in 1:88){
    df.new[i, 1]<-NAMES[i]
    df.new[i, 2]<-"AmyloidTau"
    df.new[i, 3]<-AmyloidTau_r[i, i]
    df.new[i, 4]<-AmyloidTau_p[i, i]
  }
  df.new$AmyloidTauSignificance<-as.factor(ifelse(df.new$p1 < 0.001, 1, 0))
  return(df.new)}


FancyPlot<-function(df, i){
  df.new.f<-SetUpBarCharts(df[df$GENDER == "female",])
  df.new.m<-SetUpBarCharts(df[df$GENDER == "male",])
  
  pMain<-ggplot(data = df, aes(x = df[,i], y = df[,i+119], colour= GENDER, shape = GENDER))+
    geom_point() + geom_smooth(method = "lm", aes(group = GENDER)) +xlab(names(df[i])) + 
    ylab(names(df[i+119]))+theme_classic()+xlim(c(1, 3))+ylim(c(0, 3.3))+
    annotate("label", x = 1.5, y = 3, size = 3,
             label = paste("Female \n r=", round(df.new.f[i+28, 3], 3), "\n p=", formatC(df.new.f[i+28,4], format = "e", digits = 2)))+
    annotate("label", x = 2.5, y = 3, size = 3,
             label = paste("Male \n r=", round(df.new.m[i+28, 3], 3), "\n p=", formatC(df.new.m[i+28,4], format = "e", digits = 2))) 
  
  
  pTop<-ggplot(df, aes(x=df[,i], fill=GENDER)) +
    geom_density(alpha=0.4, adjust = 2)+xlab("")+ggtitle("Amyloid")+
    theme_classic()+ guides(fill=FALSE)
  pRight<-ggplot(df, aes(x=df[,i+119], fill=GENDER)) +
    geom_density(alpha=0.4, adjust = 2)+xlab("")+ggtitle("Tau")+
    theme_classic()+ guides(fill=FALSE) +coord_flip()
  
  pEmpty <- ggplot(data = df, aes(x = df[,i], y = df[,i+119])) +
    geom_blank() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          line = element_blank(),
          panel.background = element_blank())
  
  plot<-grid.arrange(pTop, pEmpty, pMain, pRight,
                     ncol = 2, nrow = 2, widths = c(3, 1), heights = c(1, 3))
  return(plot)}

TestCSFSignificance<-function(DF){
  result<-data.frame("Region" = character(), "Amyloid_HippoControl" = numeric(), "VolRofC_HippoControl" = numeric(), "Sex_HippoControl" = numeric(),
                     "Amyloid_CerebControl" = numeric(), "VolRofC_CerebControl" = numeric(), "Sex_CerebControl" = numeric(), stringsAsFactors = FALSE)
  result[1, 1]<-"CSF"
  Model.null<-lm(E_ptau ~ E_ab42 + VolRateofChange + GENDER + apoe4 + EDUC, data = DF)
  Model.new<-lm(E_ptau ~  VolRateofChange + GENDER + apoe4 + EDUC, data = DF)
  result[1, 2]<-anova(Model.null, Model.new)$`Pr(>F)`[2]
  Model.new<-lm(E_ptau ~ E_ab42 +  GENDER + apoe4 + EDUC, data = DF)
  result[1, 3]<-anova(Model.null, Model.new)$`Pr(>F)`[2]
  Model.new<-lm(E_ptau ~ E_ab42 + VolRateofChange + apoe4 + EDUC, data = DF)
  result[1, 4]<-anova(Model.null, Model.new)$`Pr(>F)`[2]
  
  Model.null<-lm(E_ptau ~ E_ab42 + CerebVolRateofChange + GENDER + apoe4 + EDUC, data = DF)
  Model.new<-lm(E_ptau ~  CerebVolRateofChange + GENDER + apoe4 + EDUC, data = DF)
  result[1, 5]<-anova(Model.null, Model.new)$`Pr(>F)`[2]
  Model.new<-lm(E_ptau ~ E_ab42 +  GENDER + apoe4 + EDUC, data = DF)
  result[1, 6]<-anova(Model.null, Model.new)$`Pr(>F)`[2]
  Model.new<-lm(E_ptau ~ E_ab42 + CerebVolRateofChange + apoe4 + EDUC, data = DF)
  result[1, 7]<-anova(Model.null, Model.new)$`Pr(>F)`[2]
  return(result)
}



MakePlot<-function(DF, REGION, TITLE){
  p<-ggplot(DF, aes_string(x = paste("PUP_fSUVR_rsf_", REGION, ".x", sep = ""), y = paste("PUP_fSUVR_rsf_", REGION, ".y", sep = ""), group = "GENDER", colour = "GENDER", shape = "GENDER")) + geom_point() + 
    geom_smooth(method = "lm")+ scale_colour_manual("sex", values = c( "#F8766D", "#00BFC4")) + 
    scale_shape_manual("sex", values = c(2, 1)) + ylab("Tau")+ xlab("Amyloid") + theme(legend.position="none")+ 
    guides(fill=guide_legend(title="sex")) + ggtitle(TITLE)
  return(p)}

Boxplot<-function(NETWORK, NETWORKTITLE){
 p<- ggplot(FC.all, aes(x=GENDER, y=NETWORK, fill=GENDER)) +  
    geom_boxplot(alpha = 0.8, outlier.shape = NA) + geom_point(aes(fill = GENDER, shape = GENDER), size = 1,  position = position_jitterdodge()) +
    xlab("Sex") + ylab(NETWORKTITLE)+ scale_shape_manual(values=c(21, 24))+theme(legend.position = "none")
  return(p)
}


MakeBarPlots<-function(Component, ComponentTitle){
  theme_set(theme_bw())  
  pcaPlot<-as.data.frame(Component)
  pcaPlot<-cbind(rownames(pcaPlot), pcaPlot[,1])
  colnames(pcaPlot)<-c("region", "contribution")
  pcaPlot<-as.data.frame(pcaPlot)
  pcaPlot$contribution<-as.numeric(as.character(pcaPlot$contribution))
  pcaPlot <- pcaPlot[order(-pcaPlot$contribution), ]  # sort
  # Diverging Barcharts
  p<-ggplot(pcaPlot, aes(x=reorder(region, contribution), y = contribution, label=contribution)) + xlab("") +
    geom_bar(stat='identity', width=.5)  +
    #coord_cartesian(ylim = c(-0.4, 0.4)) +
    scale_y_continuous(limits = c(0, 0.125))+
    labs(subtitle="28% of variance", 
         title= ComponentTitle) + 
    coord_flip() + 
    theme(legend.position = "none")
  return(p)}
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev}
contrib <- function(var.cos2, comp.cos2){var.cos2/comp.cos2}

#Making Figure 1 - all the regions where women have more tau than men
GenerateProjectionPlot<-function(index, MODEL, XLAB = "You forgot", YLAB = "You Forgot"){
  df.predict <- cbind(df.log[,c("GENDER", "Age", "apoe4", "EDUC", names(df.log[index]), names(df.log[index+92]))], predict(Model1, interval = 'confidence'))
  dummydata1<-data.frame("Tau" = c(seq(from = min(df.log[,index]), to = max(df.log[,index]), length.out = round(length(df.log[,index])/2, 0))),
                         "Amyloid" = c(seq(from = min(df.log[,index+92]), to = max(df.log[,index+92]), length.out = round(length(df.log[,index+92])/2, 0))),
                         "apoe4" = rep(1, round(length(df.log[,index])/2, 0)), 
                         "EDUC" = rep(16, round(length(df.log[,index])/2, 0)), 
                         "Age" = rep(mean(df.log$Age), round(length(df.log[,index])/2, 0)), 
                         "GENDER" = rep("male", round(length(df.log[,index])/2, 0)) )  
  dummydata2<-data.frame("Tau" = c(seq(from = min(df.log[,index]), to = max(df.log[,index]), length.out = -1+round(length(df.log[,index])/2, 0))),
                         "Amyloid" = c(seq(from = min(df.log[,index+92]), to = max(df.log[,index+92]), length.out = -1+round(length(df.log[,index+92])/2, 0))),
                         "apoe4" = rep(1, -1+round(length(df.log[,index])/2, 0)), 
                         "EDUC" = rep(16, -1+round(length(df.log[,index])/2, 0)), 
                         "Age" = rep(mean(df.log$Age), -1+round(length(df.log[,index])/2, 0)), 
                         "GENDER" = rep("female", -1+round(length(df.log[,index])/2, 0)) )  
  dummydata<-data.frame(rbind(dummydata1, dummydata2))
  dummydata$apoe4<-as.factor(dummydata$apoe4)
  rm(dummydata1, dummydata2)
  
  df.predict<-data.frame(cbind(df.log$GENDER, df.log[,index], df.log[,index+92], dummydata$Tau, dummydata$Amyloid, dummydata$GENDER,
                               predict(MODEL, dummydata, interval = 'confidence')))
  colnames(df.predict)<-c("ActualGender", "ActualTau", "ActualAmyloid", "ProjectedTau", "ProjectedAmyloid",
                          "ProjectedGender", "fit", "lwr", "upr")
  df.predict$ActualGender<-ifelse(df.predict$ActualGender == 1, "female", "male")
  df.predict$ProjectedGender<-ifelse(df.predict$ProjectedGender == 1, "male", "female")
  p<-ggplot(df.predict, aes(x = ActualAmyloid, y = ActualTau, color = ActualGender, shape = ActualGender)) + geom_point() +
    geom_smooth(aes(x = ProjectedAmyloid, y = fit,  color = ProjectedGender, shape = ProjectedGender), method = "lm") + scale_shape_manual(values=c(1, 2))+
    labs(color = "Gender", shape = "Gender")+xlab(XLAB) + ylab(YLAB)
  return(p)}

PACCplot<-function(REGION, XLAB){
  p<-ggplot(df.log, aes_string(x = REGION, y = "PACC", group = "GENDER", colour = "GENDER", shape = "GENDER")) + geom_point() + 
    geom_smooth(method = "lm")+ scale_colour_manual("sex", values = c( "#F8766D", "#00BFC4")) + 
    scale_shape_manual("sex", values = c(2, 1)) + ylab("PACC")+ xlab(XLAB) + theme(legend.position="none")+ 
    guides(fill=guide_legend(title="sex"))
  return(p)}



