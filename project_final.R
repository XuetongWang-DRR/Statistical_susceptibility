#Read the data
library(readxl)
setwd("D:/ITC/Q2.2/Project")

my_data <- read_excel("D:/ITC/Q2.2/Project/landslide_data_absolute.xlsx", sheet = "ALL")
Data.Mat = my_data

#load packages of INLA and MASS
library(INLA)
library(MASS)

#Make histgram plots for some factors before rescaled
par(mfrow=c(3,4))
truehist(Data.Mat$SlopeAVG, xlab = "Mean Slope [degree]") 
truehist(Data.Mat$PlanAVG, xlab = "Mean Plan Curvature [1/m]")
truehist(Data.Mat$ProfileAVG, xlab = "Mean Profile Curvature [1/m]")
truehist(Data.Mat$LU, xlab = "Landuse Class")
truehist(Data.Mat$Soil_Type, xlab = "Soil_Type Class")
truehist(Data.Mat$Geology, xlab = "Geology Class")
truehist(Data.Mat$BuildCuts, xlab = "BuildCuts Area [meter^2]")
truehist(Data.Mat$NorthAVG, xlab = "mean Northness")
truehist(Data.Mat$EastAVG, xlab = "mean Eastness")
truehist(Data.Mat$IR, xlab = "Internal Relief [m/km^2]")
truehist(Data.Mat$RoadCuts, xlab = "Road Cuts Area [meter^2]")
truehist(Data.Mat$DEMAVG, xlab = "mean elevation [meter]")


#Rescaled factors
data_scaled=Data.Mat
vars2scale=c("Area","SlopeAVG", "SlopeSD","PlanAVG","PlanSD","ProfileAVG",
             "ProfileSD","BuildCuts","NorthAVG","NorthSD","EastAVG","EastSD",
             "IR","RoadCuts","DEMAVG","DEMSD")  
data_scaled[,vars2scale]=apply(data_scaled[,vars2scale],2,scale)

#change some numerical factors to categorical factors
Slope.CLASS = inla.group(data_scaled$SlopeAVG, n = 20)
IR.CLASS = inla.group(data_scaled$IR, n = 20)

#check rescaled and categorical covariate
par(mfrow=c(2,2))
truehist(data_scaled$SlopeAVG)
truehist(Slope.CLASS)
plot(Data.Mat$SlopeAVG,data_scaled$SlopeAVG)

#make Dataset1, which include addition information that will be used in addition to the linear covariates by the model
Dataset1 = data.frame(intercept = 1, Ntrials = 1, y = data_scaled$LS_Pre, Slope.CLASS, IR.CLASS,
                      LU.CLASS = Data.Mat$LU, Soil_Type.CLASS = Data.Mat$Soil_Type, Geology.CLASS = Data.Mat$Geology)

#make Dataset2, which include all linear covariates
match(c("ID","LS_Pre","LS_Maria","LS_Maria_Scarp","SlopeAVG","IR","LU","Soil_Type","Geology"),names(data_scaled))
Dataset2 = data_scaled[,-c(1,3,9,10,11,17,21,22,23)]
#Check
match(c("ID","LS_Pre","LS_Maria","LS_Maria_Scarp","SlopeAVG","IR","LU","Soil_Type","Geology"),names(Dataset2))
Dataset2 = as.matrix(Dataset2)

#We now need to set some parameters, or hyperparameters for the nonlinear covariates. 
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))

#make the formula for model, and apply it with INLA function
Formula1.COURSE = y~ -1 + intercept + Dataset2 +
  f(Slope.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
  f(IR.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) + 
  f(LU.CLASS, model="iid",hyper=hyper.iid, constr=T) +
  f(Soil_Type.CLASS, model="iid",hyper=hyper.iid, constr = T) +
  f(Geology.CLASS, model = "iid",hyper = hyper.iid, constr = T)


Mod1.COURSE = inla(formula = Formula1.COURSE, family = "binomial", data=c(as.list(Dataset1), list(Dataset2=Dataset2)),
                   control.fixed=list(prec=.1),
                   num.threads = 2, 
                   Ntrials = Ntrials,
                   control.family = list(control.link = list(model = "logit")), 
                   control.mode = list(restart = T),
                   control.inla = list(int.strategy = "eb"), 
                   control.predictor = list(compute = T, link = 1), verbose = T)

#check
check = cbind(data_scaled$LS_Pre,data_scaled$LS_Maria,Mod1.COURSE$summary.fitted.values)
check1 = cbind(data_scaled$LS_Pre,data_scaled$LS_Maria,Mod1.COURSE$summary.fitted.values$mean)
check_scarp = cbind(data_scaled$LS_Pre,data_scaled$LS_Maria_Scarp,Mod1.COURSE$summary.fitted.values)
check1_scarp = cbind(data_scaled$LS_Pre,data_scaled$LS_Maria_Scarp,Mod1.COURSE$summary.fitted.values$mean)

#simple check with boxplot
boxplot(Mod1.COURSE$summary.fitted.values$mean~data_scaled$LS_Maria, 
        col = c("green", "red"), ylab = "Susceptibility")
boxplot(Mod1.COURSE$summary.fitted.values$mean~data_scaled$LS_Maria_Scarp, 
        col = c("green", "red"), ylab = "Susceptibility")

#check for ROC
library(pROC)
SusceptibilityMod1 = Mod1.COURSE$summary.fitted.values$mean

library(verification)
Test.ROC = roc.plot(data_scaled$LS_Maria, SusceptibilityMod1, threshold = seq(0,1, 0.1))
Test.ROC$plot.data
Test.ROC_scarp = roc.plot(data_scaled$LS_Maria_Scarp, SusceptibilityMod1, threshold = seq(0,1, 0.1))
Test.ROC_scarp$plot.data
Test.ROC_pre = roc.plot(data_scaled$LS_Pre, SusceptibilityMod1, threshold = seq(0,1, 0.1))
Test.ROC_pre$plot.data

ROC.pred = roc(data_scaled$LS_Maria~SusceptibilityMod1)
ROC.pred
ROC.pred_scarp = roc(data_scaled$LS_Maria_Scarp~SusceptibilityMod1)
ROC.pred_scarp
ROC.pred_pre = roc(data_scaled$LS_Pre~SusceptibilityMod1)
ROC.pred_pre 

#plot roc
plot(roc(data_scaled$LS_Maria~SusceptibilityMod1))
plot(1-ROC.pred$specificities,ROC.pred$sensitivities,
     type = "l", xlab = "1-Specificity", ylab = "Sensitivity")
plot(roc(data_scaled$LS_Maria_Scarp~SusceptibilityMod1))
plot(1-ROC.pred$specificities,ROC.pred$sensitivities,
     type = "l", xlab = "1-Specificity", ylab = "Sensitivity")

#Susceptibility mapping
#calculate uncertainty for each fitted values
SusceptibilityUncertainty = Mod1.COURSE$summary.fitted.values$`0.975quant`-
  Mod1.COURSE$summary.fitted.values$`0.025quant`

#make the final dataset and write the dataset in the file
data4map = cbind(Data.Mat$ID,
                 Mod1.COURSE$summary.fitted.values$mean,
                 SusceptibilityUncertainty)
colnames(data4map) <- c("ID", "MeanSusceptibility", "SusceptibilityUncertainty")
write.table(data4map, "Susceptibility4FRACTION.txt", sep = "\t", col.names = T, row.names = F)

#plot the fitted vallues with its uncertainty
par(mfrow=c(1,1))
plot(density(Mod1.COURSE$summary.fitted.values$mean))
plot(Mod1.COURSE$summary.fitted.values$mean,SusceptibilityUncertainty, pch = 19)
smoothScatter(Mod1.COURSE$summary.fitted.values$mean,SusceptibilityUncertainty)
dev.off()