#Read the data
library(readxl)
setwd("D:/ITC/Q2.2/Project")
my_data <- read_excel("D:/ITC/Q2.2/Project/landslide_data_absolute.xlsx", sheet =
                          "ALL")
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

par(mfrow=c(2,2))
truehist(data_scaled$IR)
truehist(IR.CLASS)
plot(Data.Mat$IR,data_scaled$IR)

#make Dataset1, which include addition information that will be used in addition to the linear covariates by the model
Dataset1 = data.frame(intercept = 1, Ntrials = 1, y = data_scaled$LS_Pre,Slope.CLASS, IR.CLASS,
                          LU.CLASS = Data.Mat$LU, Soil_Type.CLASS = Data.Mat$Soil_Type,Geology.CLASS = Data.Mat$Geology)

#make Dataset2, which include all linear covariates
match(c("ID","LS_Pre","LS_Maria","LS_Maria_Scarp","SlopeAVG","IR","LU","Soil_Type","Geology"),names(data_scaled))
Dataset2 = data_scaled[,-c(1,3,9,10,11,17,21,22,23)]

#Check
match(c("ID","LS_Pre","LS_Maria","LS_Maria_Scarp","SlopeAVG","IR","LU","Soil_Type","Geology"),names(Dataset2))
Dataset2 = as.matrix(Dataset2)

#We now need to set some parameters, or hyperparameters for the nonlinear covariates.
hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))

#declare the formula that will pass to INLA
Formula1.COURSE = y~ -1 + intercept + Dataset2 +
f(Slope.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal= 1E-4) +
f(IR.CLASS, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal=1E-4) +
f(LU.CLASS, model="iid",hyper=hyper.iid, constr=T) +
f(Soil_Type.CLASS, model="iid",hyper=hyper.iid, constr = T) +
f(Geology.CLASS, model = "iid",hyper = hyper.iid, constr = T)

#run the model and check plot
Mod1.COURSE = inla(formula = Formula1.COURSE, family = "binomial",data=c(as.list(Dataset1), list(Dataset2=Dataset2)),
                       control.fixed=list(prec=.1),
                       num.threads = 2,
                       Ntrials = Ntrials,
                       control.family = list(control.link = list(model = "logit")),
                       control.mode = list(restart = T),
                       control.inla = list(int.strategy = "eb"),
                       control.predictor = list(compute = T, link = 1), verbose = T)
plot(Mod1.COURSE)

#pass INLA with slope and IR as categorical instead of ordinals
Formula2.COURSE = y~ -1 + intercept + Dataset2 +
  f(Slope.CLASS, model="iid",hyper=hyper.iid, constr=T) +
  f(IR.CLASS, model="iid",hyper=hyper.iid, constr=T) +
  f(LU.CLASS, model="iid",hyper=hyper.iid, constr=T) +
  f(Soil_Type.CLASS, model="iid",hyper=hyper.iid, constr = T) +
  f(Geology.CLASS, model = "iid",hyper = hyper.iid, constr = T)
Mod2.COURSE = inla(formula = Formula2.COURSE, family = "binomial",data=c(as.list(Dataset1), list(Dataset2=Dataset2)),
                      control.fixed=list(prec=.1),
                      num.threads = 2,
                      Ntrials = Ntrials,
                      control.family = list(control.link = list(model = "logit")),
                      control.mode = list(restart = T),
                      control.inla = list(int.strategy = "eb"),
                      control.predictor = list(compute = T, link = 1), verbose = T)
plot(Mod2.COURSE)

#Examine the actual numerical output
Mod1.COURSE$summary.fixed

#check only the significant covariates
extractSignificantEffects=function(results.df){
 print(results.df[sign(results.df$`0.025quant`)==sign(results.df$`0.975quant`),])
}

#which looks for the 2.5 and 97.5 percentiles of the regression coefficients
 #and extracts those that share the same sign.
extractSignificantEffects(Mod1.COURSE$summary.fixed)

#check which is model is better from graphs
Slope.Class.Mod1 = Mod1.COURSE$summary.random$Slope.CLASS
Slope.Class.Mod2 = Mod2.COURSE$summary.random$Slope.CLASS
IR.Class.Mod1 = Mod1.COURSE$summary.random$IR.CLASS
IR.Class.Mod2 = Mod2.COURSE$summary.random$IR.CLASS

par(mfrow=c(1,2))
plot(Slope.Class.Mod1$ID, Slope.Class.Mod1$mean, type = "l", col = "blue", lwd = 4,ylim = c(-1.5,1),
xlab = "Slope Classes", ylab = "Regression Coefficients")
lines(Slope.Class.Mod1$ID, Slope.Class.Mod1$`0.025quant`, col = "blue", lwd = 4)
lines(Slope.Class.Mod1$ID, Slope.Class.Mod1$`0.975quant`, col = "blue", lwd = 4)
points(Slope.Class.Mod2$ID, Slope.Class.Mod2$mean, col = "red", pch = 19)
points(Slope.Class.Mod2$ID, Slope.Class.Mod2$`0.025quant`, col = "red", pch = 19)
points(Slope.Class.Mod2$ID, Slope.Class.Mod2$`0.975quant`, col = "red", pch = 19)

plot(IR.Class.Mod1$ID,IR.Class.Mod1$mean, type = "l", col = "blue", lwd = 4, ylim = c(-1.5,1),
         xlab = "IR Classes", ylab = "Regression Coefficients")
lines(IR.Class.Mod1$ID, IR.Class.Mod1$`0.025quant`, col = "blue", lwd = 4)
lines(IR.Class.Mod1$ID, IR.Class.Mod1$`0.975quant`, col = "blue", lwd = 4)
points(IR.Class.Mod2$ID, IR.Class.Mod2$mean, col = "red", pch = 19)
points(IR.Class.Mod2$ID, IR.Class.Mod2$`0.025quant`, col = "red", pch = 19)
points(IR.Class.Mod2$ID, IR.Class.Mod2$`0.975quant`, col = "red", pch = 19)

#make plots for regression coefficient
FixedEffects = data.frame(Mod1.COURSE$summary.fixed)
Dataset3 = data_scaled[,-c(1,3,9,10,11,17,21,22,23)]
LabelFixed = names(Dataset3)

par(mfrow=c(3,2))
plot(1:14,FixedEffects$mean[2:15],
         pch = 19, cex =1, col = "blue", ylim = c(-1.5,1.5),
         ylab = "Regression Coefficient", xaxt = 'n', xlab = "")
axis(side = 1, at = seq(1,14, by = 1),labels = LabelFixed,
         las =3)
abline(h = 0, col = "grey", lty = 2)
points(1:14,FixedEffects$X0.025quant[2:15],pch = 19, col = "black")
points(1:14,FixedEffects$X0.975quant[2:15],pch = 19, col = "black")
segments(1:14,FixedEffects$X0.025quant[2:15],1:14,FixedEffects$X0.975quant[2:15],
             lwd=2)

Slope.Class.Mod1 = Mod1.COURSE$summary.random$Slope.CLASS
plot(Slope.Class.Mod1$ID, Slope.Class.Mod1$mean, type = "l", col = "blue",
         lwd = 4, ylim = c(-1.5,1.5),xlab = "Slope Classes", ylab = "Regression Coefficients")
lines(Slope.Class.Mod1$ID, Slope.Class.Mod1$`0.025quant`, col = "grey", lwd = 4)
lines(Slope.Class.Mod1$ID, Slope.Class.Mod1$`0.975quant`, col = "grey", lwd = 4)
abline(h = 0, col = "grey", lty = 2)


IR.Class.Mod1 = Mod1.COURSE$summary.random$IR.CLASS
plot(IR.Class.Mod1$ID, IR.Class.Mod1$mean, type = "l", col = "blue",
         lwd = 4, ylim = c(-1.5,1.5),xlab = "IR Classes", ylab = "Regression Coefficients")
lines(IR.Class.Mod1$ID, IR.Class.Mod1$`0.025quant`, col = "grey", lwd = 4)
lines(IR.Class.Mod1$ID, IR.Class.Mod1$`0.975quant`, col = "grey", lwd = 4)
abline(h = 0, col = "grey", lty = 2)

LU.Class.Mod1 = Mod1.COURSE$summary.random$LU.CLASS
plot(1:8,LU.Class.Mod1$mean,
         pch = 19, cex =1, col = "blue", ylim = c(-1.5,1.5),
         ylab = "Regression Coefficient",xlab = "LU Classes")
abline(h = 0, col = "grey", lty = 2)
points(1:8,LU.Class.Mod1$`0.025quant`,pch = 19, col = "black")
points(1:8,LU.Class.Mod1$`0.975quant`,pch = 19, col = "black")
segments(1:8,LU.Class.Mod1$`0.025quant`,1:8,LU.Class.Mod1$`0.975quant`,
             lwd=2)

Geology.Class.Mod1 = Mod1.COURSE$summary.random$Geology.CLASS
plot(1:20,Geology.Class.Mod1$mean,
         pch = 19, cex =1, col = "blue", ylim = c(-1.5,1.5),
         ylab = "Regression Coefficient",xlab = "Geology Classes")
abline(h = 0, col = "grey", lty = 2)
points(1:20,Geology.Class.Mod1$`0.025quant`,pch = 19, col = "black")
points(1:20,Geology.Class.Mod1$`0.975quant`,pch = 19, col = "black")
segments(1:20,Geology.Class.Mod1$`0.025quant`,1:20,Geology.Class.Mod1$`0.975quant`,
             lwd=2)

Soil.Class.Mod1 = Mod1.COURSE$summary.random$Soil_Type.CLASS
plot(1:15,Soil.Class.Mod1$mean,
         pch = 19, cex =1, col = "blue", ylim = c(-1.5,1.5),
         ylab = "Regression Coefficient",xlab = "Soil Classes")
abline(h = 0, col = "grey", lty = 2)
points(1:15,Soil.Class.Mod1$`0.025quant`,pch = 19, col = "black")
points(1:15,Soil.Class.Mod1$`0.975quant`,pch = 19, col = "black")
segments(1:15,Soil.Class.Mod1$`0.025quant`,1:15,Soil.Class.Mod1$`0.975quant`,
             lwd=2)


#check goodness of fit, and first check probabilities
plot(density(Mod1.COURSE$summary.fitted.values$mean))

#check with boxplot
boxplot(Mod1.COURSE$summary.fitted.values$mean~data_scaled$LS_Pre, col = c("green","red"), ylab = "Susceptibility")


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

##make model for fraction method
setwd("D:/ITC/Q2.2/Project")
my_dataF <- read_excel("D:/ITC/Q2.2/Project/landslide_data.xlsx", sheet = "ALL")
Data.MatF = my_dataF

data_scaledF=Data.MatF
vars2scale=c("Area","SlopeAVG", "SlopeSD","PlanAVG","PlanSD","ProfileAVG",
                 "ProfileSD","BuildCuts","NorthAVG","NorthSD","EastAVG","EastSD",
                 "IR","RoadCuts","DEMAVG","DEMSD")

data_scaledF[,vars2scale]=apply(data_scaledF[,vars2scale],2,scale)

Slope.CLASSF = inla.group(data_scaledF$SlopeAVG, n = 20)
IR.CLASSF = inla.group(data_scaledF$IR, n = 20)
Dataset1F = data.frame(intercept = 1, Ntrials = 1, y = data_scaledF$LS_Pre, Slope.CLASSF, IR.CLASSF,
                           LU.CLASSF = Data.MatF$LU, Soil_Type.CLASSF = Data.MatF$Soil_Type, Geology.CLASSF = Data.MatF$Geology)

match(c("ID","LS_Pre","LS_Maria","LS_Maria_Scarp","SlopeAVG","IR","LU","Soil_Type","Geology"),names(data_scaledF))
Dataset2F = data_scaledF[,-c(1,3,9,10,11,17,21,22,23)]
Dataset2F = as.matrix(Dataset2F)

hyper.rw = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))
hyper.iid = list(theta1 = list(prior="pc.prec", param=c(0.1, 0.5)))

Formula1F.COURSE = y~ -1 + intercept + Dataset2F +
 f(Slope.CLASSF, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
 f(IR.CLASSF, model="rw1", hyper=hyper.rw, constr=T, scale.model = TRUE, diagonal = 1E-4) +
 f(LU.CLASSF, model="iid",hyper=hyper.iid, constr=T) +
 f(Soil_Type.CLASSF, model="iid",hyper=hyper.iid, constr = T) +
 f(Geology.CLASSF, model = "iid",hyper = hyper.iid, constr = T)

Mod1F.COURSE = inla(formula = Formula1F.COURSE, family = "binomial", data=c(as.list(Dataset1F), list(Dataset2=Dataset2F)),
                        control.fixed=list(prec=.1),
                        num.threads = 2,
                        Ntrials = Ntrials,
                        control.family = list(control.link = list(model = "logit")),
                        control.mode = list(restart = T),
                        control.inla = list(int.strategy = "eb"),
                        control.predictor = list(compute = T, link = 1), verbose = T)
plot(Mod1F.COURSE)


#check for ROC of absolute method
library(pROC)
SusceptibilityMod1 = Mod1.COURSE$summary.fitted.values$mean

library(verification)

Test.ROC = roc.plot(data_scaled$LS_Maria, SusceptibilityMod1, threshold = seq(0,1,0.1))
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

#check ROC for fraction method
SusceptibilityMod1F = Mod1F.COURSE$summary.fitted.values$mean
Test.ROCF = roc.plot(data_scaledF$LS_Maria, SusceptibilityMod1F, threshold = seq(0,1, 0.1))
Test.ROCF$plot.data
Test.ROC_scarpF = roc.plot(data_scaledF$LS_Maria_Scarp, SusceptibilityMod1F, threshold = seq(0,1, 0.1))
Test.ROC_scarpF$plot.data
Test.ROC_preF = roc.plot(data_scaledF$LS_Pre, SusceptibilityMod1F, threshold = seq(0,1, 0.1))
Test.ROC_preF$plot.data

ROC.predF = roc(data_scaledF$LS_Maria~SusceptibilityMod1F)
ROC.predF
ROC.pred_scarpF = roc(data_scaledF$LS_Maria_Scarp~SusceptibilityMod1F)
ROC.pred_scarpF
ROC.pred_preF = roc(data_scaledF$LS_Pre~SusceptibilityMod1F)
ROC.pred_preF

#Susceptibility mapping
#calculate uncertainty for each fitted values
SusceptibilityUncertainty = Mod1.COURSE$summary.fitted.values$`0.975quant`-
Mod1.COURSE$summary.fitted.values$`0.025quant`

#make the final dataset and write the dataset in the file
data4map = cbind(Data.Mat$ID,
                     Mod1.COURSE$summary.fitted.values$mean,
                     SusceptibilityUncertainty)
colnames(data4map) <- c("ID", "MeanSusceptibility", "SusceptibilityUncertainty")
#write.table(data4map, "Susceptibility4GIS.txt", sep = "\t", col.names = T, row.names = F)

#plot the fitted vallues with its uncertainty
par(mfrow=c(1,1))
plot(density(Mod1.COURSE$summary.fitted.values$mean))
plot(Mod1.COURSE$summary.fitted.values$mean,SusceptibilityUncertainty, pch = 19)
smoothScatter(Mod1.COURSE$summary.fitted.values$mean,SusceptibilityUncertainty)
#dev.off()
#plot roc
par(mfrow = c(1,3))
plot(1-ROC.pred$specificities,ROC.pred$sensitivities,
         col = "red",type = "l", xlab = "FPR", ylab = "TPR")
lines(1-ROC.predF$specificities,ROC.predF$sensitivities,col = "blue", lwd=2)
plot(1-ROC.pred$specificities,ROC.pred$sensitivities,
         col = "red",type = "l", xlab = "FPR", ylab = "TPR")
lines(1-ROC.pred_scarp$specificities,ROC.pred_scarp$sensitivities,
          col = "green",lwd=2)
plot(1-ROC.pred_scarp$specificities,ROC.pred_scarp$sensitivities,
         col = "green",type = "l", xlab = "FPR", ylab = "TPR")
lines(1-ROC.pred_pre$specificities, ROC.pred_pre$sensitivities,
          col = "black",lwd=2)
