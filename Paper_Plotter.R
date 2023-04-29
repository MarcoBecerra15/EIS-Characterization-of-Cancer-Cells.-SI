#                                                    DiazLab Microfluidics Laboratory
#                                                               Malcom Díaz García
#                                                                 05/May/2022
#                                                    University of Puerto Rico at Mayaguez 
#                                                                 Any Cell
# Main Functions:
#   1) Run Random Forest Model and Plot based on given resolution and dimensions for paper  
#   2) Save Plot in Directory
#   3) Make Raw Data Plots
#   4) Save Raw Data Plots in Directory
#-------LIBRARIES-------
library (cowplot)
library (randomForest)
library (gghighlight)
library(RColorBrewer)
library(xlsx)
library(ggplot2)
library(R.matlab)
library(reshape2)
library(gridExtra)
library(greekLetters)
library(lemon)
library(mdatools)
library(prospectr)
library(ellipse)
library(plotrix)
library(tiger)
library(utils)
library(base)
library(stats)
library(ChemoSpec)
library(HotellingEllipse)
library(purrr)
library(ggforce)
#----------------RF DATA Plotting------------------------
#----------------PLOT 1 ----------------------------

setwd("D:/Data Review Thesis/Random Forests/HELAMCFMix") 
location = 'ZmodRb_EIS_HELA_MCF_60'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)
char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#Recall RF Function with Optimal mtry & mtree
#optional
set.seed(42) 
model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                     importance=TRUE, keep.inbag=TRUE, mtry=12, ntree = 300)

model

distmatrix = as.dist(1-model$proximity)
scaleddistmtrx = cmdscale(distmatrix, eig = TRUE, x.ret = TRUE)
mdsaxisvar = round(abs(scaleddistmtrx$eig)/sum(abs(scaleddistmtrx$eig))*100,1)
mdsvals = scaleddistmtrx$points
Samples = rownames(mdsvals)
x = mdsvals[,1]
y = mdsvals[,2]
Cell = spectranf$Type
mdsdata = data.frame(Samples, x, y, Cell)


RFplot1 = ggplot(mdsdata, aes(x = x, y = y, label = Cell))+
  geom_point(aes(color= Cell, shape = Cell))+theme_bw()+
  xlab(paste("Axis 1 - ", mdsaxisvar[1], "%", sep = ""))+
  ylab(paste("Axis 2 - ", mdsaxisvar[2], "%", sep = ""))+theme(legend.position = "none")+
  scale_shape(solid = FALSE)+ggtitle("A.2")

RFplot1
#----------------PLOT 2 -------------------------------
setwd("D:/Data Review Thesis/Random Forests/HELAMCFMix") 
location = 'ZmodRb_EIS_HELA_MCF_90'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)
char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#Recall RF Function with Optimal mtry & mtree
#optional
set.seed(42) 
model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                     importance=TRUE, keep.inbag=TRUE, mtry=28, ntree = 300)

model

distmatrix = as.dist(1-model$proximity)
scaleddistmtrx = cmdscale(distmatrix, eig = TRUE, x.ret = TRUE)
mdsaxisvar = round(abs(scaleddistmtrx$eig)/sum(abs(scaleddistmtrx$eig))*100,1)
mdsvals = scaleddistmtrx$points
Samples = rownames(mdsvals)
x = mdsvals[,1]
y = mdsvals[,2]
Cell = spectranf$Type
mdsdata = data.frame(Samples, x, y, Cell)


RFplot2 = ggplot(mdsdata, aes(x = x, y = y, label = Cell))+
  geom_point(aes(color= Cell, shape = Cell))+theme_bw()+
  xlab(paste("Axis 1 - ", mdsaxisvar[1], "%", sep = ""))+
  ylab(paste("Axis 2 - ", mdsaxisvar[2], "%", sep = ""))+theme(legend.position = "none")+
  scale_shape(solid = FALSE)+ggtitle("A.3")

RFplot2
#----------------PLOT 3 -------------------------------
setwd("D:/Data Review Thesis/Random Forests/HELAMCFMix") 
location = 'ZphzRb_EIS_HELA_MCF_60'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)
char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#Recall RF Function with Optimal mtry & mtree
#optional
set.seed(42) 
model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                     importance=TRUE, keep.inbag=TRUE, mtry=28, ntree = 1000)

model

distmatrix = as.dist(1-model$proximity)
scaleddistmtrx = cmdscale(distmatrix, eig = TRUE, x.ret = TRUE)
mdsaxisvar = round(abs(scaleddistmtrx$eig)/sum(abs(scaleddistmtrx$eig))*100,1)
mdsvals = scaleddistmtrx$points
Samples = rownames(mdsvals)
x = mdsvals[,1]
y = mdsvals[,2]
Cell = spectranf$Type
mdsdata = data.frame(Samples, x, y, Cell)


RFplot3 = ggplot(mdsdata, aes(x = x, y = y, label = Cell))+
  geom_point(aes(color= Cell, shape = Cell))+theme_bw()+
  xlab(paste("Axis 1 - ", mdsaxisvar[1], "%", sep = ""))+
  ylab(paste("Axis 2 - ", mdsaxisvar[2], "%", sep = ""))+theme(legend.position = "none")+
  scale_shape(solid = FALSE)+ggtitle("B.2")

RFplot3
#----------------PLOT 4 -------------------------------
setwd("D:/Data Review Thesis/Random Forests/HELAMCFMix") 
location = 'ZphzRb_EIS_HELA_MCF_90'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)
char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#Recall RF Function with Optimal mtry & mtree
#optional
set.seed(42) 
model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                     importance=TRUE, keep.inbag=TRUE, mtry=28, ntree = 1000)

model

distmatrix = as.dist(1-model$proximity)
scaleddistmtrx = cmdscale(distmatrix, eig = TRUE, x.ret = TRUE)
mdsaxisvar = round(abs(scaleddistmtrx$eig)/sum(abs(scaleddistmtrx$eig))*100,1)
mdsvals = scaleddistmtrx$points
Samples = rownames(mdsvals)
x = mdsvals[,1]
y = mdsvals[,2]
Cell = spectranf$Type
mdsdata = data.frame(Samples, x, y, Cell)


RFplot4 = ggplot(mdsdata, aes(x = x, y = y, label = Cell))+
  geom_point(aes(color= Cell, shape = Cell))+theme_bw()+
  xlab(paste("Axis 1 - ", mdsaxisvar[1], "%", sep = ""))+
  ylab(paste("Axis 2 - ", mdsaxisvar[2], "%", sep = ""))+theme(legend.position = "none")+
  scale_shape(solid = FALSE)+ggtitle("B.3")

RFplot4

#----------------PLOT 5 -------------------------------
setwd("D:/Data Review Thesis/Random Forests/HELAMCFMix") 
location = 'ZmodRb_EIS_HELA_MCF_30'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)
char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#Recall RF Function with Optimal mtry & mtree
#optional
set.seed(42) 
model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                     importance=TRUE, keep.inbag=TRUE, mtry=28, ntree = 300)

model

distmatrix = as.dist(1-model$proximity)
scaleddistmtrx = cmdscale(distmatrix, eig = TRUE, x.ret = TRUE)
mdsaxisvar = round(abs(scaleddistmtrx$eig)/sum(abs(scaleddistmtrx$eig))*100,1)
mdsvals = scaleddistmtrx$points
Samples = rownames(mdsvals)
x = mdsvals[,1]
y = mdsvals[,2]
Cell = spectranf$Type
mdsdata = data.frame(Samples, x, y, Cell)


RFplot5 = ggplot(mdsdata, aes(x = x, y = y, label = Cell))+
  geom_point(aes(color= Cell, shape = Cell))+theme_bw()+
  xlab(paste("Axis 1 - ", mdsaxisvar[1], "%", sep = ""))+
  ylab(paste("Axis 2 - ", mdsaxisvar[2], "%", sep = ""))+theme(legend.position = "none")+theme(plot.margin = unit(c(5.5, 6.3, 5.5, 5.5), "pt"))+
  scale_shape(solid = FALSE)+ggtitle("A.1")

RFplot5
#----------------PLOT 6 -------------------------------
setwd("D:/Data Review Thesis/Random Forests/HELAMCFMix") 
location = 'ZphzRb_EIS_HELA_MCF_30'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)
char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#Recall RF Function with Optimal mtry & mtree
#optional
set.seed(42) 
model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                     importance=TRUE, keep.inbag=TRUE, mtry=28, ntree = 1000)

model

distmatrix = as.dist(1-model$proximity)
scaleddistmtrx = cmdscale(distmatrix, eig = TRUE, x.ret = TRUE)
mdsaxisvar = round(abs(scaleddistmtrx$eig)/sum(abs(scaleddistmtrx$eig))*100,1)
mdsvals = scaleddistmtrx$points
Samples = rownames(mdsvals)
x = mdsvals[,1]
y = mdsvals[,2]
Cell = spectranf$Type
mdsdata = data.frame(Samples, x, y, Cell)


RFplot6 = ggplot(mdsdata, aes(x = x, y = y, label = Cell))+
  geom_point(aes(color= Cell, shape = Cell))+theme_bw()+
  xlab(paste("Axis 1 - ", mdsaxisvar[1], "%", sep = ""))+
  ylab(paste("Axis 2 - ", mdsaxisvar[2], "%", sep = ""))+theme(legend.position = "none")+scale_shape(solid = FALSE)+
  ggtitle("B.1")

RFplot6
#---------- Merge RF Plots to Make Paper Figure --------------------------------------------
setwd("D:/Data Review Thesis/Plots for Paper")
tiff("rftest3Reduced.tiff", width=7480, height=4500, res=1000)
grid_arrange_shared_legend(RFplot5, RFplot1, RFplot2, RFplot6, RFplot3, RFplot4, nrow = 2, ncol = 3, position = 'bottom')
dev.off()
#-----------RAW DATA Plotting---------------------------------------------------------
#------- Plot 1 -----------------------------------------------------------

Zphznouts = readMat('D:/Data Review Thesis/HELA/60mV/Matlab Tables/HELANouts25R_Zphz_Wc_60.mat')
Zphznouts = Zphznouts$Zphznouts

Freqs = readMat('D:/Data Review Thesis/Frequency/Freq_.mat')
Freqs = Freqs$Freq
Freqs = t(Freqs)
data = data.frame(Freqs, Zphznouts)
cellnames = paste("Cell ", 1:(ncol(data)-1), sep = '')
cols = c("Freq", cellnames)
colnames(data) = cols

df <- melt(data ,  id.vars = 'Freq', variable.name = 'Sample')
rawplot1 = ggplot(df, aes(log(Freq),value)) + geom_line(aes(colour = Sample))+theme_bw() +theme(legend.position = "none")+
  xlab("log[Frequency] (Hz)") + ylab("Phase Angle (degrees)")+ggtitle(" A")

rawplot1
#------- Plot 2 -----------------------------------------------------------

Zmodnouts = readMat('D:/Data Review Thesis/HELA/60mV/Matlab Tables/HELANouts25R_Zmod_Wc_60.mat')
Zmodnouts = Zmodnouts$Zmodnouts

Freqs = readMat('D:/Data Review Thesis/Frequency/Freq_.mat')
Freqs = Freqs$Freq
Freqs = t(Freqs)
data = data.frame(log(Freqs), log(Zmodnouts))
cellnames = paste("Cell ", 1:(ncol(data)-1), sep = '')
cols = c("Freq", cellnames)
colnames(data) = cols

df <- melt(data ,  id.vars = 'Freq', variable.name = 'Sample')
rawplot2 = ggplot(df, aes(Freq,value)) + geom_line(aes(colour = Sample))+theme_bw()+ theme(legend.position = "none")+
  xlab("log[Frequency] (Hz)") + ylab(paste("log[ Z ] (", greeks("Omega"), ")", sep = ''))+ggtitle(" B")

rawplot2

#---------- Merge Raw Data Plots to Make Paper Figure --------------------------------------------
setwd("D:/Data Review Thesis/Plots for Paper")
tiff("rawtestReduced.tiff", width=7480, height=2600, res=1000)
grid.arrange(rawplot1, rawplot2, nrow = 1) 

dev.off()


#-----------PRINCIPAL COMPONENT ANALYSIS SCORES PLOTTING-----------------------------------------------------
#------------PLOT 1 ------------------------------------------------------------
setwd("D:/Data Review Thesis/HELAMCFMix/Reduced Files")
spectranf = read.csv('RedZmodRb_HELA_MCF_30.csv', header= FALSE)

colnames(spectranf)=spectranf[1,]
rownames(spectranf)=spectranf[,1]
spectranf=spectranf[-1,-1]

pca.res = pca(spectranf[,1:ncol(spectranf)], ncomp = 7, method = 'nipals', center = TRUE , scale = TRUE)
Cell = rep(c("HELA", "MCF"),times=c(nrow(spectranf)/2, nrow(spectranf)/2))
scores = data.frame(pca.res$res$cal$scores)
scplotdata = data.frame(Cell, pca.res$res$cal$scores)

ellpca <- ellipseParam(data = scores, k = 2, pcx = 1, pcy = 2)
a <- pluck(ellpca, "Ellipse", "a.99pct")
b <- pluck(ellpca, "Ellipse", "b.99pct")
pcaplot1 = ggplot(scplotdata, aes(x = Comp.1, y = Comp.2 ))+geom_point(aes(color= Cell, shape = Cell))+ 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a, b = b, angle = 0), size = 0.5, linetype = "dotted")+ggtitle("A.1")+scale_shape(solid = FALSE)+
  xlab(paste("PC 1 - ", round(pca.res$calres$expvar[1], 2), "%", sep = ""))+theme_bw()+
  ylab(paste("PC 2 - ", round(pca.res$calres$expvar[2], 2), "%", sep = ""))+ theme(legend.position = "none")

pcaplot1

#------------PLOT 2 ------------------------------------------------------------
setwd("D:/Data Review Thesis/HELAMCFMix/Reduced Files")
spectranf = read.csv('RedZmodRb_HELA_MCF_60.csv', header= FALSE)

colnames(spectranf)=spectranf[1,]
rownames(spectranf)=spectranf[,1]
spectranf=spectranf[-1,-1]

pca.res = pca(spectranf[,1:ncol(spectranf)], ncomp = 7, method = 'nipals', center = TRUE , scale = TRUE)
Cell = rep(c("HELA", "MCF"),times=c(nrow(spectranf)/2, nrow(spectranf)/2))
scores = data.frame(pca.res$res$cal$scores)
scplotdata = data.frame(Cell, pca.res$res$cal$scores)

ellpca <- ellipseParam(data = scores, k = 2, pcx = 1, pcy = 2)
a <- pluck(ellpca, "Ellipse", "a.99pct")
b <- pluck(ellpca, "Ellipse", "b.99pct")
pcaplot2 = ggplot(scplotdata, aes(x = Comp.1, y = Comp.2 ))+geom_point(aes(color= Cell, shape = Cell))+ 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a1, b = b1, angle = 0), size = 0.5, linetype = "dotted")+ggtitle("A.2")+scale_shape(solid = FALSE)+
  xlab(paste("PC 1 - ", round(pca.res$calres$expvar[1], 2), "%", sep = ""))+theme_bw()
  ylab(paste("PC 2 - ", round(pca.res$calres$expvar[2], 2), "%", sep = "")) + theme(legend.position = "none")

pcaplot2
#------------PLOT 3 ------------------------------------------------------------
setwd("D:/Data Review Thesis/HELAMCFMix/Reduced Files")
spectranf = read.csv('RedZmodRb_HELA_MCF_90.csv', header= FALSE)

colnames(spectranf)=spectranf[1,]
rownames(spectranf)=spectranf[,1]
spectranf=spectranf[-1,-1]

pca.res = pca(spectranf[,1:ncol(spectranf)], ncomp = 7, method = 'nipals', center = TRUE , scale = TRUE)
Cell = rep(c("HELA", "MCF"),times=c(nrow(spectranf)/2, nrow(spectranf)/2))
scores = data.frame(pca.res$res$cal$scores)
scplotdata = data.frame(Cell, pca.res$res$cal$scores)

ellpca <- ellipseParam(data = scores, k = 2, pcx = 1, pcy = 2)
a <- pluck(ellpca, "Ellipse", "a.99pct")
b <- pluck(ellpca, "Ellipse", "b.99pct")
pcaplot3 = ggplot(scplotdata, aes(x = Comp.1, y = Comp.2 ))+geom_point(aes(color= Cell, shape = Cell))+ 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a1, b = b1, angle = 0), size = 0.5, linetype = "dotted")+ggtitle("A.3")+scale_shape(solid = FALSE)+
  xlab(paste("PC 1 - ", round(pca.res$calres$expvar[1], 2), "%", sep = ""))+theme_bw()+
  ylab(paste("PC 2 - ", round(pca.res$calres$expvar[2], 2), "%", sep = "")) + theme(legend.position = "none")

pcaplot3
#------------PLOT 4 ------------------------------------------------------------
setwd("D:/Data Review Thesis/HELAMCFMix/Reduced Files")
spectranf = read.csv('RedZphzRb_HELA_MCF_30.csv', header= FALSE)

colnames(spectranf)=spectranf[1,]
rownames(spectranf)=spectranf[,1]
spectranf=spectranf[-1,-1]

pca.res = pca(spectranf[,1:ncol(spectranf)], ncomp = 7, method = 'nipals', center = TRUE , scale = TRUE)
Cell = rep(c("HELA", "MCF"),times=c(nrow(spectranf)/2, nrow(spectranf)/2))
scores = data.frame(pca.res$res$cal$scores)
scplotdata = data.frame(Cell, pca.res$res$cal$scores)

ellpca <- ellipseParam(data = scores, k = 2, pcx = 1, pcy = 2)
a <- pluck(ellpca, "Ellipse", "a.99pct")
b <- pluck(ellpca, "Ellipse", "b.99pct")
pcaplot4 = ggplot(scplotdata, aes(x = Comp.1, y = Comp.2 ))+geom_point(aes(color= Cell, shape = Cell))+ 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a1, b = b1, angle = 0), size = 0.5, linetype = "dotted")+ggtitle("B.1")+scale_shape(solid = FALSE)+
  xlab(paste("PC 1 - ", round(pca.res$calres$expvar[1], 2), "%", sep = ""))+theme_bw()+
  ylab(paste("PC 2 - ", round(pca.res$calres$expvar[2], 2), "%", sep = "")) + theme(legend.position = "none")

pcaplot4
#------------PLOT 5 ------------------------------------------------------------
setwd("D:/Data Review Thesis/HELAMCFMix/Reduced Files")
spectranf = read.csv('RedZphzRb_HELA_MCF_60.csv', header= FALSE)

colnames(spectranf)=spectranf[1,]
rownames(spectranf)=spectranf[,1]
spectranf=spectranf[-1,-1]

pca.res = pca(spectranf[,1:ncol(spectranf)], ncomp = 7, method = 'nipals', center = TRUE , scale = TRUE)
Cell = rep(c("HELA", "MCF"),times=c(nrow(spectranf)/2, nrow(spectranf)/2))
scores = data.frame(pca.res$res$cal$scores)
scplotdata = data.frame(Cell, pca.res$res$cal$scores)

ellpca <- ellipseParam(data = scores, k = 2, pcx = 1, pcy = 2)
a <- pluck(ellpca, "Ellipse", "a.99pct")
b <- pluck(ellpca, "Ellipse", "b.99pct")
pcaplot5 = ggplot(scplotdata, aes(x = Comp.1, y = Comp.2 ))+geom_point(aes(color= Cell, shape = Cell))+ 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a1, b = b1, angle = 0), size = 0.5, linetype = "dotted")+ggtitle("B.2")+scale_shape(solid = FALSE)+
  xlab(paste("PC 1 - ", round(pca.res$calres$expvar[1], 2), "%", sep = ""))+theme_bw()+
  ylab(paste("PC 2 - ", round(pca.res$calres$expvar[2], 2), "%", sep = "")) + theme(legend.position = "none")

pcaplot5
#------------PLOT 6 ------------------------------------------------------------
setwd("D:/Data Review Thesis/HELAMCFMix/Reduced Files")
spectranf = read.csv('RedZphzRb_HELA_MCF_90.csv', header= FALSE)

colnames(spectranf)=spectranf[1,]
rownames(spectranf)=spectranf[,1]
spectranf=spectranf[-1,-1]

pca.res = pca(spectranf[,1:ncol(spectranf)], ncomp = 7, method = 'nipals', center = TRUE , scale = TRUE)
Cell = rep(c("HELA", "MCF"),times=c(nrow(spectranf)/2, nrow(spectranf)/2))
scores = data.frame(pca.res$res$cal$scores)
scplotdata = data.frame(Cell, pca.res$res$cal$scores)

ellpca <- ellipseParam(data = scores, k = 2, pcx = 1, pcy = 2)
a <- pluck(ellpca, "Ellipse", "a.99pct")
b <- pluck(ellpca, "Ellipse", "b.99pct")
pcaplot6 = ggplot(scplotdata, aes(x = Comp.1, y = Comp.2 ))+geom_point(aes(color= Cell, shape = Cell))+ 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = a1, b = b1, angle = 0), size = 0.5, linetype = "dotted")+ggtitle("B.3")+scale_shape(solid = FALSE)+
  xlab(paste("PC 1 - ", round(pca.res$calres$expvar[1], 2), "%", sep = ""))+theme_bw()+
  ylab(paste("PC 2 - ", round(pca.res$calres$expvar[2], 2), "%", sep = "")) + theme(legend.position = "none")

pcaplot6
#---------- Merge PCA Plots to Make Paper Figure --------------------------------------------
setwd("D:/Data Review Thesis/Plots for Paper")
tiff("pcatest3Reduced.tiff", width=7480, height=4500, res=1000)
grid_arrange_shared_legend(pcaplot1, pcaplot2, pcaplot3, pcaplot4, pcaplot5, pcaplot6, nrow = 2, ncol = 3, position = 'bottom')

dev.off()