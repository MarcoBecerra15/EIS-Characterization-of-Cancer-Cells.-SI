#                                                    DiazLab Microfluidics Laboratory
#                                                               Malcom Díaz García
#                                                                 05/May/2022
#                                                    University of Puerto Rico at Mayaguez 
#                                                                 Any Cell
# Main Functions:
#   1) Perform Random Forest Model 
#   2) Evaluate its Basic Parameters
#   3) Plot Variable Importance Based On Accuracy
#   4) Generate MDS Plot to Visualize Clusters


#----------------Pt. 1: Specify Libraries, Call Data & Format it------------------------


library (cowplot)
library (randomForest)
library (gghighlight)
library(RColorBrewer)
library(xlsx)

setwd("C:/Users/MALCOMDIAZGARCIA/OneDrive - University of Puerto Rico/Documents/Investigacion Bioelectra Charact Cells/cd4-jURKAT") 
location = 'Jurkat & CD4 - Phase Angle (45 mV) {Phase Angle as a Function of Frequency}'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)

#_____Convert Z from character to numeric________

char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#-----------------Pt. 2: Initial Call Random Forest Function----------------------------------------------

  set.seed(42) #optional
  model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                       importance=TRUE, keep.inbag=TRUE, ntree = 1000)
  
  model



#-----------------Part. 3: Plot Trees vs Errors and save to PC (To optimize ntree)--------------------

oob.error.data <- data.frame( 
  Trees=rep(1:nrow(model$err.rate), times=3),
  Type=rep(c("OOB", "Jurkat", "CD4"), each=nrow(model$err.rate)),
  Error=c(model$err.rate[,"OOB"], 
          model$err.rate[,"Jurkat"], 
          model$err.rate[,"CD4"]))

ggplot(data=oob.error.data, aes(x=Trees, y=Error)) +
  geom_line(aes(color=Type))+ ggtitle("Error Rates After Computing N-Amount of Trees")
#ggsave("Error_Rates_1000_Trees.pdf") #Optional




#-----------------Pt. 4: Check Optimal Number of Variables (To optimize mtry)-------------------------

oob.values <- vector(length=28) #length and p max must be equal
#Modify ntree Based on Pt. 3 Beforehand !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for (p in 1:28) {
  set.seed(42)
  temp.model = randomForest(Type ~ ., data=numspectra, mtry=p, ntree=750,proximity=TRUE, importance=TRUE, keep.inbag=TRUE)
  oob.values[p] = temp.model$err.rate[nrow(temp.model$err.rate),1]
}
optvarit = which(oob.values == min(oob.values))
optmtry = max(optvarit) #Personal Choice: Go w/ highest numb. if >1 have same min. OOBER
looptOOB = min(oob.values)
oob.values


#Recall RF Function with Optimal mtry & mtree
#optional
set.seed(42) 
model = randomForest(Type ~ ., data = numspectra, proximity=TRUE, 
                     importance=TRUE, keep.inbag=TRUE, mtry=27, ntree = 750)

model

#-----------------Pt. 5: Variable Importance Plot----------------------------------------------------

impvarmodel=as.data.frame(varImpPlot(model))
impvarmodel = data.frame(rownames(impvarmodel),impvarmodel)
rownames(impvarmodel) = 1:nrow(impvarmodel)
colnames(impvarmodel) = c("Frequency", "MeanDecreaseAccuracy", "GiniImpurity")
str(impvarmodel) 

impvarmodel = impvarmodel[order(impvarmodel$MeanDecreaseAccuracy, decreasing = TRUE),]
impvarmodel$Frequency = factor(impvarmodel$Frequency, levels = impvarmodel$Frequency)
top=head(impvarmodel, n=10)
ggplot(impvarmodel, aes(x=MeanDecreaseAccuracy, y=Frequency))+
  geom_bar(stat='identity', fill = "tomato3")+
  labs(y="Frequency", x="Importance")+
  gghighlight(impvarmodel$MeanDecreaseAccuracy == head(impvarmodel$MeanDecreaseAccuracy, n=10))+
  ggtitle("Variable Importance by Mean Decrease in Accuracy")



#----------------Part. 6: MDS Plot----------------------------------------------------------

distmatrix = as.dist(1-model$proximity)
scaleddistmtrx = cmdscale(distmatrix, eig = TRUE, x.ret = TRUE)
mdsaxisvar = round(abs(scaleddistmtrx$eig)/sum(abs(scaleddistmtrx$eig))*100,1)
mdsvals = scaleddistmtrx$points
Samples = rownames(mdsvals)
x = mdsvals[,1]
y = mdsvals[,2]
Celltype = spectranf$Type
mdsdata = data.frame(Samples, x, y, Celltype)

ggplot(mdsdata, aes(x = x, y = y, label = Celltype))+
  geom_text(aes(color= Celltype))+theme_bw()+
  xlab(paste("MDS1 - ", mdsaxisvar[1], "%", sep = ""))+
  ylab(paste("MDS2 - ", mdsaxisvar[2], "%", sep = ""))+
  ggtitle("MDS Plot Using Random Forest Proximities")


scaleddistmtrx$eig
#---------------- {OPTIONAL} Part. 7: Save Array for MDS Plot----------------------------------------------------------

mdsdata$AxisVariation = mdsaxisvar 
write.xlsx(mdsdata, file="C:/Users/MALCOMDIAZGARCIA/OneDrive - University of Puerto Rico/Documents/Investigacion Bioelectra Charact Cells/cd4-jURKAT/RFMDS_PlotJDdata.xlsx", row.names = FALSE, col.names = TRUE, sheetName = location, append = TRUE)
