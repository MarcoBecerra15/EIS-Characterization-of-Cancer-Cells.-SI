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
#--Call overall set-----------------------------------------------------------------
setwd("D:/Data Review Thesis/Random Forests/HELAMDAMix") 
location = 'ZmodRb_EIS_HELA_MDA_150'
spectranf = read.csv(paste(location, ".csv", sep = ""), header= FALSE)
colnames(spectranf) = spectranf[1,]
rownames(spectranf) = spectranf[,1]
spectranf = spectranf[-1,-1]
spectranf$Type = as.factor(spectranf$Type)

#__Convert Z from character to numeric________

char_columns <- sapply(spectranf, is.character)             
numspectra <- spectranf                              
numspectra[ , char_columns] <- as.data.frame(   
  apply(numspectra[ , char_columns], 2, as.numeric))

#----Enter optimal parameters---------------------------------------------------
opttrees = as.numeric(readline(prompt = "Enter optimal number of trees of tunned model \n"))
optmtry = as.numeric(readline(prompt = "Enter optimal variables per split of tunned model \n"))

#----LOOCV LOOP --------------------------------------------------------------

type = factor(levels = levels(numspectra$Type))
PerformanceUlt=data.frame() #Only execute once when rerunning the code consequently
#~~Make submodel and Prediction Sets~~~~~~~~~
for (n in 1:nrow(numspectra)) {
  subspect = numspectra[-n,]
  rownames(subspect) = c(1:nrow(subspect))
  predictset = numspectra[n,]
  rownames(predictset) = 1
  
  #~Generate RF submodel~~~~~~~~~~~~~~~
  set.seed(42) #optional
  submodel = randomForest(Type ~ ., data = subspect, proximity=TRUE, 
                       importance=TRUE, keep.inbag=TRUE, ntree = opttrees, mtry = optmtry)
  
  #Use Model to Predict 
  pred = predict(submodel, predictset, type="response", predict.all=TRUE, proximity = TRUE)
  type[n] = pred$predicted$aggregate

} #End of Loop~~~~~~~~~~~~~~~~~~~~~~~~~~

#------Compile results into Dataframe-----------------------------------------------
 predicted = data.frame( Cell = 1:length(type), Prediction = type)
 predicted$Real_Type = numspectra$Type
 
confmatrix = confusionMatrix(predicted$Real, predicted$Prediction)
Results = confmatrix$overall

Component = substr(location, 1, 4) #Zmod or Zphz
Voltage = substr(location, nchar(location)-2, nchar(location))

Performance = data.frame(Component, Voltage, Accuracy=Results[1]*100, Conf_95Int_Low= Results[3]*100, Conf_95Int_High=Results[4]*100, P_value=Results[6])
rownames(Performance) = 1

PerformanceUlt = rbind(Performance, PerformanceUlt)
Performance=data.frame() #reset value



#----Continue when done passing results over the loop------------------------------------------


if (grepl("HELA_MCF", location) == 'TRUE') {
  Mix = "HELA_MCF"
  
  } else if (grepl("HELA_MDA", location) == 'TRUE') {
  Mix = "HELA_MDA"
  
  } else if (grepl("MCF_MDA", location) == 'TRUE') {
  Mix = "MCF_MDA"
}


#---SAVE RESULTS TO CSV FILE (1 SHEET PER MIX !!!)--------------------------------

write.xlsx(PerformanceUlt, file="D:/Data Review Thesis/Random Forests/RF_LOOCV.xlsx", row.names = FALSE, col.names = TRUE, sheetName = Mix, append = TRUE)  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
