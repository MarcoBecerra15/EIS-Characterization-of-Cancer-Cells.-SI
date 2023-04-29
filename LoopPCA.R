library(utils)
library(base)
library(stats)
library(mdatools)
library(prospectr)
library(ellipse)
library(plotrix)
library(ggplot2)
library(tiger)

datafldr = "D:/Data Review Tesis/"
cells <- readline(prompt= "Choose Mix for PCA: \na) HELA-MDA \nb) MCF-MDA \nc) HELA-MCF \n")

#Define variables based on mix choice------------------------------------------
if (cells == "a") { 
  cells = "HELAMDAMix"
  loadtitle = "HELA-MDA"
  c1 = "HELA"
  c2 = "MDA"
} else if (cells == "b") {
  cells = "MCFMDAMix"
  loadtitle = "MCF-MDA"
  c1 = "MCF"
  c2 = "MDA"
} else if  (cells == "c") {
  cells = "HELAMCFMix"
  loadtitle = "HELA-MCF"
  c1 = "HELA"
  c2 = "MCF"
}

#Select and set Impedance components to iterate on-----------------------------
zcomp <- c("Zmod", "Zphz")
askcomp <- readline(prompt ="Enter impedance component for PCA (Zmod, Zphz or both seperated by a space) \n")
if (askcomp == "both") {
  askcomp = zcomp
} 
itcomps = intersect(askcomp,zcomp)

type <- readline(promp ="Choose data type: \nRnb=no baseline with outliers, Rnob=no baseline nor outliers, Rb=with baseline no outliers, Rbo=with baseline & outliers \n")

#Specify voltages to iterate on---------------------------------------------------------
voltages <- c("30", "60", "90", "120", "150")
askvolts <- readline(prompt ="Enter voltage(s) for PCA (30, 60, 90, 120, 150)mV. \nSeperate with single space. For all, type all. \n")

if (askvolts == "all") {
  askvolts = voltages
} else if (nchar(askvolts) > 2) {
  askvolts = unlist(strsplit(askvolts, " "))
}
itvolts = intersect(askvolts,voltages)

#Set Data Folder
setwd(paste(datafldr, cells, sep=""))

#LOOP FOR Z COMPONENT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~---
for (zcount in 1:length(itcomps)) {
  #LOOP FOR VOLTAGE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~--
  for (vcount in 1:length(itvolts)) {
    
    spectra = read.csv(paste(itcomps[zcount], type, "_EIS_", c1, "_", c2, "_", itvolts[vcount], ".csv", sep=""), header = FALSE)
    
    colnames(spectra)=spectra[1,]
    rownames(spectra)=spectra[,1]
    spectra=spectra[-1,-1]
    q = nrow(spectra) / 2
    
    #------SCORES PLOT----------------------------------------------------------
    pca.res = pca(spectra[,1:56], ncomp = 7, method = 'nipals', center = TRUE , scale = TRUE)
    loadings = pca.res$loadings
    
    scorestitle = paste(c1, "_", c2, "_", itcomps[zcount], itvolts[vcount], "_", type, " Scores.pdf", sep="") 
    pdf(scorestitle)
    p = plotScores(pca.res$res$cal, xlim = c(-20,20), ylim = c(-13,13), show.labels = FALSE, col=c(rep('red', q),rep('blue', q)))
    plotHotellingEllipse(p, conf.lim = 0.99, col = "black")
    legend(-20,13,legend = c( c1, c2), pch = c(16,16), cex=0.6, col = c('red','blue'))
    dev.off()
    
    #------LOADINGS PLOTS-------------------------------------------------------
    pca.res = pca(spectra[,1:56], ncomp = 7, method = 'nipals')
    
    loadings = pca.res$loadings #loading not centered
    pcompsfilename = paste("PrinComps_", itcomps[zcount], type, "_", c1, "_", c2, "_", itvolts[vcount], ".csv", sep="")
    write.csv(loadings, pcompsfilename)
    
    loadings1title = paste(c1, "_", c2, "_", itcomps[zcount], itvolts[vcount], "_", type, " Loadings C1.pdf", sep="")
    loadings2title = paste(c1, "_", c2, "_", itcomps[zcount], itvolts[vcount], "_", type, " Loadings C2.pdf", sep="")
    plotc1 = paste(c1, "-", c2, " ", itvolts[vcount], "mV ", itcomps[zcount], type, " Loadings C1")
    plotc2 = paste(c1, "-", c2, " ", itvolts[vcount], "mV ", itcomps[zcount], type, " Loadings C2")
    
    pdf(loadings1title)
    plot(rownames(loadings), loadings[,1], xlab = 'Frequency (Hz)', ylab = 'Loading', , cex.lab = 1, main = plotc1, type = 'b')
    dev.off()
    
    pdf(loadings2title)
    plot(rownames(loadings), loadings[,2], xlab = 'Frequency (Hz)', ylab = 'Loading', , cex.lab = 1, main = plotc2, type = 'b')
    dev.off()
  }
  }


































