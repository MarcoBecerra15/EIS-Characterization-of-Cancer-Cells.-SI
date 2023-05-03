library(utils)
library(base)
library(stats)
library(mdatools)
library(prospectr)
library(ellipse)
library(plotrix)
library(ggplot2)
library(tiger)

setwd("C:\\Users\\Eugenio J. Caraballo\\Desktop\\R Analysis\\Electric Tests for Chips\\CSV Files")


#----------------------------------------------------- Scores Plot -----------------------------------------------------


spectra = read.csv('Reduced Chips Impedance Magnitude (90 mV) {Impedance Magnitude as a Function of Frequency}.csv')

rownames(spectra) = spectra[,2]
colnames(spectra)[3:7] = spectra[1,3:7]
spectra = spectra[-1,-1]
spectra = spectra[,-1]

pca.res = pca(spectra[,2:5], ncomp = 7, method = 'nipals', center = TRUE, scale = TRUE)

p = plotScores(pca.res$res$cal, xlim = c(-10,10), ylim = c(-2,2), show.labels = FALSE, col=c(rep('red',1),rep('blue',1),rep('darkgreen',1),rep('darkgrey',1),rep('brown',1),rep('aquamarine',1),rep('chartreuse',1),rep('antiquewhite',1),rep('chocolate',1),rep('black',1),rep('cyan',1),rep('darkgoldenrod1',1)))
plotHotellingEllipse(p, conf.lim = 0.99, col = "black")
legend(-10,2,legend = c('Chip 1','Chip 2','Chip 3','Chip 4','Chip 5','Chip 6','Chip 7','Chip 8','Chip 9','Chip 10','Chip 11','Chip 12'), pch = c(16,16), col = c('red','blue','darkgreen','darkgrey','brown','aquamarine','chartreuse','antiquewhite','chocolate','black','cyan','darkgoldenrod1'))


spectra = read.csv('Reduced Chips Phase Angle (90 mV) {Phase Angle as a Function of Frequency}.csv')

rownames(spectra) = spectra[,2]
colnames(spectra)[3:46] = spectra[1,3:46]
spectra = spectra[-1,-1]
spectra = spectra[,-1]

pca.res = pca(spectra[,2:44], ncomp = 7, method = 'nipals', center = TRUE, scale = TRUE)

p = plotScores(pca.res$res$cal, xlim = c(-25,25), ylim = c(-20,20), show.labels = FALSE, col=c(rep('red',1),rep('blue',1),rep('darkgreen',1),rep('darkgrey',1),rep('brown',1),rep('aquamarine',1),rep('chartreuse',1),rep('antiquewhite',1),rep('chocolate',1),rep('black',1),rep('cyan',1),rep('darkgoldenrod1',1)))
plotHotellingEllipse(p, conf.lim = 0.99, col = "black")
legend(-25,20,legend = c('Chip 1','Chip 2','Chip 3','Chip 4','Chip 5','Chip 6','Chip 7','Chip 8','Chip 9','Chip 10','Chip 11','Chip 12'), pch = c(16,16), col = c('red','blue','darkgreen','darkgrey','brown','aquamarine','chartreuse','antiquewhite','chocolate','black','cyan','darkgoldenrod1'))


#---------------------------------------------------- Loadings Plot ----------------------------------------------------


spectra = read.csv('Chips Impedance Magnitude (90 mV) {Impedance Magnitude as a Function of Frequency}.csv')

rownames(spectra) = spectra[,2]
colnames(spectra)[3:58] = spectra[1,3:58]
spectra = spectra[-1,-1]
spectra = spectra[,-1]

pca.res = pca(spectra[,2:56], ncomp = 7, method = 'nipals')

loadings = pca.res$loadings

plot(rownames(loadings), loadings[,2], xlab = 'Frequency (Hz)', ylab = 'Loading', , cex.lab = 1, main = 'Loadings', type = 'b')


spectra = read.csv('Chips Phase Angle (90 mV) {Phase Angle as a Function of Frequency}.csv')

rownames(spectra) = spectra[,2]
colnames(spectra)[3:58] = spectra[1,3:58]
spectra = spectra[-1,-1]
spectra = spectra[,-1]

pca.res = pca(spectra[,2:56], ncomp = 7, method = 'nipals')

loadings = pca.res$loadings

plot(rownames(loadings), loadings[,2], xlab = 'Frequency (Hz)', ylab = 'Loading', , cex.lab = 1, main = 'Loadings', type = 'b')
