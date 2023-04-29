#                                                    DiazLab Microfluidics Laboratory
#                                                               Malcom D?az Garc?a
#                                                                 16/April/2023
#                                                    University of Puerto Rico at Mayaguez 
#                                                                 Any Cell
# Main Functions:
#   1) Make Publication Looking Table 

library(rempsyc)
library(xlsx)
library(gtsummary)
library(gtExtras)
library(gt)
library(flextable)
library(officer)
setwd("C:/Users/MALCOMDIAZGARCIA/OneDrive - University of Puerto Rico/Documents/Investigacion Bioelectra Charact Cells/Paper")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test2 = read.xlsx('RF_LOOCV.xlsx', 8)
test2[,2:ncol(test2)] = round(test2[,2:ncol(test2)], 1) 
for (i in 1:nrow(test2) ) {
  
  if (test2$Component[i] == "Zphz") {
    test2$Component[i] = "Phase Angle"
    
  } else if (test2$Component[i] == "Zmod") {
    test2$Component[i] = "Magnitude"
  }
    
}
colnames(test2)[1] = "Impedance Component"


xd = flextable(test2)
xd

xd = delete_part(xd, part = "header")
xd = add_header_row(xd, values = c("Impedance Component", "Voltage (mV)", rep(c("HeLa vs MCF", "HeLa vs MDA", "MCF vs MDA"), each = 2)))
xd = add_header_row(xd, values = c("Impedance Component", "Voltage (mV)", rep("Accuracy", each = 6)))
xd = add_header_row(xd, values = c("Impedance Component", "Voltage (mV)", rep(c("LOOCV", "OOBER"), times = 3)), top = FALSE)
xd = merge_h(xd, part = "header")
xd = merge_v(xd, part = "header")
xd = align(xd, align = "center", part = "all")
xd = merge_v(xd, j = "Impedance Component")
xd = hline(xd, part = "header")

xd = border(xd, i = 6, border.top = fp_border(width = 2))
xd = border(xd, part = "header", border.bottom = fp_border(width = 2))
xd = border(xd, part = "body", i = 10, border.bottom = fp_border(width = 2))

xd = font(xd, fontname = "Times New Roman", part = "all")
xd = fontsize(xd, size = 12, part = "all")
xd = line_spacing(xd, space = 2, part = "all")
xd

xd = border(xd, part = "header", border.top = fp_border(width = 2) )
#border(xd, i=10, j=1, border.top = fp_border(width = 2))
xd

xd = add_body_row(xd, top = FALSE, values  = rep("", each = 8))
xd = border(xd, i=10, border.bottom = fp_border(width = 0))
xd
