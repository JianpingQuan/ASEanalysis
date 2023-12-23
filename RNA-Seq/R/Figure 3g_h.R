setwd("F:\\Manuscript\\RNA\\Ori\\figure\\AGE_POE_TissuesVenn")

library(VennDiagram)

##POE
Brain <- read.table("POE_Brain_bias2.txt", header=F)[,1]
Liver <- read.table("POE_Liver_bias2.txt", header=F)[,1]
Muscle <-read.table("POE_Muscle_bias2.txt", header=F)[,1]
Placenta <-read.table("POE_Placenta_bias2.txt", header=F)[,1]


venn.plot <- venn.diagram(
  x = list(Brain = Brain, Muscle = Muscle, Liver = Liver,Placenta=Placenta),
  filename = "Tissues_Venn_poe.tiff",
  imagetype = "tiff",
  col = "black",
  lty = 5,
  lwd = 1,
  fill = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF"),
  alpha = 0.5,
  #4 group
  label.col = c("black", "white", "black", "white","white", "white","white", "white", "black", "white","white", "white", "white", "black", "white"), cex = 2.0,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("black","black","black","black"),
  cat.cex = 1.8,
  cat.fontfamily = "serif")


##AGE
Brain <- read.table("AGE_Brain_bias2.txt", header=F)[,1]
Liver <- read.table("AGE_Liver_bias2.txt", header=F)[,1]
Muscle <-read.table("AGE_Muscle_bias2.txt", header=F)[,1]
Placenta <-read.table("AGE_Placenta_bias2.txt", header=F)[,1]

venn.plot <- venn.diagram(
  x = list(Brain = Brain, Muscle = Muscle, Liver = Liver,Placenta=Placenta),
  filename = "Tissues_Venn_age.tiff",
  imagetype = "tiff",
  col = "black",
  lty = 5,
  lwd = 1,
  fill = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF"),
  alpha = 0.5,
  #4 group
  label.col = c("black", "white", "black", "white","white", "white","white", "white", "black", "white","white", "white", "white", "black", "white"), cex = 2.0,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("black","black","black","black"),
  cat.cex = 1.8,
  cat.fontfamily = "serif")
