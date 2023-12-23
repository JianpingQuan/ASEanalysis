## 3D scatter

library(scatterplot3d)
pc <- read.table("C:\\Users\\qjp\\Desktop\\BS\\PCA\\data\\SamplesPCA_prcomp_x.txt", header = T, row.names = 1)

Tissue <- c()
tis <- substr(rownames(pc),7,7)

for(i in 1:nrow(pc)){
  if(tis[i]=="B"){
    Tissue[i] <- "Br"
  }else if(tis[i]=="L"){
    Tissue[i] <- "Li"
  }else if(tis[i]=="M"){
    Tissue[i] <- "Mu"
  }else if(tis[i]=="P"){
    Tissue[i] <- "Pl"
  }
}

pc$Tissue <- Tissue


shapes = c(15, 16, 17, 18)
shapes <- shapes[factor(pc$Tissue)]
colors <- c("#3B4992", "#EE2200", "#008B45", "#631779")
colors <- colors[factor(pc$Tissue)]

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

s3d <- scatterplot3d(dat[,1:3],
              main=NULL,
              pch = shapes,
              color=colors,
              grid=FALSE,
              box=FALSE,
              xlab = "PC1 (52.9%)",
              ylab = "PC2 (7.1%)",
              zlab = "PC3 (5.9%)")

# 3. Add grids
addgrids3d(dat[,1:3], grid = c("xy", "xz", "yz"))

legend("top", legend = levels(factor(pc$Tissue)), col =  c("#3B4992", "#EE2200", "#008B45", "#631779"), pch = c(15, 16, 17, 18),inset =-0.4 , xpd = TRUE, horiz = TRUE)