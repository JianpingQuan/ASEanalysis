##
cd /mnt/home/quanjian/resources/ASElist

grep -f AGE.list <(awk '$3=="gene"' /mnt/home/quanjian/resources/Sus_scrofa.Sscrofa11.1.98.gtf|cut -f1,4,5,9)|cut -d ';' -f1|sed 's/gene_id "//g'|sed 's/"//g'|sort -k1,1n -k2,2n > AGE.bed

grep -f POE2.list <(awk '$3=="gene"' /mnt/home/quanjian/resources/Sus_scrofa.Sscrofa11.1.98.gtf|cut -f1,4,5,9)|cut -d ';' -f1|sed 's/gene_id "//g'|sed 's/"//g'|sort -k1,1n -k2,2n > POE2.bed

python POEgapLength.py|cut -d ' ' -f8,9|sed 's/://g'|tr ' ' '\t' > POE.gapLength2.txt

python AGEgapLength.py|cut -d ' ' -f8,9|sed 's/://g'|tr ' ' '\t' > AGE.gapLength.txt


##Set all distances less than 0 to 0

##
setwd("F:\\Manuscript\\RNA\\Ori\\figure\\POE_AGE_density")
library(dplyr)

poe <- read.table("Shuf179.gapLength.txt", header = F)
colnames(poe) <- c("Chr","GapLength")

for(i in 1:nrow(poe)){
  if(poe[i,2]<0){
    poe[i,2]=abs(poe[i,2])
  }
}

poe$GapLength <- poe$GapLength+0.1

hist(log10(poe$GapLength),breaks = seq(0,9,1))

a <- hist(log10(poe$GapLength),breaks = seq(0,9,1))
b <- data.frame(Breaks=a$breaks[-length(a$breaks)], Counts=a$counts, CumCount=cumsum(a$counts))
b$Pro <- b$CumCount/179*100

##
age <- read.table("AGE.gapLength.txt", header = F)
colnames(age) <- c("Chr","GapLength")

for(i in 1:nrow(age)){
  if(age[i,2]<0){
    age[i,2]=abs(age[i,2])
  }
}

age$GapLength <- age$GapLength+0.1

hist(log10(age$GapLength),breaks = seq(0,9,1))

c <- hist(log10(age$GapLength),breaks = seq(0,9,1))
d <- data.frame(Breaks=c$breaks[-length(c$breaks)], Counts=c$counts,CumCount=cumsum(c$counts))
d$Pro <- d$CumCount/394*100



##211 genes and 394 genes were randomly selected from the autosomes for the above analysis
awk '$3=="gene"' /mnt/home/quanjian/resources/Sus_scrofa.Sscrofa11.1.98.gtf|cut -f1,4,5,9|awk '$1 ~ /^[1-9]/'|cut -d ';' -f1|sed 's/gene_id "//g'|sed 's/"//g'|shuf -n 394|sort -k1,1n -k2,2n > shuf394.bed

python ShufgapLength.py|cut -d ' ' -f8,9|sed 's/://g'|tr ' ' '\t' > Shuf394.gapLength.txt

setwd("F:\\Manuscript\\RNA\\Ori\\figure\\POE_AGE_density")
library(dplyr)

shuf <- read.table("Shuf394.gapLength.txt", header = F)
colnames(shuf) <- c("Chr","GapLength")

for(i in 1:nrow(shuf)){
  if(shuf[i,2]<0){
    shuf[i,2]=abs(shuf[i,2])
  }
}

shuf$GapLength <- shuf$GapLength+0.1

hist(log10(shuf$GapLength), breaks = seq(0, 9, by = 1))

a <- hist(log10(shuf$GapLength),breaks = seq(0, 9, by = 1))
b <- data.frame(Breaks=a$breaks[-length(a$breaks)], Counts=a$counts, CumCount=cumsum(a$counts))
b$Pro <- b$CumCount/394*100



##