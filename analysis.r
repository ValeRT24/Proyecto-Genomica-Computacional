# Utilizaremos DNAcopy para realizar la segmentación
#if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#BiocManager::install("DNAcopy")

library(DNAcopy)

# leer archivo
TCGA_B1_2 <- read.csv("TCGA_B1_2.v2.txt", sep="\t", header=TRUE)
summary(TCGA_B1_2)

CNA.object <- CNA(cbind(TCGA_B1_2$Segment_Mean),
                   + TCGA_B1_2$Start,TCGA_B1_2$Num_Probes,
                   ,sampleid="TCGA_B1_2")

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)

plot(segment.smoothed.CNA.object, plot.type="w")