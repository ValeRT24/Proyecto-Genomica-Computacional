# Utilizaremos DNAcopy para realizar la segmentación
#if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#BiocManager::install("DNAcopy")

library(DNAcopy)

data("coriell")

CNA.object <- CNA(cbind(coriell$Coriell.05296),
                   + coriell$Chromosome,coriell$Position,
                   ,sampleid="c05296")

smoothed.CNA.object <- smooth.CNA(CNA.object)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)


plot(segment.smoothed.CNA.object, plot.type="w")

plot(segment.smoothed.CNA.object, plot.type="s")

plot(segment.smoothed.CNA.object, plot.type="p")