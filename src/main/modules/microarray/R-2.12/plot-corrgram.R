# TOOL plot-corrgram.R: Correlogram (Plots the correlations between samples in a graph.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS 
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC 
# OUTPUT corrgram.pdf: corrgram.pdf 
# PARAMETER image.width: image.width TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the plotted network image)
# PARAMETER image.height: image.height TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the plotted network image)


# Correlogram
# JTT 18.10.2007

# Parameter settings (default) for testing purposes
#image.width<-c(600)
#image.height<-c(600)

# Renaming variables
w<-image.width
h<-image.height

# Loads libraries
library(corrgram)

# Loads the normalized data
file<-c("normalized.tsv")
dat<-read.table(file, header=T, sep="\t", row.names=1)

# Separates expression values and flags
calls<-dat[,grep("flag", names(dat))]
dat2<-dat[,grep("chip", names(dat))]

# Loads phenodata
phenodata<-read.table("phenodata.tsv", header=T, sep="\t")

# Generating new labels
colnames(dat2)<-gsub(" ", "", phenodata$description)

# Plotting
pdf(file="corrgram.pdf", width=w/72, height=h/72)
corrgram(x=dat2, lower.panel=panel.shade, upper.panel=panel.pts, pch=".")
dev.off()


