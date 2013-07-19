# TOOL dexseq-exoncount.R: "Count aligned reads per exons for DEXSeq" (Given mapped reads in a BAM file, this tool counts the reads that fall into each non-overlapping exonic part using the script dexseq-count.py. In order to use the output in DEXSeq, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT exon-counts.tsv
# PARAMETER paired: "Does the alignment file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER organism: "Organism" TYPE [hg19: "Human (hg19.72)", mm10: "Mouse (mm10.68)", rn4: "Rat (rn4.68)"] DEFAULT hg19 (Which organism is your data from.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1: "chr1", 1: "1"] DEFAULT chr1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)

# 18.9.2012 TH and EK 
# 16.7.2013 EK, BAM sorting changed
# 19.7.2013 AMS, added parameter for chr1/1 

# convert bam to sam, sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
samtools.sort <- ifelse(paired == "yes", paste(samtools.binary, "sort -on alignment.bam sorted-by-name"), "cat alignment.bam")

samtools.view <- paste(samtools.binary, "view alignment.bam")

# gtf
annotation.file <- ""
if (organism == "hg19") {
	if (chr == 1){
		annotation.file <- "Homo_sapiens.GRCh37.72.DEXSeq.gtf"
	}else {
		annotation.file <- "Homo_sapiens.GRCh37.72.chr.DEXSeq.gtf"
	}		
}
if (organism == "mm9") {
	if (chr == 1){
		annotation.file <- "Mus_musculus.NCBIM37.62.DEXSeq.gtf"
	}else {
		annotation.file <- "Mus_musculus.NCBIM37.62.chr.DEXSeq.gtf"
	}
}
if (organism == "mm10") {
	if (chr == 1){
		annotation.file <- "Mus_musculus.GRCm38.68.DEXSeq.gtf"
	}else{
		annotation.file <- "Mus_musculus.GRCm38.68.chr.DEXSeq.gtf"
	}
}
if (organism == "rn4") {
	if (chr == 1){
		annotation.file <- "Rattus_norvegicus.RGSC3.4.68.DEXSeq.gtf"
	}else{
		annotation.file <- "Rattus_norvegicus.RGSC3.4.68.chr.DEXSeq.gtf"
	}
}
gtf <- file.path(chipster.tools.path, "genomes", "gtf", annotation.file)


# exoncount
dexseq.binary <- file.path(chipster.tools.path, "dexseq-exoncounts", "dexseq_count.py")
paired.end.data <- ifelse(paired == "yes", paste("-p yes"), "")
dexseq.command <- paste("python", dexseq.binary, paired.end.data, gtf,"- exon-counts.tsv")

# run 
command <- paste(samtools.view, "|", dexseq.command)
system(command)


# bring in file to R environment for formating
file <- c("exon-counts.tsv")
dat <- read.table(file, header=F, sep="\t")
names(dat) <- c("id", "count")

# write result table to output
write.table(dat, file="exon-counts.tsv", col.names=T, quote=F, sep="\t", row.names=F)


# EOF


