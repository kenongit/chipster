# TOOL dexseq-exoncount.R: "Count aligned reads per exons for DEXSeq" (Given mapped reads in a BAM file, this tool counts the reads that fall into each non-overlapping exonic part using the script dexseq-count.py. In order to use the output in DEXSeq, you need to select all samples and run the tool \"Utilities - Define NGS experiment\".)
# INPUT alignment.bam: "BAM alignment file" TYPE GENERIC
# OUTPUT exon-counts.tsv
# PARAMETER paired: "Does the alignment file contain paired-end data" TYPE [yes, no] DEFAULT no (Does the alignment data contain paired end or single end reads?)
# PARAMETER organism: "Annotation GTF" TYPE [Homo_sapiens.GRCh37.72.DEXSeq.gtf: "Human (hg19.72\)", Mus_musculus.NCBIM37.62.DEXSeq.gtf: "Mouse (mm9.62\)", Mus_musculus.GRCm38.68.DEXSeq.gtf: "Mouse (mm10.68\)", Rattus_norvegicus.RGSC3.4.68.DEXSeq.gtf: "Rat (rn4.68\)"] DEFAULT Homo_sapiens.GRCh37.72.DEXSeq.gtf (You can use own GTF file or one of those provided on the server.)
# PARAMETER chr: "Chromosome names in my BAM file look like" TYPE [chr1: "chr1", 1: "1"] DEFAULT chr1 (Chromosome names must match in the BAM file and in the reference annotation. Check your BAM and choose accordingly.)

# 18.9.2012 TH and EK 
# 16.7.2013 EK, BAM sorting changed
# 19.7.2013 AMS, added parameter for chr1/1 

# convert bam to sam, sort bam if the data is paired-end
samtools.binary <- file.path(chipster.tools.path, "samtools", "samtools")
samtools.sort <- ifelse(paired == "yes", paste(samtools.binary, "sort -on alignment.bam sorted-by-name"), "cat alignment.bam")

samtools.view <- paste(samtools.binary, "view alignment.bam")

# GTF: If chr version is required, we create a local copy with chr ont he fly, otherwise we use the existing file 
annotation.file <- file.path(chipster.tools.path, "genomes", "gtf", organism)
if (chr == 1){
	gtf <- annotation.file
}else{ 
	source(file.path(chipster.common.path, "AddChr.R"))
	addChrToGtf(annotation.file, organism)
	gtf <- organism
}

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


