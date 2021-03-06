# ANALYSIS Utilities/"Merge tables" (Merge two tables using row, i.e., gene names.)
# INPUT GENERIC normalized.tsv, GENERIC normalized-too.tsv OUTPUT combined.tsv
# PARAMETER include.everything [yes, no] DEFAULT no (include not matching lines into result)

# Combines two different tables using gene names
# JTT 22.10.2007

# Name of the first table
name1<-c("normalized.tsv")

# Name of the second table
name2<-c("normalized-too.tsv")

# Loads the tables
table1<-read.table(file=name1, sep="\t", header=T, row.names=1)
table2<-read.table(file=name2, sep="\t", header=T, row.names=1)

for (i in 1:ncol(table1)) {
	table1[,i] <- gsub("\t+", " ", table1[,i], perl=T)
	table1[,i] <- gsub("\n+", " ", table1[,i], perl=T)
	table1[,i] <- gsub(" +", " ", table1[,i], perl=T)
	
	table1[,i] <- gsub("\"+", ",", table1[,i], perl=T)
	table1[,i] <- gsub("\'+", ",", table1[,i], perl=T)
}

# Combines tables using row names
include=FALSE
if( include.everything == "yes" ) include=TRUE
table3<-merge(table1, table2, by.x="row.names", by.y="row.names", all=include)
row.names(table3)<-table3$Row.names
table3<-table3[-1]

# Writes out the combined table
write.table(table3, "combined.tsv", sep="\t", row.names=T, col.names=T, quote=F)
