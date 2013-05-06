# Utilities for adding "chr" to chromosome names int GTF and FASTA


addChrToGtf <- function(file.name) {
	
	base.name <- sub("\\.[[:alnum:]]*$", "", file.name)
	new.name <-
	system(paste("awk '{ if ($1 ~ /^[0-9,X,Y,MT,Z,W]/) {print "chr"$0} else {print $0} }'", file.name, ">", new.name))
	
}


#isGZipFile <- function(file.name) {
#	
#	# get file type with the unix file command
#	file.type = system(paste("file -Lb --mime", file.name), intern=TRUE)
#	
#	if (!is.na(pmatch("application/x-gzip", file.type))) {
#		return(TRUE);
#	} else { 
#		return(FALSE);
#	}
#}
