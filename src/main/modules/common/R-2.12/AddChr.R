# Utilities for adding "chr" to chromosome names int GTF and FASTA


addChrToGtf <- function(oldfile.name, newfile.name) {
	# Add chr to chromosomes (assumed to be any name that is plain number or X,Y,Z or W)
	# Note: awk selected for robustness in detecting whitespaces.
	system(paste("awk '{ if ($1 ~ /^[0-9,X,Y,Z,W]/) {print \"chr\"$0} else {print $0} }'", oldfile.name, ">", newfile.name))
	
}


addChrToFasta <- function(oldfile.name, newfile.name) {
	# Add chr to chromosomes (assumed to be any name that is plain number or X,Y,Z or W)
	# Note: sed selected for speed.
	system(paste("sed /^.[0-9,X,Y,Z,W]/s/\">\"/\">chr\"/  <", oldfile.name, ">", newfile.name))
}
	
