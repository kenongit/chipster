# Utilities for adding "chr" to chromosome names int GTF and FASTA


addChrToGtf <- function(oldfile.name, newfile.name) {
	# Add chr to chromosome name (assumed to be any name that is plain number or X, Y, Z, W or MT)
	# Note: awk selected for robustness in detecting whitespaces.
	#system(paste("awk '{ if ($1 ~ /^[0-9,X,Y,Z,W]/) {print \"chr\"$0} else {print $0} }'", oldfile.name, ">", newfile.name))
	system(paste("awk '{if ($1 ~ /^[0-9][0-9]*$/ || $1 == \"X\" || $1 == \"Y\" || $1 == \"Z\" || $1 == \"W\" || $1 == \"MT\") {print \"chr\"$0}  else {print $0}}'", oldfile.name, ">", newfile.name))
	
}
