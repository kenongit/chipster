# TOOL infoseq_summary.R: "Sequence file summary" (Tool to calculate basic properties of a sequence file.)
# INPUT input.txt: "Query sequences" TYPE GENERIC
# OUTPUT summary_file.txt

# K.M 28.10.2013


# pb settings
emboss.path <- file.path(chipster.tools.path, "emboss" ,"bin")
infoseq.binary <- file.path(chipster.tools.path, "blast", "/ncbi-blast-2.2.28+", "bin", "infoseq_summary.bash")
command.full <- paste(infoseq.binary, emboss.path, "input.txt > summary_file.txt 2>&1" )
system(command.full)
