# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))
# Set relative path from root
rel_path_from_root <- find_root_file("results", "tables", criterion = has_file("README.md"))

# List files in this folder to make sure you're in the right spot.
list.files(rel_path_from_root) 
# Set working directory as path from root
setwd(rel_path_from_root)


###
library(foreach)

# Break up the number of contigs (i.e. 27922) into 50 contig chunks
seq(from=1, to=19014, by= 50) -> vec.num

genome.chunks=foreach(i=vec.num, .combine = "rbind")%do%{
  data.frame(start=i, end=i+49)
  }

write.table(genome.chunks, file = "genome.chunks.txt", append = FALSE, quote = FALSE, sep = "\t",
eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")