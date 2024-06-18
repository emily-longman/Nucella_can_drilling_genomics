# PCAngsd

# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

dir(find_root_file("results", criterion = has_file("README.md")))
# Set relative path from root
rel_path_from_root <- find_root_file("results", "output", "pcangsd", criterion = has_file("README.md"))

# List files in this folder to make sure you're in the right spot.
list.files(rel_path_from_root)
# Set working directory as path from root
setwd(rel_path_from_root)