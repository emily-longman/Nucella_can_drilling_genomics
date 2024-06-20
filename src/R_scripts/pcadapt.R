# pcadapt and OutFLANK

# Install packages:
install.packages("devtools")
install.packages(c("BiocManager","vcfR","pcadapt"))
BiocManager::install("qvalue") # When prompted, type "a" - this will take some time
library(devtools)
install_github("whitlock/OutFLANK")

# Load packages:
library(pcadapt)
library(qvalue)
