# This R script will create a bamfile lists based on binary drilling ability

argv <- commandArgs(T)
bamlist <- argv[1]
metadata <- argv[2]

# Read in bamlist
bamlist<-read.table(bamlist, header=F)

# Read in metadata
metadata<-read.csv(metadata, header=T)

# Extract bamlist for Drilled
Drilled.bamlist <- bamlist[which(metadata$Drilled.Binary == "Drilled"),]

# Extract bamlist for Not.Drilled
Not.Drilled.bamlist <- bamlist[which(metadata$Drilled.Binary == "Not.Drilled"),]


# Write bam lists 
write.table(Drilled.bamlist, "Nucella_bam_Drilled.list", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(Not.Drilled.bamlist, "Nucella_bam_Not.Drilled.list", col.names = FALSE, row.names = FALSE, quote = FALSE)