###
library(foreach)

# Break up the number of contigs (i.e. 27922) into 50 contig chunks
seq(from=1, to=27922, by= 50) -> vec.num
seq(from=1, to=88155, by= 50) -> vec.num

dat.win=foreach(i=vec.num, .combine = "rbind")%do%{
  data.frame(start=i, end=i+49)
  }

write.table(dat.win, file = "dat.win.partitions.txt", append = FALSE, quote = FALSE, sep = "\t",
eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")