library(httr)

mitomap <- "https://mitomap.org/mitomaster/websrvc.cgi"
r = POST(mitomap, body = list(file = upload_file("mySequences.fasta"), fileType = "snvlist", output = "detail"), encode = "multipart")
df <- read.table(text=content(r, "text"), sep = "\t", header = T)
