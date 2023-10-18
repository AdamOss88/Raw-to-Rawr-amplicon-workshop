```R
#################libs
library("dada2")

#################data
path = "./raw/" #set path to raw data

rawF <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE, recursive = F)) # change to recursive if multi-folder
rawR <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE, recursive = F))


################# reads quality
#make folder
if (!dir.exists(paste0("/reports/"))){ dir.create(paste0("/reports/")) }else{}

pdf(file = paste0("./reports/","qualityplots_rawdata.pdf")) # saving as pdf
for (i in 1:length(rawF)) { print(dada2::plotQualityProfile(c(rawF[i],rawR[i]), aggregate = F, n = 5000)) }
dev.off() # end saving as pdf

#wait ~5-10 min


```