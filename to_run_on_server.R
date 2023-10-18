```R
#libs
library("dada2")
library("Biostrings")
library("ShortRead")
library("magrittr")
library("phyloseq")

#functions
source("./amplicon_functions.R")

#################data; set path to raw data
path = "./rawdata/" #set path to raw data

rawF <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE, recursive = F)) # change to recursive if multi-folder
rawR <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE, recursive = F))

#get sample names
sample.names <- sapply(strsplit(basename(rawF), "_raw_"), `[`, 1)

#################trimming

#primers
primers = Biostrings::readDNAStringSet(filepath = "./primers.fasta") # path to the file with TWO primers

path_trim = "./trimmed_reads/"
trimF <- sort(list.files(path_trim, pattern="trim_1.fq.gz", full.names = TRUE, recursive = F))
trimR <- sort(list.files(path_trim, pattern="trim_2.fq.gz", full.names = TRUE, recursive = F))

#################


#################filtering; provide path to trimmed reads - path trim

#paths to filtering output
filtF <- file.path("filtered", paste0(sample.names, "_filt_1.fastq.gz")) ## here you can add a tested parameter
filtR <- file.path("filtered", paste0(sample.names, "_filt_2.fastq.gz"))

filtered_report <- dada2::filterAndTrim(fwd = trimF, filt = filtF, 
                                        rev = trimR, filt.rev =  filtR, 
                     minLen = 30, # usually 20
                     rm.phix = T,
                     maxN=0,
                     truncQ=2, # optimize this parameter
                     maxEE=c(2,2), # optimize this parameter
                     compress=T, 
                     multithread=T,
                     verbose=T)

################# learning error rates for novaseq sequencing
##make sure "magrittr" is loaded

set.seed(33) # makes the error learning consistent

errF <- dada2::learnErrors(filtF,
                    nbases = 1e8,
                    errorEstimationFunction = loessErrfun_mod4,# skip for NOT novaseq
                    randomize = T,
                    MAX_CONSIST = 15,
                    multithread =T,
                    verbose = TRUE)

errR <- dada2::learnErrors(filtR,
                    nbases = 1e8,
                    errorEstimationFunction = loessErrfun_mod4, # skip for NOT novaseq
                    randomize = T,
                    MAX_CONSIST = 15,
                    multithread =F,
                    verbose = TRUE)

#wait ~2 min each


##checkpoint1
save.image(file = "dada2.RData")
##

################# dereplication
derepF <- dada2::derepFastq(filtF, verbose=T)
names(derepF) <- sample.names
derepR <- dada2::derepFastq(filtR, verbose=T)
names(derepR) <- sample.names


################# sample interference
dadaF <- dada2::dada(derepF, err=errF, multithread=T, pool=F) # pool sensitivity vs. time
dadaR <- dada2::dada(derepR, err=errR, multithread=T, pool=F)

################# merging reads
mergers <- dada2::mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = 12, verbose=T)

################# Construct ASV table
seqtab <- dada2::makeSequenceTable(mergers)

################# Chimeras removal
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

##checkpoint2
save.image(file = "dada2.RData")
##


```
