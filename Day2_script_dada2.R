```R
#libs
library("dada2")
library("Biostrings")
library("ShortRead")
library("magrittr")
library("phyloseq")

#set working directory if needed
getwd()
#setwd("somewhere")

#functions
source("./amplicon_functions.R")

#################data; set path to raw data
path = "./raw/" #set path to raw data

rawF <- sort(list.files(path, pattern="_1.fq.gz", full.names = TRUE, recursive = F)) # change to recursive if multi-folder
rawR <- sort(list.files(path, pattern="_2.fq.gz", full.names = TRUE, recursive = F))

#get sample names
sample.names <- sapply(strsplit(basename(rawF), "_raw_"), `[`, 1)

#################trimming

#path to cutadapt windows release
cutapath = "./cutadapt/cutadapt.exe" # set where the cutadapt.exe is
#output folder
if (!dir.exists("./trimmed_reads/")){ dir.create("/trimmed_reads/") }else{} 
#primers
primers = Biostrings::readDNAStringSet(filepath = "./primers.fasta") # path to the file with TWO primers
primers_all = c(primers,Biostrings::complement(primers),Biostrings::reverse(primers), Biostrings::reverseComplement(primers))
#trimming arguments and trimming
for (i in 1:length(rawF)) {
  cutadargs = c("-g",toString(primers_all[1]),
                "-G",toString(primers_all[2]),
                "-a",toString(primers_all[7]),
                "-A",toString(primers_all[8]),
                "-n",5,
                "-o",paste0("./trimmed_reads/", gsub("raw","trim",basename(rawF[i]))),
                "-p",paste0("./trimmed_reads/", gsub("raw","trim",basename(rawR[i]))),
                rawF[i],
                rawR[i],
                "--minimum-length", 30, # minimum length , can be longer
                "--cores",0, #all available cores
                "--report=minimal")
  system2(cutapath, args=cutadargs) } # this line runs it 

#wait ~1 min

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
#wait ~3 min 

################# learning error rates for novaseq sequencing
##make sure "magrittr" is loaded
require("magrittr")
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
save.image(file = "rawr_amplicons.RData")
##

################# dereplication
derepF <- dada2::derepFastq(filtF, verbose=T)
names(derepF) <- sample.names
derepR <- dada2::derepFastq(filtR, verbose=T)
names(derepR) <- sample.names

#wait ~1 min

################# sample interference
dadaF <- dada2::dada(derepF, err=errF, multithread=T, pool=F) # pool sensitivity vs. time
dadaR <- dada2::dada(derepR, err=errR, multithread=T, pool=F)
#wait ~2 min

################# merging reads
mergers <- dada2::mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = 12, verbose=T)

################# Construct ASV table
seqtab <- dada2::makeSequenceTable(mergers)

################# Chimeras removal
seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

##checkpoint2
save.image(file = "rawr_amplicons.RData")
##

################# assign taxonomy
# find the database online and download !
taxonomy <- dada2::assignTaxonomy(seqtab.nochim, "./SILVA138/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#wait ~2-3 min

#saveRDS(taxonomy, file = "MamestraP3_taxonomy.rds")
#taxa = dada2::addSpecies(taxa, "../databases/SILVA138_1/silva_species_assignment_v138.1.fa.gz") #need more memory, use server
#taxonomy_sp = readRDS(file = "taxonomy/MamestraP3_taxonomy_sp.rds")

##checkpoint3
save.image(file = "rawr_amplicons.RData")
##
#################################################################################################################################





####################################################### export phyloseq #######################################################

##metadata
#libs to load
library("phyloseq")
library("Biostrings")

metadata = read.csv(file = "./amplicons_metadata.csv", sep = ";", row.names = 1)

MagicCreatures_16S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                              sample_data(metadata), 
                              tax_table(taxonomy))

# change sequences to ASVn
ASVsDNA = Biostrings::DNAStringSet(taxa_names(MagicCreatures_16S))
names(ASVsDNA) = phyloseq::taxa_names(MagicCreatures_16S)
MagicCreatures_16S = phyloseq::merge_phyloseq(MagicCreatures_16S, ASVsDNA)
taxa_names(MagicCreatures_16S) = paste0("ASV", seq(ntaxa(MagicCreatures_16S)))

#see whats there using View(command(dataset))

saveRDS(MagicCreatures_16S, file= "./MagicCreatures_16S_notnormalized.rds")

###############################################################################################################################






####################################################### reports and checks #######################################################
#make folder
if (!dir.exists(paste0("/reports/"))){ dir.create(paste0("/reports/")) }else{}

################# check for primers 
source("./amplicon_functions.R")
primers_hits_raw = primers_in_fastq("primers.fasta", rawF, rawR)# takes "primers.fasta" file with TWO primers
primers_hits_trim = primers_in_fastq("primers.fasta", trimF, trimR)
primers_hits_filtered = primers_in_fastq("primers.fasta", filtF, filtR)
primers_hits_final = primers_in_strings("primers.fasta", refseq(MagicCreatures_16S))

primers_summary = cbind(primers_hits_raw,primers_hits_trim[,4:5],primers_hits_filtered[,4:5],primers_hits_final[,4])
colnames(primers_summary)[4:10] = c("hits_F_raw","hits_R_raw","hits_F_trim","hits_R_trim","hits_F_filt","hits_R_filt","hits_final")

write.csv2(primers_summary, paste0("./reports/","primers.csv") )

#wait it takes a LONG time

################# errors
pdf(file = paste0("./reports/","dada2_errors.pdf"))
dada2::plotErrors(errF, nominalQ=TRUE, err_in = T)
dada2::plotErrors(errR, nominalQ=TRUE, err_in = T)
dev.off()


################# filtering, dechimering
filtered_report = as.data.frame(filtered_report)
filtered_report$survived = round(filtered_report$reads.out / filtered_report$reads.in * 100,2)
filtered_report$seq = rowSums(seqtab)
filtered_report$seq.nochim = rowSums(seqtab.nochim)
filtered_report$seq_pr_rem = round(filtered_report$seq.nochim / filtered_report$seq * 100 , 2)

write.csv2(filtered_report, paste0("./reports/","filtered_report.csv") )

################# length distribution
pdf(file = paste0("./reports/","length_distribution.pdf"))
hist(nchar(phyloseq::refseq(MagicCreatures_16S)),breaks = 50 )
dev.off()

################# reads quality - filtered
pdf(file = paste0("./reports/","qualityplots_filtered.pdf")) # saving as pdf
for (i in 1:length(filtF)) { print(dada2::plotQualityProfile(c(filtF[i],filtR[i]), aggregate = F, n = 5000)) }
dev.off() # end saving as pdf

```
