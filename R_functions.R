```R
###start primers_in_fastq
primers_in_fastq = function(primers_file, fastq_path_F, fastq_path_R) {
  require("Biostrings")
  require("phyloseq")
  require("ShortRead")
  primers = Biostrings::readDNAStringSet(filepath = primers_file)
  primers_all = c(primers,Biostrings::complement(primers),Biostrings::reverse(primers), Biostrings::reverseComplement(primers))
  
  #make results table
  primer_count = data.frame(primer = names(primers), 
                            orientation = c(rep("forward",length(primers)), rep("Comp",length(primers)), rep("Rev",length(primers)),
                                            rep("RC",length(primers))), sequence = primers_all,
                            hits_F = rep(NA,length(primers)), hits_R = rep(NA,length(primers)) )
  
  #count
  for (i in 1:length(primers_all)) {
    
    primer_count[i,"hits_F"] = sum(Biostrings::vcountPattern(toString(primers_all[i]),ShortRead::sread(ShortRead::readFastq(fastq_path_F)), fixed = F) >0 ) 
    primer_count[i,"hits_R"] = sum(Biostrings::vcountPattern(toString(primers_all[i]),ShortRead::sread(ShortRead::readFastq(fastq_path_R)), fixed = F) >0 )}
  return(primer_count)
}
###end primers_in_fastq

###start primers_in_strings
primers_in_strings = function(primers_file, sequence) {
  require("Biostrings")
  require("ShortRead")
  primers = Biostrings::readDNAStringSet(filepath = primers_file)
  primers_all = c(primers,Biostrings::complement(primers),Biostrings::reverse(primers), Biostrings::reverseComplement(primers))
#make results table
  primer_count = data.frame(primer = names(primers), 
                            orientation = c(rep("forward",length(primers)), rep("Comp",length(primers)), rep("Rev",length(primers)),
                                            rep("RC",length(primers))), sequence = primers_all,
                            hits = rep(NA,length(primers)) )
  
  #count
  for (i in 1:length(primers_all)) {
    primer_count[i,"hits"] = sum(Biostrings::vcountPattern(toString(primers_all[i]),sequence, fixed = F) >0 )}
  return(primer_count)
}


quality_graphs = function(forward,reverse){
require(dada2)
#make folder
if (!dir.exists("./qualityplots/")){ dir.create("/qualityplots/") }else{}
#plot  
pdf(file = paste0("/qualityplots/",gsub(".fq.gz","",basename(trimF[1])),"_to_",gsub(".fq.gz","",basename(trimF[length(trimF)])),"_qualityplots.pdf"))
for (i in 1:length(trimF)) {
  print(dada2::plotQualityProfile(c(forward[i],reverse[i]), aggregate = F)) 
}
dev.off()
}


filtering_tests = function(forward,reverse){
  
}


```
