```R
#libs
library("dada2")
library("Biostrings")
library("ShortRead")
library("magrittr")
library("phyloseq")
library("vegan")

#functions
source("./amplicon_functions.R")

####################################################### export to phyloseq ######################################################

##metadata
#libs to load
library("phyloseq")
library("Biostrings")

metadata = read.csv(file = "./amplicons_metadata.csv", sep = ";", row.names = 1)

MagicCreatures_16S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                              sample_data(metadata), 
                              tax_table(taxonomy))

# change sequences to ASVs to make the table clear
ASVsDNA = Biostrings::DNAStringSet(taxa_names(MagicCreatures_16S))
names(ASVsDNA) = phyloseq::taxa_names(MagicCreatures_16S)
MagicCreatures_16S = phyloseq::merge_phyloseq(MagicCreatures_16S, ASVsDNA)
taxa_names(MagicCreatures_16S) = paste0("ASV", seq(ntaxa(MagicCreatures_16S)))

#save the results
saveRDS(MagicCreatures_16S, file= "./objects/MagicCreatures_16S_notnormalized.rds")

###############################################################################################################################


####################################################### normalize #############################################################
##rarefy
library(vegan)
asv_table = as.matrix(as.data.frame(phyloseq::otu_table(MagicCreatures_16S)))
set.seed(456) # make it "random"

pdf("reports/rarecurve.pdf")rare_curve = vegan::rarecurve(asv_table, step = 100 , cex=0.5)
dev.off()

MagicCreatures_16S_rare = phyloseq::rarefy_even_depth(MagicCreatures_16S, rngseed=456, sample.size=0.9*min(sample_sums(MagicCreatures_16S)), replace=F)

saveRDS(MagicCreatures_16S_rare, file= "objects/MagicCreatures_16S_rare.rds")
###############################################################################################################################

####################################################### reports and checks #######################################################

################# length distribution
pdf(file = "reports/length_distribution.pdf")
hist(nchar(phyloseq::refseq(MagicCreatures_16S)),breaks = 100 )
dev.off()


####################################################### visualizations #######################################################
library("phyloseq")
#diversity
plot_richness(MagicCreatures_16S_rare, measures=c("Observed","Chao1", "Shannon"))
#TASK!plot diversity per habitat including ACE metric, save both plots in one pdf file

#ordination
Magic_ordination <- ordinate(MagicCreatures_16S_rare, method = "NMDS", distance = "bray")
plot_ordination(MagicCreatures_16S_rare, Magic_ordination, type="samples", 
                color="organism", shape="habitat") +
geom_point(size=5) + ggtitle("Magic Creatures NMDS, Bray-Curtis") + theme_bw()

#TASK!Do Principal coordinates analysis using Bray-Curtis distances and
# Redundancy analysis using Jaccard distances
#save both plots in one pdf file



##########################
#playing with ampvis2
install.packages("remotes")
remotes::install_github("kasperskytte/ampvis2")
library("ampvis2")
#convert to ampvis format
MagicCreatures_16S_rare_ampvis = phyloseq_to_ampvis2(MagicCreatures_16S_rare)
#diversity with ampvis
amp_alpha_diversity(MagicCreatures_16S_rare_ampvis, rarefy = NULL)
#heatmaps with ampvis2
amp_heatmap(data = MagicCreatures_16S_rare_ampvis, tax_aggregate = "Family",
            group_by = "organism")

#TASK! plot a heatmap with facet using habitat on Genus level and different color scale

#octave plot
amp_octave( data = MagicCreatures_16S_rare_ampvis, tax_aggregate = "Genus", 
            group_by = "organism", scales = "fixed" )
#venn diagram
amp_venn(data = MagicCreatures_16S_rare_ampvis, group_by= "organism", 
                  normalise = F)

##################################################################################################################################
```
