```R
#libs
library("phyloseq")
library("microbiome")
library("DESeq2") # instal via BiocManager


#load data 
MagicCreatures_16S = readRDS("objects/MagicCreatures_16S_notnormalized.rds")
MagicCreatures_16S_rare = readRDS("objects/MagicCreatures_16S_rare.rds")

###

##################  Differential abundance DESeq2 ###################################################### 

#convert phyloseq to deseq2 object
sample_data(MagicCreatures_16S) # check your data to choose the "design"
MagicCreatures_16S_deseq2 = phyloseq_to_deseq2(MagicCreatures_16S, ~ habitat) # it will convert to factors
#do the testing
MagicCreatures_deseq2_diffabd = DESeq(MagicCreatures_16S_deseq2, test="Wald", fitType="parametric", minReplicatesForReplace=Inf)
#get results
MagicCreatures_deseq2_results = results(MagicCreatures_deseq2_diffabd, cooksCutoff=FALSE)
#filter out significant value(s) and make a table with taxonomy
significant_table = MagicCreatures_deseq2_results[which(MagicCreatures_deseq2_results$padj < 0.001),]
significant_table = cbind(as(significant_table, "data.frame"), as(tax_table(MagicCreatures_16S)[rownames(significant_table), ], "matrix"))
#plot a graph...
deseq2_sig_counts = otu_table(MagicCreatures_16S)[,colnames(otu_table(MagicCreatures_16S)) %in% row.names(significant_table)]
deseq2_sig_counts = cbind(deseq2_sig_counts, sample_data(MagicCreatures_16S)[,c("organism","habitat")])
deseq2_sig_counts_melt = reshape2::melt(deseq2_sig_counts, id.vars = c("organism","habitat"))
library(ggplot2)
ggplot(deseq2_sig_counts_melt, aes(x=habitat, y=value)) + 
  geom_jitter(aes(colour = organism), width = 0.2, size = 2) + facet_wrap(~variable) + theme_minimal()

#what do you think about the results ?

#TASK!# plot it with normalized data 
#TASK!# use taxonomic names instead of ASV etc..

#######################################################################################################

#libs
library("phyloseq")
library("mia") # instal via BiocManager
library("ANCOMBC") # instal via BiocManager



##################  Differential abundance ANCOM-BC ################################################### 
#convert the phyloseq to the object ancom can read using mia package
MagicCreatures_pseq = mia::makeTreeSummarizedExperimentFromPhyloseq(MagicCreatures_16S)

set.seed(353)

MagicCreatures_ancom = ancombc2(
  data = mia::makeTreeSummarizedExperimentFromPhyloseq(MagicCreatures_16S), # use mia package to convert data
  assay_name = "counts", 
  tax_level = "Genus", 
  fix_formula = "organism", 
  p_adj_method = "fdr", 
  lib_cut = 0,
  prv_cut = 0,
  group = "organism", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
 
#plot the results
library("dplyr")
library(ggplot2)
res_prim = MagicCreatures_ancom$res
res_prim$taxon = gsub(".*_","",rownames(res_prim))


df_unicorn = dplyr::select(res_prim,taxon, ends_with("unicorn")) 

df_fig_unicorn = df_unicorn %>%
  filter(diff_organismunicorn == 1) %>% 
  arrange(desc(lfc_organismunicorn)) %>%
  mutate(direct = ifelse(lfc_organismunicorn > 0, "Positive LFC", "Negative LFC"))
df_fig_unicorn$taxon = factor(df_fig_unicorn$taxon, levels = df_fig_unicorn$taxon)
df_fig_unicorn$direct = factor(df_fig_unicorn$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

fig_unicorn = df_fig_unicorn %>%
  ggplot(aes(x = taxon, y = lfc_organismunicorn, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_organismunicorn - se_organismunicorn, ymax = lfc_organismunicorn + se_organismunicorn), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of unicorn") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
fig_unicorn

#it is kind of bullshit analysis  on a small dataset but the tool is powerfull


#######################################################################################################

#boxplots - one of the ways

#source: https://github.com/joey711/phyloseq/issues/1501

#Original phyloseq object
MagicCreatures_16S_rare
#Drop taxa not present in at least 10% of samples
MagicCreatures_16S_trimmed = filter_taxa(MagicCreatures_16S_rare, function(x) sum(x>0) >= (0.1*length(x)), TRUE)
#Convert to relative abundance
MagicCreatures_16S_trimmed_RA = transform_sample_counts(MagicCreatures_16S_trimmed, function(x) round(x / sum(x), 5))
#Create boolean table (TRUE = >1% relative abundance)
MagicCreatures_16S_trimmed_RA_boolean = filter_taxa(MagicCreatures_16S_trimmed_RA, function(x) sum(x) >= 0.01, FALSE)
#Drop taxa with <1% relative abundance
MagicCreatures_16S_filtered = prune_taxa(MagicCreatures_16S_trimmed_RA_boolean, MagicCreatures_16S_trimmed)
#Plot
MagicCreatures_16S_filtered_phylum <- MagicCreatures_16S_filtered %>% tax_glom("Phylum")
phyloseq::taxa_names(MagicCreatures_16S_filtered_phylum) <- phyloseq::tax_table(MagicCreatures_16S_filtered_phylum)[, "Phylum"]

pdf("reports/boxplot_prevalence.pdf")
  ggplot(data = phyloseq::psmelt(MagicCreatures_16S_filtered_phylum), aes(x = organism, y = Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") + theme_minimal()
dev.off()
#!TASK - include more taxa by changing the <1% taxa parameter
#!TASK - Plot Genus instead of phylum
#!TASK - change "free" to "fixed" in   facet_wrap(~ OTU, scales = "free") and see whats gonna happen


###barplots easy way
pdf("reports/easy_barplots.pdf")
plot_bar(MagicCreatures_16S_rare, fill="Phylum") 
plot_bar(MagicCreatures_16S_rare, x="organism", fill="Phylum")
plot_bar(MagicCreatures_16S_rare, x="organism", fill="Phylum",facet_grid=~habitat)
dev.off()

#Question:whats wrong with those graphs ?
#how to fix it ?
```