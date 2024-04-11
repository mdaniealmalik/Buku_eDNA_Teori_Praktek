setwd("D:/DMBP_Undip/Bioinformatics_book")

library(tidyverse)
library(phyloseq)


#Read ASV table
asv_table <- read.delim("D:/DMBP_Undip/Bioinformatics_book/statistic/asv_table.txt", header=FALSE, row.names=1, comment.char="#")
colnames(asv_table) <- c("S2406", "S2408", "S2419", "S2421")
View(asv_table)

#Read Taxonomy table
taxonomy <- read.delim("D:/DMBP_Undip/Bioinformatics_book/statistic/taxonomy.tsv", row.names=1) %>% select(Taxon)
View(taxonomy)
taxonomy <- taxonomy %>% separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ';')
taxonomy <- as.matrix(taxonomy)

#Metadata
metadata <- read.csv("map.csv", header=T, sep="\t", row.names=1)

#Phyloseq Object
OTU = otu_table(asv_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
META = sample_data(metadata)

physeq1 = phyloseq(OTU, TAX, META)

#Rarefraction
plot(ranacapa::ggrare(physeq1, color = "InputFileName", se = FALSE, parallel=TRUE)) +
  geom_line(size= 1.5) +
  theme_bw() + 
  guides(color = guide_legend(title = "Sample")) +
  theme(axis.text.x = element_text(size=10, angle=0, hjust=0.5, vjust=0, family="serif"), 
        axis.text.y = element_text(size=10, family="serif"),
        legend.text = element_text(size=11, family="serif"),
        legend.title = element_text(size=12, family="serif"))

plot(ranacapa::ggrare(physeq1, color = "InputFileName", se = FALSE, parallel=TRUE)) +
  geom_line(size= 1.5) 


#Unique Taxa
get_taxa_unique(physeq1, "Phylum")
get_taxa_unique(physeq1, "Genus")

### Stacked Barplot (phylum)
## note: prune out low abundance taxa and only include Phylum that contribute more than 2% of the relative abundance of each sample
physeq_phylum <- physeq1 %>%
  tax_glom(taxrank = "Phylum") %>%                         # agglomerate at Class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%    # Transform to rel. abundance
  psmelt() %>%                                            # Melt to long format
  filter(Abundance > 0.02) %>%                            # Filter out low abundance taxa
  arrange(Phylum)

unique(physeq_phylum$Phylum)

Class_colors <- c("coral1", "chartreuse2", "cadetblue2", "burlywood2", "brown1",
                  "green3", "azure3", "darkturquoise", "darkorchid2", "darkorange2")

#between sample
ggplot(physeq_phylum, aes(x = InputFileName,  y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =Class_colors) +
  theme(axis.text.x = element_text(size=10, angle=0, hjust=0.5, vjust=0, family="serif"), 
        axis.text.y = element_text(family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"),
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16)) +
  theme(axis.title.y = element_text(size=12, vjust = +3)) +
  theme(axis.title.x = element_text(size=12, vjust = -2)) +
  ylab("Relative Abundance") +
  xlab("") +
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) + 
  theme_bw()

ggplot(physeq_phylum, aes(x = InputFileName,  y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =Class_colors) 

#between site
ggplot(physeq_phylum, aes(x = Site,  y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =Class_colors) +
  theme(axis.text.x = element_text(size=10, angle=0, hjust=0.5, vjust=0, family="serif"), 
        axis.text.y = element_text(family="serif"),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"),
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme(legend.position="right") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16)) +
  theme(axis.title.y = element_text(size=12, vjust = +3)) +
  theme(axis.title.x = element_text(size=12, vjust = -2)) +
  ylab("Relative Abundance") +
  xlab("") +
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) + 
  theme_bw()

ggplot(physeq_phylum, aes(x = Site,  y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values =Class_colors) 

##Alpha diversity
#####
p <- plot_richness(physeq1, x="InputFileName") +  
  geom_point(size=7, alpha=0.5) + xlab("")+ ylab("") + theme_bw()

estimate_richness(physeq1, measures = c("Shannon", "Simpson"))
