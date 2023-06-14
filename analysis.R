#ANALYSIS OF FUNGAL SEQUENCING DATA (VIA R STUDIO)

#PREPROCESSING 

# Importing data 
library(phyloseq)
library(ggplot2)
library(ape)

#import OTU table (from .csv file)
otu<-read.csv(file.choose(), header=TRUE, check.names=FALSE)

rownames(otu)<-otu$OTUID
otu<-otu[-c(1)]

taxonomy <- read.csv("taxonomy.csv", sep = ",", row.names = 1)
taxonomy <- as.data.frame(taxonomy)
taxmat = matrix(nrow = nrow(taxonomy), ncol = 7)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

library(stringr), che
split<-str_split_fixed(taxonomy$Taxon, ";", 7)
rownames(split)<-rownames(taxonomy)
colnames(split)<-colnames(taxmat)
taxonomy2<-as.matrix(split)


#import as csv – changed rownames to sample IDs 
metadata4<-read.csv(file.choose(metatdata.tsv), header=TRUE)
rownames(metadata4)<-metadata4$sample.id
metadata4<-metadata4[-c(1)]
META<-sample_data(metadata4)
phy_tree <- read_tree("tree.nwk")



OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy2)
META <- sample_data(metadata4)

# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

# Same sample names
sample_names(OTU)
sample_names(META)

# Finally merge!
ps <- phyloseq(OTU, TAX, META, phy_tree)
ps


#ALPHA DIVERSITY 

pruned = prune_taxa(taxa_sums(ps) > 0, ps)
#subsetted to Clinical Activity vs Remission 
lod<-subset_samples(pruned, Remission.State!= "NA")
sample_data(lod)$Remission.State<-as.character(sample_data(lod)$Remission.State)
#plot alpha diversity
plot_richness(lod, x="Remission.State", measures=c("Observed", "Shannon"), color="Remission.State")+geom_point(size=5, alpha=0.7) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_discrete(c("0","1")) + scale_color_manual(values = c("navyblue", "firebrick")) + geom_boxplot(width=0.5)+ geom_jitter(height = 0, width = 0.25, size=3)

results = estimate_richness(lod, measures = c(“Shannon”, “Observed”)
d = sample_data(lod)
remission = results[d[,'Remission.State'] == '0',]
activity = results[d[,'Remission.State'] == '1',]
wilcox.test(remission,activity)
summary(remission)
summary(activity)


#BETA DIVERSITY 
pruned = prune_samples(sample_sums(ps) > 0, ps)
sample_sums(pruned)
lod<-subset_samples(pruned, Remission_cat2!= "NA")
sample_data(lod)$Remission_cat2<-as.character(sample_data(lod)$Remission_cat2)

For weighted unifrac + NMDS
wunifrac_dist = phyloseq::distance(lod, method="unifrac", weighted=T)
NMDS.lod <- ordinate(lod, "NMDS", wunifrac_dist, weakties=FALSE)
plot_ordination(lod, NMDS.lod, color="Remission_cat2") + theme_bw() +theme(aspect.ratio=1)+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#with ellipses
ordplot<-plot_ordination(lod, NMDS.lod, color="Remission_cat2") + theme_bw() +theme(aspect.ratio=1)+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

library("ggplot2")
ordplot + 
+     stat_ellipse(type = "norm", linetype = 2) +
+     theme_bw()


#For unweighted unifrac + PCOA
#Unifrac
wunifrac_dist = phyloseq::distance(lod, method="unifrac", weighted=F)
ordination = ordinate(lod, method="PCoA", distance=wunifrac_dist)
plot_ordination(lod, ordination, color="Remission.State") + theme_bw() +theme(aspect.ratio=1)+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#Plotted Beta Diversity
plot_ordination(lod, ordination, color="Remission.State") + theme_bw() +theme(aspect.ratio=1)+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

results2<-adonis2(wunifrac_dist ~ sample_data(lod.rarefied)$Remission.State)
summary(results2)
write_tsv(results2, "betadivadonis.tsv")

#FRACTIONAL PREVALENCE OF PHYLUM AND GENUS

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")
library(metagMisc)
ps3<-phyloseq_filter_taxa_tot_fraction(ps, frac = 0.01)
ps.genus = tax_glom(ps3, taxrank="Genus", NArm=FALSE)

prevdf2 <- apply(X = otu_table(ps.genus), MARGIN = ifelse(taxa_are_rows(ps.genus), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
prevdf2 <- data.frame(Prevalence = prevdf2, TotalAbundance = taxa_sums(ps.genus), tax_table(ps.genus))
plyr::ddply(prevdf2, "Genus", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
ggplot(prevdf2, aes(TotalAbundance, Prevalence / nsamples(ps.genus),color=Genus)) +geom_point(size = 7, alpha = 0.7) +geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +scale_x_log10() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]")

## then looking at relative abundance box plots for these organisms 
## go back to ps then remove rare organisms, then transform to RA
ps3<-phyloseq_filter_taxa_tot_fraction(ps, frac = 0.01)
ps3 = tax_glom(ps3, taxrank="Genus", NArm=FALSE)

GPr  = transform_sample_counts(ps3, function(x) x / sum(x) )
ps.phylum = tax_glom(ps, taxrank="Genus", NArm=FALSE)
taxa_names(ps.genus)<-tax_table(ps.genus)[,"Genus"]
p<-psmelt(ps.genus)
ggplot(data = p, aes(x = OTU, y = Abundance)) +
+     geom_boxplot(outlier.shape = NA) +
+     geom_jitter(aes(color = OTU), height = 0, width = .2) +
+     labs(x = "", y = "Abundance\n") +
+     facet_wrap(~ OTU, scales = "free")+theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#DIFFERENTIAL ABUNDANCE ANALYSIS

#adding 1 pseudocount to all OTU counts and remaking physeq object 
Otu2<-otu+1
OTU2 <- otu_table(Otu2, taxa_are_rows = TRUE)
Ps.differential<-phyloseq(OTU2,TAX,META,phy_tree)

#removed rare organisms
ps3<-phyloseq_filter_taxa_tot_fraction(ps.differential, frac = 0.01)
ps.genus = tax_glom(ps3, taxrank="Genus", NArm=FALSE)

head(sample_data(ps)$Remission.State, 25)
ps = subset_samples(ps, Remission.State != "NA")
head(sample_data(ps)$Remission.State, 25)
sample_data(ps)$Remission.State<-as.character(sample_data(ps)$Remission.State) 
# do not need this as.character step if you want to model “2” as a higher fold for more severe disease (2 vs 1 vs 0)

diagdds = phyloseq_to_deseq2(ps, ~ Remission.State)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

write_tsv(sigtab, “deseqresults.tsv”)


## then looking at relative abundance box plots for these organisms 
## go back to ps then remove rare organisms, then transform to RA
Ps4<-phyloseq_filter_taxa_tot_fraction(ps.differential, frac = 0.01)
Ps4 = tax_glom(ps4, taxrank="Genus", NArm=FALSE)
total = median(sample_sums(ps4))
standf = function(x, t=total) round(t * (x / sum(x)))
gps = transform_sample_counts(ps4, standf)





#Adjusting for confounders in relative abundance

ps.na = subset_samples(ps, Remission.State != "NA")
ps.combo = subset_samples(ps.na, biologics!=”NA”)
ps.combo2 = tax_glom(ps.combo, taxrank="Genus", NArm=FALSE)
diagdds = phyloseq_to_deseq2(ps.combo2, design=~biologics+Genderscore+Agescore+Remission.State) 

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.combo2)[rownames(sigtab), ], "matrix"))
head(sigtab)

x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#ASV ANALYSIS

library("Biostrings")
rep.seqs <- Biostrings::readDNAStringSet(file.choose(), format = "fasta")

install.packages("remotes")
remotes::install_github("jfq3/RDPutils")

library(RDPutils)
ps <- phyloseq(OTU, TAX, META, phy_tree, rep.seqs)
Ps4<-phyloseq_filter_taxa_tot_fraction(ps, frac = 0.01)
GP.candida <- subset_taxa(Ps4, Genus=="g__Candida")
writeXStringSet(refseq(GP.candida), filepath=<file>, append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

#CORRELATION ANALYSIS

candata2<-Remissiontop30candida
candata2<-as.data.frame(candata2)
rownames(candata2)<-candata2$Sample_ID
candata2<-candata2[-c(1)]
cor_mat <- cor(candata2)
corrplot(cor_mat, method = 'circle', is.corr=FALSE, col=colorRampPalette(c("blue","white","red"))(6), type = 'upper', diag = FALSE)

# to determine p values
cor_test_mat <- corr.test(candata2)$p
cor_test_mat


M<-cor(candata2)
corrplot(M, method="circle", tl.cex=0.5)
