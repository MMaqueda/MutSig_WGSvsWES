
# This script if or analyzing the effect of the 3mer correction relative contr in the mutational profiles

setwd("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/WGSvsWES_kmerCorrection/Analysis_kmer_CorrectionEffect/")
load("k3mer_CorrectionEffect.RData")

# Compute the rel weight per k-3mer seq in the different regions. No chrom differentiation.

k3mer_WGS <- colSums(rel_3mer_WGS)
k3mer_CDS <- colSums(rel_3mer_CDS)
k3mer_Exon <- colSums(rel_3mer_Exon)

#  And now we can compute the correction factors (Exon -> WGS, CDS -> WGS):
corr_factorCDS <- k3mer_WGS / k3mer_CDS
corr_factorExon <- k3mer_WGS / k3mer_Exon

# We could compute the correction factor per chromosome

corr_factorCDS_chrom <- rel_3mer_WGS / rel_3mer_CDS
corr_factorExon_chrom <- rel_3mer_WGS / rel_3mer_Exon

# Plot the results:
# First prepare the data to be plotted

toplot <- melt(as.data.frame(corr_factorCDS_chrom))
colnames(toplot) <- c("k3mer","Corr_cds2wg")

corrExon <- melt(as.data.frame(corr_factorExon_chrom))
toplot <- cbind(toplot, Corr_exon2wg = corrExon$value)

toplot <- melt(toplot)
colnames(toplot) <- c("k3mer","Region","Correction")

# We wanna include the abs correction values applied (without considering chromosome)

abs_corr <- data.frame(k3mer = names(corr_factorCDS), cds2wg = corr_factorCDS, exon2wg = corr_factorExon)
abs_corr <- melt(abs_corr)
colnames(abs_corr) <- c("k3mer","Region","Correction")

# And now the boxplot

boxplot <- ggplot(toplot, aes(y=Correction, x=k3mer, fill=Region)) + 
  geom_boxplot() + 
  geom_point(data = abs_corr,
             aes(x=k3mer, y=Correction, shape=Region)) +
  scale_shape_manual(values=c(3, 4)) +
  ggtitle("Correction factor values - Exon or CDS to WG - considering chromosome dispersion") +
  theme_bw() + theme(axis.text.x = element_text(angle=90,size=11,vjust=0.5),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) 

pdf("Comparison_correctionXchrom.pdf", width=20)
boxplot
dev.off()