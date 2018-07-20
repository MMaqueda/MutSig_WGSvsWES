
# This function is from Mutational Patterns package but it's been impossible to install!!!
source("/project/devel/PCAWG/mmaqueda/Rscripts/cos_sim_matrix.R") 
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")

# Plot Delta-Cos Similarity between (Corrected - No corrected) versus #SSMs
# Two cases: CDS and Exon

# Computations were already done during descriptive analysis. Let's retrieve that data-

load("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/Comparison_WGS_WES.RData")

# Prepare data
toplot <- data.frame(Tumour = sapply(colnames(mut96_mat_allCDS), function(x) unlist(strsplit(x, "[.]"))[2]),
                     SSMs_CDS = colSums(mut96_mat_allCDS),
                     SSMs_Exon = colSums(mut96_mat_allExon),
                     Delta_CosSim_CDS = cos_sim_WGS_WEScds_corr$cos_sim - cos_sim_WGS_WEScds$cos_sim,
                     Delta_CosSim_Exon = cos_sim_WGS_WESexon_corr$cos_sim - cos_sim_WGS_WESexon$cos_sim)

############
# Let's plot
cols_tumors <- pcawg.colour.palette(x = levels(toplot$Tumour),scheme = "tumour.subtype")
names(cols_tumors) <- levels(toplot$Tumour)

pdf("CDSvsCosSimXTumor.pdf",width=12,height=10)
ggplot(toplot, aes(y=SSMs_CDS, x=Delta_CosSim_CDS,color=Tumour)) +
  facet_wrap( ~ Tumour,scale="free") +
  geom_point() + 
  scale_color_manual(values=cols_tumors) +
  scale_y_continuous(trans = "log10") + 
  geom_vline(xintercept=0, linetype='dashed') +
  geom_hline(yintercept=100, linetype = 'dashed') +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=8),
                     strip.text = element_text(size=7),
                     legend.position="none") +
  #annotate(geom="text", label="100 SSMs", x=-0.25, y=100, vjust=-1) +
  #annotate("label", x = -0.15, y = 1, label = "Less similar \n after correction", size=5, colour = "red") +
  #annotate("label", x = 0.15, y = 1, label = "More similar \n after correction", size=5, colour = "green") +
  ggtitle("Number of SSMs vs Delta Cos Sim (Corr-Nocorr) - CDS") 
dev.off()


