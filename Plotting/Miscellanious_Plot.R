
# To plot the SSM96 mutational profile of all the sample for that particular tumour type
# and to plot the signature contributions of all the samples in that particular tumour type (boxplot).

source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_profile_SSM96_normalized.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/normalizeMotifs.R")
require(Biostrings)


mut_counts <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_","Lung-SCC",".txt",sep=""),header=TRUE,check.names=FALSE) 
mut_counts <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_","Lung-SCC",".txt",sep=""),header=TRUE,check.names=FALSE) 
mut_norm <- sweep(mut_counts, 2, colSums(mut_counts), `/`)

# Load the correction factor from cds2wgs and exon2wgs
load("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Correction_Factors.RData")
# In order to use the correction function we need to renamed the rownames as in the SomSign package
subs <-  c('CA','CG','CT','TA','TC','TG')
context <- c('A.A','A.C','A.G','A.T','C.A','C.C','C.G','C.T',
             'G.A','G.C','G.G','G.T','T.A','T.C','T.G','T.T')
ssm96_somsign <- paste(rep(subs,each=16), rep(context,times=4), sep=" ")
ssm96_std_rownames <- rownames(mut_counts)

# Apply the correction
# Second, correction is applied to Exon or CDS region -> WGS. No manipulation to COSMIC
rownames(mut_counts) <- ssm96_somsign
mut_counts_corr <- normalizeMotifs(mut_counts, corr_factorCDS2WGS)

# Coming back to the standard nomenclature
rownames(mut_counts) <- ssm96_std_rownames
rownames(mut_counts_corr) <- ssm96_std_rownames

mut_corr_norm <- sweep(mut_counts_corr, 2, colSums(mut_counts_corr), `/`)

pdf("Lung-SCC_WGS.pdf",height=18)
plot_profile_SSM96_normalized(mut_norm, ymax = 0.2, condensed = FALSE)
dev.off()


# Signature Contribution

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Lung-SCC",".RData",sep=""))
source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_summary_SignaturesContribution.R")

for(i in 1: length(fits_m1))
{
  pdf(paste0("SignContr_",names(fits_m1)[i],".pdf",sep=""))
  p <- plot_summary_SignaturesContribution(fits_m1[[i]]$contribution, whichtumour="Lung-SCC", norm=TRUE,data=names(fits_m1)[i]) 
  plot(p)
  dev.off()
}
