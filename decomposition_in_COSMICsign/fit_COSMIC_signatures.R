
require(MutationalPatterns)
require(deconstructSigs)
require(SomaticSignatures)
require(ggplot2)
require(reshape)

# Identify the COSMIC signatures contribution in a SSM96 mutational profile

# In a prior step, the 30 COSMIC signatures were already obtained and prepared (matrix 96 x30)
# The needed object is called 'cancer_signatures'
load("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/COSMICsign.RData")

# Select the specific mut matrix per tumour from the complete SSM96 matrix
source("/project/devel/PCAWG/mmaqueda/Rscripts/select_tumour.R")

mut_mat_WGS <- select_tumour(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/mut96_mat_allWGS.txt",
     resultsPath = "/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/WGS/Lung-AdenoCA/",
     whichtumour = c("Lung-AdenoCA"))

# In case  of Exon and CDS we already have the specific mut_matrix

mut_mat_Exon <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_Lung-AdenoCA.txt",header=TRUE,check.names=FALSE)
mut_mat_CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_Lung-AdenoCA.txt",header=TRUE,check.names=FALSE)

# Method ONE: fit_to_signatures function from MutationalPatterns ----------
# This method asks for counts (not freq)

# First, no correction is applied in Exon or CDS data
fit_COSMIC_WGS <- fit_to_signatures(mut_mat_WGS, cancer_signatures)
fit_COSMIC_Exon <- fit_to_signatures(mut_mat_Exon, cancer_signatures)
fit_COSMIC_CDS <- fit_to_signatures(mut_mat_CDS, cancer_signatures)

# Second, correction is applied to COSMIC signatures to Exon or CDS region
# We need the 3mer freq values (k3mer_Exon/WGS/CDS) for computing correction factors

#load("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/WGSvsWES_kmerCorrection/Analysis_kmer_CorrectionEffect/k3mer_CorrectionEffect.RData")
corr_factorCDS2WGS <- k3mer_WGS / k3mer_CDS
corr_factorExon2WGS <- k3mer_WGS / k3mer_Exon
corr_factorWGS2CDS <- k3mer_CDS / k3mer_WGS
corr_factorWGS2Exon <- k3mer_Exon / k3mer_WGS

# Correct the signatures thru normalizeMotifs (SomaticSignatures)
# In order to use the correction function we need to renamed the rownames as in the SomSign package

subs <-  c('CA','CG','CT','TA','TC','TG')
context <- c('A.A','A.C','A.G','A.T','C.A','C.C','C.G','C.T',
             'G.A','G.C','G.G','G.T','T.A','T.C','T.G','T.T')

ssm96_somsign <- paste(rep(subs,each=16), rep(context,times=4), sep=" ")
ssm96_std_rownames <- rownames(mut96_mat_allCDS_norm)

rownames(cancer_signatures) <- ssm96_somsign

cancer_signa_corr2Exon <- normalizeMotifs(cancer_signatures, corr_factorWGS2Exon)
cancer_signa_corr2CDS <- normalizeMotifs(cancer_signatures, corr_factorWGS2CDS)

# Return to our standard rownames and also the new calculated
rownames(cancer_signatures) <- ssm96_std_rownames
rownames(cancer_signa_corr2Exon) <- ssm96_std_rownames
rownames(cancer_signa_corr2CDS) <- ssm96_std_rownames

# Fit the corrected signatures
fit_COSMICcorr_Exon <- fit_to_signatures(mut_mat_Exon, cancer_signa_corr2Exon)
fit_COSMICcorr_CDS <- fit_to_signatures(mut_mat_CDS, cancer_signa_corr2CDS)

# Third, correction is applied to Exon or CDS region -> WGS. No manipulation to COSMIC
rownames(mut_mat_Exon) <- ssm96_somsign
rownames(mut_mat_CDS) <- ssm96_somsign

mut_mat_Exon_corr <- normalizeMotifs(mut_mat_Exon, corr_factorExon2WGS)
mut_mat_CDS_corr <- normalizeMotifs(mut_mat_CDS, corr_factorCDS2WGS)

rownames(mut_mat_Exon) <- ssm96_std_rownames
rownames(mut_mat_CDS) <- ssm96_std_rownames
rownames(mut_mat_Exon_corr) <- ssm96_std_rownames
rownames(mut_mat_CDS_corr) <- ssm96_std_rownames

# Fit the corrected signatures
fit_COSMIC_Exon_corr <- fit_to_signatures(mut_mat_Exon_corr, cancer_signatures)
fit_COSMIC_CDS_corr <- fit_to_signatures(mut_mat_CDS_corr, cancer_signatures)

############
# Results Inspection
############

# Let's put all results on a list so we do not have to repeat the same commands again and again

allfits <- list(WGS = fit_COSMIC_WGS, 
                CDS = fit_COSMIC_CDS, CDS_corr = fit_COSMIC_CDS_corr, COSMICcorr_CDS = fit_COSMICcorr_CDS,
                Exon = fit_COSMIC_Exon, Exon_corr = fit_COSMIC_Exon_corr, COSMICcorr_Exon = fit_COSMICcorr_Exon)

corresponding_mutmats <- list(mut_mat_WGS = mut_mat_WGS,
                              mut_mat_CDS = mut_mat_CDS, mut_mat_CDS_corr = mut_mat_CDS_corr, mut_mat_CDS_COSMcorr = mut_mat_CDS,
                              mut_mat_Exon = mut_mat_Exon, mut_mat_Exon_corr = mut_mat_Exon_corr, mut_mat_Exon_COSMcorr = mut_mat_Exon)

corresponding_cancer_sign <- list(COSMIC_WGS = cancer_signatures,
                              COSMIC_CDS = cancer_signatures, COSMIC_CDS_corr = cancer_signatures, CDS_COSMcorr = cancer_signa_corr2CDS,
                              COSMIC_Exon = cancer_signatures, COSMIC_Exon_corr = cancer_signatures, Exon_COSMcorr = cancer_signa_corr2Exon)

# First, we can compute the RSS and cos_sim for all samples in each case in order to assess
# which of the three diff methods for Exon and CDS is offering best results.

# Calculate all the cos sim between all original and reconstructed mutational profiles

cos_sim_ori_rec <- mapply(function(x,y) as.data.frame(diag(cos_sim_matrix(x, y$reconstructed))),
                          x = corresponding_mutmats,
                          y = allfits)

names(cos_sim_ori_rec) <- names(allfits)
cos_sim_ori_rec <- as.data.frame(cos_sim_ori_rec)
rownames(cos_sim_ori_rec) <- colnames(mut_mat_WGS)
cos_sim_ori_rec$samples <- rownames(cos_sim_ori_rec)

toplot <- melt(cos_sim_ori_rec)

# Let's plot results

pdf("Cos_Sim_Original_Recons.pdf")

ggplot(toplot, aes(y=value, x=samples)) +
  facet_wrap(~ variable) +
  geom_bar(stat="identity", fill = "skyblue4") +
  coord_cartesian(ylim=c(0.8, 1)) +
  ylab("Cosine similarity\n original VS reconstructed") +
  theme_bw() +
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank(),
        axis.text.x = element_text(angle = 90, size=4, vjust=0.5)) +
  geom_hline(aes(yintercept=.95))

dev.off()

# Calculate all the RSS between all original and reconstructed mutational profiles

RSSs <- mapply(function(x,y) {
  s1_relative <- sweep(as.matrix(x), 2, colSums(as.matrix(x)), `/`)
  s2_relative <- sweep(y$reconstructed, 2, colSums(y$reconstructed), `/`)
  diff <- s1_relative - s2_relative
  RSS <- colSums(diff^2) #RSS por muestra
  #RSS <- format(RSS, scientific = TRUE, digits = 3)
  return(RSS)
  }, x = corresponding_mutmats,
  y = allfits)

RSSs <- as.data.frame(RSSs)
RSSs$samples <- rownames(RSSs)
toplot <- melt(RSSs)

pdf("RSS_Original_Recons.pdf")

ggplot(toplot, aes(y=value, x=samples)) +
  facet_wrap(~ variable, scales="free") +
  geom_bar(stat="identity", fill = "skyblue4") +
  ylab("Error - Sum of Squares\n original VS reconstructed") +
  theme_bw() +
  theme(panel.grid.minor.y=element_blank(),panel.grid.major.y=element_blank(),
        axis.text.x = element_text(angle = 90, size=4, vjust=0.5)) +
  geom_hline(aes(yintercept=0.001), color="red")

dev.off()

# Some other plots regarding signatures contribution

for(i in 1: length(allfits))
{
  pdf(paste0("Heatmap_",names(allfits)[i],".pdf",sep=""))
  p <- plot_contribution_heatmap(allfits[[i]]$contribution, cluster_samples = FALSE)
  plot(p)
  dev.off()
}

# Is there any signature that contributes in all samples?

lapply(allfits, function(fit) which(apply(fit$contribution,1, function(row) all(row!=0))))

# Compute the Signatures Contribution Summary: not a good idea, values are not normal so plot
# mean and sd is not meaningful. A boxplot will be graphed

source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_summary_SignaturesContribution.R")

for(i in 1: 1)#length(allfits))
{
  pdf(paste0("SignContr_",names(allfits)[i],".pdf",sep=""))
  p <- plot_summary_SignaturesContribution(allfits[[i]]$contribution, whichtumour="Lung-AdenoCA", norm=TRUE,data=names(allfits)[i]) 
  plot(p)
  dev.off()
}

# Plot the results per signature
require(ggpubr)
source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_multiple_Sign_fits.R")

pdf("SignContribution_facetSignature.pdf",height=20, width=20)
plot_multiple_Sign_fits(allfits,"Lung-AdenoCA",norm=TRUE)
dev.off()


# Method TWO: whichSignatures function from deconstructSigs ---------------

decS_fit_COSMICcorr_CDS <- sapply(colnames(mut_mat_CDS), function(sample)  
  {compute <- whichSignatures(tumor.ref = as.data.frame(t(mut_mat_CDS)), #Expected to have samples in rows 
                  sample.id = sample, #Ponerlo en un loop
                  #Package already has the cosmic signatures, but can be included (we'll use ours)
                  signatures.ref = as.data.frame(t(cancer_signa_corr2CDS)), #Signatures in rows in case of inclusion
                  associated = c(), #If we wanna narrow the list of signatures
                  signatures.limit = NA, 
                  signature.cutoff = 0, #Weight less than this amount (this value is the one by default)
                  contexts.needed = TRUE, #Since the data input is in counts => This is for normalizing
                  tri.counts.method = "default")
   return(t(compute$weights))})

rownames(decS_fit_COSMICcorr_CDS) <- c(paste(rep("Signature.",30), seq(1,30,1), sep=""))


# For tri.counts.method we could use i.e 'default' (no scaling factor) or 'exome2genome' or 'genome2exome'
# We could also pass a data frame with the proper scaling factor, in any case it should be:
# a data frame should match the format of ‘tri.counts.exome‘ or ‘tri.counts.genome‘

allfits_decSigs <- list(WGS = decS_fit_WGS, 
                CDS = decS_fit_CDS, CDS_corr =  decS_fit_CDScorr, COSMICcorr_CDS = decS_fit_COSMICcorr_CDS,
                Exon = decS_fit_Exon, Exon_corr = decS_fit_Exoncorr, COSMICcorr_Exon = decS_fit_COSMICcorr_Exon)

for(i in 1: length(allfits_decSigs))
{
  pdf(paste0("SignContr_decSigs_",names(allfits_decSigs)[i],".pdf",sep=""))
  p <- plot_summary_SignaturesContribution(allfits_decSigs[[i]], whichtumour="Lung-AdenoCA", norm=FALSE,data=names(allfits_decSigs)[i]) 
  plot(p)
  dev.off()
}

#Individual Signature contribution

pdf("SignContri_decSigs_facetSignature.pdf",height=20, width=20)
plot_multiple_Sign_fits(allfits_decSigs,"Lung-AdenoCA", norm=FALSE)
dev.off()


# Method THREE: findSigExposures function from SignatureEstimation -------------

# We will use QP method only (not SA)

source("/project/devel/PCAWG/mmaqueda/Rscripts/decomposeQP.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/findSigExposures.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/FrobeniusNorm.R")

#norm == TRUE para normalizar input
decQP_WGS <- findSigExposures(M = mut_mat_WGS, P = cancer_signatures, decomposition.method = decomposeQP, norm=TRUE) 

allfits_decQP <- list(WGS = decQP_WGS, 
                      CDS = decQP_CDS, CDS_corr =  decQP_CDScorr, COSMICcorr_CDS = decQP_COSMICcorr_CDS,
                      Exon = decQP_Exon, Exon_corr = decQP_Exoncorr, COSMICcorr_Exon = decQP_COSMICcorr_Exon)

pdf("SignContri_decQP_facetSignature.pdf",height=20, width=20)
plot_multiple_Sign_fits(allfits_decQP,"Lung-AdenoCA", norm=FALSE,decQP = TRUE)
dev.off()


# Plotting all fit results simultaneously -------------------------------------

source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_Sign_fits_methods.R")

pdf("SignCont_ALLMETHODS_facetSignature.pdf",height=20, width=20)
plot_Sign_fits_methods(fits_lsqnonneg = allfits, 
                       fits_decSigs = allfits_decSigs, 
                       fits_QP = allfits_decQP, whichtumour="Lung-AdenoCA")
dev.off()

# In case we have the results split per tumor type and each tumour type is in a RData

files <- list.files(path="/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/",
           pattern="RData", full.names = TRUE)

names <- list.files(path="/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/",
                    pattern="RData", full.names = FALSE)
names <- sapply(names, function(x) unlist(strsplit(x, "_"))[3])
names <- sapply(names, function(x) unlist(strsplit(x, ".",fixed=TRUE))[1])

for(i in 1:length(files))
{
  load(files[i])
  
  pdf(paste0("SignCont_ALLMETHODS_perSignature_",names[i],".pdf", sep=""),
      height=20, width=20)
  p <- plot_Sign_fits_methods(fits_lsqnonneg = fits_m1, 
                         fits_decSigs = fits_m2, 
                         fits_QP = fits_m3, whichtumour=names[i])
  plot(p)
  dev.off()
}


# Compute decomposition (Frobenius norm) error ----------------------------

# The objective is to compute the error for each of the three methods and plot them

# First we compute the errors and prepare the information for plotting

source("/project/devel/PCAWG/mmaqueda/Rscripts/get_errors_diff_methods.R")
errors_methods <- get_errors_diff_methods(fits_lsqnonneg = allfits,
                        fits_decSigs = allfits_decSigs,
                        fits_decQP = allfits_decQP,
                        mut_mats = corresponding_mutmats, 
                        cancer_signatures = corresponding_cancer_sign)

# Now plot the errors

#allresults$Data <- factor(allresults$Data, levels = names(allfits))

whichtumour <- "Lung-AdenoCA"

pdf("Errors_ALLMETHODS.pdf",height=14, width=16)
ggplot(errors_methods, aes(y=Error, x=Data, fill=method),show.legend=F) + 
  geom_boxplot() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12,vjust=0.5),
                     axis.text.y = element_text(size=16),
                     plot.title = element_text(face = "bold", size = (15))) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "kruskal.test", paired=TRUE) +
  labs(title=paste("Decomposition Error (Frobenius norm) in ",whichtumour, " dataset \n (Kruskal test comparison between methods)"))
dev.off()


# In case we have the results split per tumor type and each tumour type is in a RData

for(i in 1:length(files))
{
  load(files[i])
  
  pdf(paste0("Errors_ALLMETHODS_",names[i],".pdf", sep=""),
      height=10, width=10)
  p <- ggplot(errors_methods, aes(y=Error, x=Data, fill=method),show.legend=F) + 
    geom_boxplot() +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, size=15,vjust=0.5),
                       axis.text.y = element_text(size=18),
                       plot.title = element_text(face = "bold", size = (15))) +
    stat_compare_means(aes(label = ..p.signif..), size=8,
                       method = "kruskal.test", paired=TRUE) +
    labs(title=paste("Decomposition Error (Frobenius norm) in ",names[i], 
                     " dataset \n (Kruskal test comparison between methods)")) 
    
  plot(p)
  dev.off()
}


