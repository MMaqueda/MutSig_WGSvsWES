
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(R.utils)) #Para el sourceDirectory
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pls))

#Functions needed in this script
sourceDirectory("/project/devel/PCAWG/mmaqueda/Rscripts/")
load("/project/devel/PCAWG/mmaqueda/Rscripts/allRscripts.RData")

# STANDARD PIPELINE to get the signatures of a specific tumour type

# Prepare Data for specific Tumour Type -----------------------------------

#(1) Get the mutational profile (SSM96) for the tumour type of interest

mutSSM96 <- select_tumour(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/mut96_matALL.txt",
     resultsPath = "/project/devel/PCAWG/mmaqueda/Results/ColoRect/Data/",
     whichtumour = c("ColoRect-AdenoCA"))


# Plotting Tumour Type Data -----------------------------------------------
#(2) Plot a heatmap with the SSM96 mutational profile - all samples

pdf("Heatmap_SSM96.pdf")
plot_heatmap_SSM96(mutSSM96)
dev.off()

#(3) Plot the SSM6 spectrum (counts and freq) for all samples (together and separate)
#    A prior step is to get the SSM6 matrix

mutSSM6 <- SSM96toSSM6matrix(Path_mut96_mat = 
                                     "/project/devel/PCAWG/mmaqueda/Results/ColoRect/Data/mut96_mat_ColoRect-AdenoCA.txt")

pdf("Spectrum_SSM6_CountsSamples.pdf")
plot_spectrum_SSM6(mutSSM6 = mutSSM6, plotby="SAMPLE", 
                        legend = TRUE, relative=F)
dev.off()

pdf("Spectrum_SSM6_FreqSamples.pdf")
plot_spectrum_SSM6(mutSSM6 = mutSSM6, plotby="SAMPLE", 
                        legend = TRUE, relative=T)
dev.off()

#    This plot function is from the MutPatt package
pdf("Spectrum_SSM6_FreqALL.pdf")
plot_spectrum(as.data.frame(t(mutSSM6)))
dev.off()

#(4) Plot the SSM96 profile for 10 samples (randomly chosen)
#    This plot function is from the MutPatt package

pdf("SSM96_RandomSamples.pdf")
plot_96_profile(mutSSM96[,sample(1:dim(mutSSM96)[2],size=10)], condensed = TRUE)
dev.off()

#(5) Plot the first 4PC of PCA over all samples (SSM96 info)
#    A prior step is to normalized (freq not counts) the SSM96 matrix

mutSSM96_norm <- normalize_ssm96mat(Path_mut96_mat = 
                                      "/project/devel/PCAWG/mmaqueda/Results/ColoRect/Data/mut96_mat_ColoRect-AdenoCA.txt",
                                    resultsPath = "/project/devel/PCAWG/mmaqueda/Results/ColoRect/Data/")


pdf("PCA_SSM96_Norm.pdf")
plot_PCA_SSM96(mutSSM96_norm)
dev.off()

#(6) Compute the hierarchical clustering (complete method) over all samples
#    and plot the dendrogram

pdf("Dendrogram.pdf")
plot(cluster_SSM96profile(mutSSM96))
dev.off()

# Estimate NMF for different k signatures ---------------------------------
#(7) Estimate the max number of signatures that can be retrieve with
#    the number of genomes to be analyzed

kmax <- EstMaxNumOfSignatures(n.genomes = dim(mutSSM96)[2])

#(8) Estimate NMF for different k values considering kmax from previous point
# We will consider 100 iterations  

estimate <- nmf(mutSSM96, rank=2:kmax, method="brunet", nrun=100, seed=123456,.options="kvp2")

# Selection of optimal k value from NMF -----------------------------------
#(9) First, let's plot the different default measures from nmf function

pdf("Estimatek.pdf")
plot(estimate)
dev.off()

#(10a) Selection of k values based on CCC

selCCC <-  kSelectionCCC(nmf_measures = estimate$measures,cutoffCCCrel=0.85)

#(10b) Selection of k values based on RSS

selRSS <- kSelectionRSS(nmf_measures = estimate$measures,cutoffRSSrel=0.80)

#(10c) Selection of k values based on Frobenius norm stability

selFrob <- kSelectionStableFrob(mut96_mat = mutSSM96, 
                                nmf_estimation = estimate) 

# selFrob[[1]] #selected

pdf("Stability_Frobnorm.pdf")
boxplot(selFrob[[2]],use.cols=TRUE,main="Lung SCC - Frobenius Norm")  #plot Frob norm values
dev.off()

#(10d) Selection of k values based on similarity Original and Reconstructed

#First, let's do the wrap-up
NMFks_results <- wrapup_results_NMFks(selFrob,estimate)

#Cos Sim for all k's and samples
CSOriReco_ks <- kCosSimOriReco(mutSSM96, NMFks_results)

#Summary of Cos Sim values
selSimOriRec <- summaryCosSim(CSOriReco_ks,estimate)

pdf("summCS_OriRec.pdf")
plot_summaryCosSim_OrRec(selSimOriRec)
dev.off()

#(11) Selection of optimum k value - Consensus among previous metrics

optimumK <- kSelection_Consensus(selCCC,selRSS,selFrob,selSimOriRec,estimate$measures)

# optimumK[[1]] #selected

pdf("Votes4K.pdf")
barplot(optimumK[[2]],main="K votes")  
dev.off()

# Cos Sim Values per k candidates -----------------------------------------

#After visualizing previous results, we'll plot this for k=...
pdf("CS_OriRec_k2.pdf")
plot_CS_OriRec_forK(CSOriReco_ks,kvalue=2,estimate$measures)
dev.off()

# Retrieve results for optimum k ------------------------------------------
#Based on optimumK[[1]] (selected k value), in this case k=2

#(12) Plot the samples dendrogram based on the consensus matrix from NMF estimation
#     for the k selected

pdf("consensushc_k5.pdf")
plot(consensushc(object = estimate$fit$`5`, method="complete",what="consensus",dendrogram=F))
dev.off()

#(13) Plot the contribution of each signature per sample in rel and abs mode for the
#     k selected
results_k_selected <- NMFks_results[[which(estimate$measures$rank == 5)]] #optimumK[[1]])]]

pc1 <- plot_contribution(results_k_selected$contribution, 
                         results_k_selected$signatures,mode = "relative")
pc2 <- plot_contribution(results_k_selected$contribution, 
                         results_k_selected$signatures,mode = "absolute")

pdf("Contribution_K5.pdf")
grid.arrange(pc1, pc2)
dev.off()

#(14) Plot the SSM96 profile of the signatures

pdf("Signatures_K5.pdf")
plot_96_profile(results_k_selected$signatures, condensed = TRUE)
dev.off()

#(15) Plot the signatures contribution per sample in a heatmap

pdf("Contribution_Heat_K5.pdf")
plot_contribution_heatmap(results_k_selected$contribution)
dev.off()




