
# This script is to retrieve COSMIC signatures contribution (with NNLS method)
# from MSI tumors  - results are distributed in different RData files

load("/project/devel/PCAWG/stats_replTime2cluster2drivers.RData")
MSI_samples <- stats_replTime2cluster2drivers[which(stats_replTime2cluster2drivers[,"msi_riken"] %in% 1),c("sample_id","tumor_type")]

labels <- character()
for(j in 1:dim(MSI_samples)[1])
{labels[j] <- paste(MSI_samples$sample_id[j], MSI_samples$tumor_type[j],sep=".")}

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Ovary-AdenoCA",".RData",sep=""))
WGS <- fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[1])]
CDS <- fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[1])]
CDScorr <- fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[1])]

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Kidney-RCC",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[2])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[2])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[2])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Liver-HCC",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[3])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[3])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[3])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Skin-Melanoma",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[4])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[4])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[4])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Panc-AdenoCA",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[5])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[5])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[5])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "ColoRect-AdenoCA",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[6])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[6])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[6])])

WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[7])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[7])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[7])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Stomach-AdenoCA",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[8])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[8])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[8])])

WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[9])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[9])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[9])])

WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[10])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[10])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[10])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Uterus-AdenoCA",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[11])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[11])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[11])])

WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[12])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[12])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[12])])

WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[13])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[13])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[13])])

WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[14])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[14])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[14])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Liver-HCC",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[15])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[15])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[15])])

load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Biliary-AdenoCA",".RData",sep=""))
WGS <- cbind(WGS,fits_m1$WGS$contribution[,which(colnames(fits_m1$WGS$contribution) == labels[16])])
CDS <- cbind(CDS, fits_m1$CDS$contribution[,which(colnames(fits_m1$CDS$contribution) == labels[16])])
CDScorr <- cbind(CDScorr, fits_m1$CDS_corr$contribution[,which(colnames(fits_m1$CDS_corr$contribution) == labels[16])])

plot_summary_SignaturesContribution_WGS_CDS(WGS, CDS, CDScorr, "MSI", norm=TRUE)