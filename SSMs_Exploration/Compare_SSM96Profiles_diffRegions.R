
source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_spectrum_SSM96.R")
require(SomaticSignatures)
source("/project/devel/PCAWG/mmaqueda/Rscripts/normalizeMotifs.R")

# This is for comparing the SSM96 mutational profile from WGS or CDS data (this one corrected)
# and generate the plot of these profiles. We will focus on the profile based on the mean

# Define the list of tumour types to analyze
tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS")


SSM96_Comparison <- as.data.frame(matrix(ncol=3, nrow=length(tumors)))
colnames(SSM96_Comparison) <- c("Tumors","CS_CDSvsWGS","CS_CDScorrvsWGS")

for(i in 1:length(tumors))
{
  # Obtain the mutational profiles for a specific tumor type
  
  CDS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_",tumors[i],".txt",sep=""),header=TRUE,check.names=FALSE)
  WGS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_",tumors[i],".txt",sep=""),header=TRUE,check.names=FALSE)
  
  # SSM96 mean spectrum
  
  meanCDS <-apply(CDS, 1, function(SSM96type)  mean(SSM96type, na.rm=TRUE))
  meanWGS <-apply(WGS, 1, function(SSM96type)  mean(SSM96type, na.rm=TRUE))
  
  # Apply the correction to the CDS profile before comparison
  
  # Load the correction factor from cds2wgs and exon2wgs
  load("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Correction_Factors.RData")
  
  # In order to use the correction function we need to renamed the rownames as in the SomSign package
  
  subs <-  c('CA','CG','CT','TA','TC','TG')
  context <- c('A.A','A.C','A.G','A.T','C.A','C.C','C.G','C.T',
               'G.A','G.C','G.G','G.T','T.A','T.C','T.G','T.T')
  
  ssm96_somsign <- paste(rep(subs,each=16), rep(context,times=4), sep=" ")
  
  # Correction is applied to obtained processes
  names(meanCDS) <- ssm96_somsign
  meanCDS_corr <- normalizeMotifs(as.matrix(meanCDS), corr_factorCDS2WGS)
  
  # Cos Sim to compare both spectrums
  CosSim_WGSvsCDScorr <- as.data.frame(cos_sim_matrix(meanCDS_corr, as.matrix(meanWGS)))
  CosSim_WGSvsCDS <- as.data.frame(cos_sim_matrix(as.matrix(meanCDS), as.matrix(meanWGS)))
  
  # Save results in the dataframe
  SSM96_Comparison[i,] <- c(tumors[i], CosSim_WGSvsCDS, CosSim_WGSvsCDScorr)
  
  # Plot both spectrums
  p <- plot_spectrum_SSM96(WGS, legend = TRUE, relative=F, wSD=F, ymax=NA, metric = "mean", tumors[i], "WGS")
  
  pdf(paste0("WGS_mean_",tumors[i],".pdf", sep=""))
  print(p) 
  dev.off()
  
  p <- plot_spectrum_SSM96(CDS, legend = TRUE, relative=F, wSD=F, ymax=NA, metric = "mean", tumors[i], "CDS") 
  pdf(paste0("CDS_mean_",tumors[i], ".pdf",sep=""))
  print(p)
  dev.off()
}

save(SSM96_Comparison, file="/project/devel/PCAWG/mmaqueda/Results/Samples_Visualization/SSM96_Spectrum/SSM96_Comparison.RData")

