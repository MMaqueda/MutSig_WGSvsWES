
#setwd("/project/devel/PCAWG/mmaqueda/Results/Samples_Visualization/WGS/SSM/SpectrumSSM6/")
#setwd("/project/devel/PCAWG/mmaqueda/Results/Samples_Visualization/WGS/SSM/Heatmaps/")

# Data Visualization ------------------------------------------------------

# List of tumors

tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS")

# Import mutational profile for all the tumors types
#mut96_mat_WGS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)  
#mut96_mat_Exon <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_allExon.txt",header=TRUE,check.names=FALSE)  
#mut96_mat_CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE)  


# Convert the SSM96 matrix to SSM6 profile
source("/project/devel/PCAWG/mmaqueda/Rscripts/SSM96toSSM6matrix.R")

mut6_WGS <- SSM96toSSM6matrix(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt")
mut6_Exon <- SSM96toSSM6matrix(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_allExon.txt")
mut6_CDS <- SSM96toSSM6matrix(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt")

# Let's put that information in a list
#mut96_mat <- list(WGS = mut96_mat_WGS, Exon = mut96_mat_Exon, CDS = mut96_mat_CDS) 
  
####
# SSM6 Spectrum
####

source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_spectrum_SSM6.R")

pdf("SSM6_all_CDS.pdf",height = 18, width=18)
plot_spectrum_SSM6(mut6_CDS, plotby="TUMOUR",legend = TRUE, relative=T) 
dev.off()

####
# Heatmap
####

source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_heatmap_SSM96.R")

for(i in 1:length(tumors))
{
  # Import mutational profile for the specific tumor type
  WGS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_",tumors[i],".txt",sep=""),
                    header=TRUE,check.names=FALSE)  
  Exon <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_",tumors[i],".txt",sep=""),
                    header=TRUE,check.names=FALSE)  
  CDS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_",tumors[i],".txt",sep=""),
                    header=TRUE,check.names=FALSE)  
  
  # Plotting
  
  info <- list(WGS = WGS, Exon=Exon, CDS=CDS)
  
  for(j in 1:3)
  {
    p <- plot_heatmap_SSM96(info[[j]])
    
    pdf(paste0("Heatmap_",names(info)[j],"_",tumors[i],".pdf",sep=""),height=16)
    print(p)
    dev.off()  
  }

}

####
# Clustering
####

# The clustering is based on the cosine similarity
source("/project/devel/PCAWG/mmaqueda/Rscripts/cluster_SSM96profile.R")

for(i in 1:length(tumors))
{
  # Import mutational profile for the specific tumor type
  WGS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_",tumors[i],".txt",sep=""),
                    header=TRUE,check.names=FALSE)  
  Exon <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_",tumors[i],".txt",sep=""),
                     header=TRUE,check.names=FALSE)  
  CDS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_",tumors[i],".txt",sep=""),
                    header=TRUE,check.names=FALSE)  
  
  # Plotting
  
  info <- list(WGS = WGS, Exon=Exon, CDS=CDS)
  
  for(j in 1:3)
  {
    p <- cluster_SSM96profile(info[[j]], method = "complete")
    
    pdf(paste0("Cluster_",names(info)[j],"_",tumors[i],".pdf",sep=""))
    plot(p)
    dev.off()  
  }
  
}


# Detect weak Subtypes ----------------------------------------------------

# List those subtypes (from SSM96 profile) that together account for <= 1% of the mutations in all genomes (Alexandrov paper)
# We'll use removeWeak from mutSignatures

params <- data.frame(remove.weak.muttypes=0.01) #Defined like this for directly use of the function

summary_Weak <- as.data.frame(matrix(ncol=4,nrow=length(tumors)))
colnames(summary_Weak) <- c("Tumour_Type","WGS","Exon","CDS")
summary_Weak$Tumour_Type <- tumors

summary_Weak_simpl <- as.data.frame(matrix(ncol=4,nrow=length(tumors)))
colnames(summary_Weak_simpl) <- c("Tumour_Type","WGS","Exon","CDS")
summary_Weak_simpl$Tumour_Type <- tumors

for(i in 1:length(tumors))
{
  # Import mutational profile for the specific tumor type
  WGS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_",tumors[i],".txt",sep=""),
                    header=TRUE,check.names=FALSE)  
  Exon <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_",tumors[i],".txt",sep=""),
                     header=TRUE,check.names=FALSE)  
  CDS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_",tumors[i],".txt",sep=""),
                    header=TRUE,check.names=FALSE) 
  
  # Detecting the weakest
  
  info <- list(WGS = WGS, Exon=Exon, CDS=CDS)
  
  for(j in 1:3)
  {
    weak <- removeWeak(info[[j]], params)
    summary_Weak[i,j+1] <- paste0(rownames(info[[j]])[weak$removed.mutset],collapse=";")
    
    summary_Weak_simpl[i,j+1] <- length(weak$removed.mutset)
  }
}



