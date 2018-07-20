
# This script allows the computation of the mutation density ratio (mutations per Mb)

# Read the complete mutational profile for WGS, Exon and CDS
WGS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)  
Exon <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_allExon.txt",header=TRUE,check.names=FALSE) 
CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE) 


LungAdenoWGS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_Lung-AdenoCA.txt",header=TRUE,check.names=FALSE)  
LungAdenoCDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_Lung-AdenoCA.txt",header=TRUE,check.names=FALSE) 

# Remember that samples have the same order in all three matrices

# Prepare a data frame where to store all the data.

Mutation_density_ratios <- as.data.frame(matrix(ncol=11, nrow=dim(WGS)[2]))
colnames(Mutation_density_ratios) <- c("Sample_ID", "SSMs_CDS", "SSMs_Exon",
                              "SSMs_WGSwoCDS", "SSMs_WGSwoExon",
                              "SSMsxMb_CDS", "SSMsxMb_Exon",
                              "SSMsxMb_WGSwoCDS", "SSMsxMb_WGSwoExon",
                              "SSMs_CDStoWGSratio", "SSMs_ExontoWGSratio")

# Needed values

CDSregion <- 35194324
Exonregion <- 121892536

WGSwoCDSregion <- 2861327131 - CDSregion
WGSwoExonregion <- 2861327131 - Exonregion

#percentageCDS <- 1.23/100
#percentageExon <- 4.26/100

# Compute mutation density ratio per sample

for(i in 1:dim(WGS)[2])
{
  print(paste0("Sample:",i))
  
  # Obtain number of mutations in each case (SSMs_XXX columns)
  
  Mutation_density_ratios[i,c(1:5)] <- c(colnames(WGS)[i],
                                         colSums(CDS)[i],
                                         colSums(Exon)[i],
                                         colSums(WGS)[i] - colSums(CDS)[i],
                                         colSums(WGS)[i] - colSums(Exon)[i])
  
  
  # Compute the number of mutations per Mb (SSMsxMb_XXX columns)
  
  Mutation_density_ratios[i,c(6:9)] <- c(colSums(CDS)[i]/(CDSregion/1000000),
                                         colSums(Exon)[i]/(Exonregion/1000000),
                                         (colSums(WGS)[i] - colSums(CDS)[i])/(WGSwoCDSregion/1000000),
                                         (colSums(WGS)[i] - colSums(Exon)[i])/(WGSwoExonregion/1000000))
  
  # Compute the mutations density ratio
  
  Mutation_density_ratios[i,c(10:11)] <- c(Mutation_density_ratios[i,"SSMsxMb_CDS"]/Mutation_density_ratios[i,"SSMsxMb_WGSwoCDS"],
                                           Mutation_density_ratios[i,"SSMsxMb_Exon"]/Mutation_density_ratios[i,"SSMsxMb_WGSwoExon"])
  
}  

# We can add the tumor type since we'll need it for the plotting

Tumour_Type <- sapply(Mutation_density_ratios$Sample_ID, function(x) unlist(strsplit(x, "[.]"))[2])
Mutation_density_ratios <- cbind(Mutation_density_ratios[,1],Tumour_Type, Mutation_density_ratios[,2:dim(Mutation_density_ratios)[2]])

save(Mutation_density_ratios,file="/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/MutationDensityRatios.RData")

# Let's plot results
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")

# All the needed data is in Mutation_density_ratios dataframe

orderCDS <- names(tapply(Mutation_density_ratios$SSMs_CDStoWGSratio, Mutation_density_ratios$Tumour_Type, median))[order(tapply(Mutation_density_ratios$SSMs_CDStoWGSratio, Mutation_density_ratios$Tumour_Type, median))]
orderExon <- names(tapply(Mutation_density_ratios$SSMs_ExontoWGSratio, Mutation_density_ratios$Tumour_Type, median))[order(tapply(Mutation_density_ratios$SSMs_ExontoWGSratio, Mutation_density_ratios$Tumour_Type, median))]

Mutation_density_ratios$Tumour_Type <- factor(Mutation_density_ratios$Tumour_Type, levels = orderCDS)
Mutation_density_ratios$Tumour_Type <- factor(Mutation_density_ratios$Tumour_Type, levels = orderExon)

cols_pcawg <- pcawg.colour.palette(x = levels(Mutation_density_ratios$Tumour_Type),scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

boxplot <- ggplot(Mutation_density_ratios, aes(y=SSMs_ExontoWGSratio, x=Tumour_Type, fill=Tumour_Type)) + 
  geom_boxplot() +
  scale_fill_manual(values= cols_pcawg) +
  coord_cartesian(ylim = c(0, 4)) +
  ggtitle("Mutation Density Ratio (SSM per Mb) \n Exon vs WGS-Exon regions") +
  ylab("Mutation Density Ratios Exon region / WGS-Exon region") +
  geom_hline(yintercept = 1, color="red", linetype="dashed") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12,vjust=0.5),
                     axis.title.y=element_text(size=14,face="bold"),
                     axis.title.x=element_blank(),
                     axis.text.y = element_text(size=18),
                     plot.title = element_text(face = "bold", size = (15))) +
  guides(fill=FALSE)

pdf("MutationDensityRatioExon.pdf",width=14, height=14)
boxplot
dev.off()

# Include a graph with the mutation density rates from both ROIs and compared them

#toplot <- melt(Mutation_density_ratios[,c(2,11,12)])
toplot <- melt(MutRatios)
colnames(toplot) <- c("Tumour_Type","Region","Ratio")

# Let's remove outliers for better visualization (above 4 in ratio)
toplot <- toplot[-which(toplot$Ratio >4),]

#Colors
cols_pcawg <- pcawg.colour.palette(x = levels(toplot$Tumour_Type),scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

boxplot <- ggplot(toplot, aes(y=Ratio, x=Tumour_Type, color=Region, fill=Tumour_Type)) + 
  geom_boxplot() + 
  scale_fill_manual(values= rep(cols_pcawg,2)) +
  scale_color_manual(values= c("black","plum4")) +
  labs(x = "Tumor Type", y = "Mutation density ratio ROI region /WGS-ROI region") +
  geom_hline(yintercept = 1, color="red", linetype="dashed") +
  ggtitle("Mutation density Ratios") +
  stat_compare_means(aes(group = Region), label = "p.signif", hide.ns=TRUE,
                     method = "wilcox.test") + #Kruskal.Wallis test by default
  theme_bw() + theme(axis.text.x = element_text(angle=90,size=11,vjust=0.5),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) 

