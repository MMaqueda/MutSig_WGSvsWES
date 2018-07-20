
Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(parallel))

source("/project/devel/PCAWG/mmaqueda/Rscripts/get_ssm96.R")

setwd("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/")

#All tumour types
tumours <- c("Breast-DCIS", "Breast-AdenoCA", "Breast-LobularCA", "ColoRect-AdenoCA", "Prost-AdenoCA", "Eso-AdenoCA", "CNS-GBM", "Head-SCC", "Kidney-RCC", "Stomach-AdenoCA", "Biliary-AdenoCA", "Liver-HCC", "Thy-AdenoCA", "SoftTissue-Liposarc", "SoftTissue-Leiomyo", "Kidney-ChRCC", "Skin-Melanoma", "CNS-Oligo", "Lymph-BNHL", "Panc-AdenoCA", "Myeloid-AML", "Lung-AdenoCA", "Lung-SCC", "Ovary-AdenoCA", "CNS-Medullo", "CNS-PiloAstro", "Uterus-AdenoCA", "Bladder-TCC", "Cervix-SCC", "Panc-Endocrine", "Cervix-AdenoCA", "Lymph-CLL", "Bone-Epith", "Bone-Benign", "Bone-Osteosarc", "Myeloid-MPN", "Myeloid-MDS")

for(i in 1:length(tumours))
{
  get_ssm96(resultsPath = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/",
            VCFfilesPath = "/project/devel/PCAWG/mmaqueda/WES_VCFfiles/VCF_Filtered/EXON/",
            InfoVCFs.RData = "/project/devel/PCAWG/somatic_consensus_aug2016/data/t_sample2tumortype_aliquot.RData", 
            ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
            whichtumour = tumours[i],
            region="Exon")
}


