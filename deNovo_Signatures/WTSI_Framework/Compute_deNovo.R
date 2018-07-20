
require(mutSignatures)
require(MutationalPatterns)
source("/project/devel/PCAWG/mmaqueda/Rscripts/FrobeniusNorm.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/EstMaxNum_Signatures.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/decipherMutationalProcesses_wSD.R")

# General pipeline for computing deNovo signatures with mutSignatures package

# Define the list of tumour types to analyze
# tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
#             "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
#             "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
#             "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
#             "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS")

tumors <- c("Lung-AdenoCA")

# Set the initial parameters for the process
# The minimum number of processes to compute is 2: we are assuming there are at least 2 underlying processes
# If 1 has to be computed, the function evaluateStability must be modified to accept one cluster

params <- setMutClusterParams(num.processes.toextract = 2, tot.iterations = 100, 
                              tot.cores = 4, remove.weak.muttypes = 0.01, 
                              remove.last.percent = 0.07, process.distance = "cosine", 
                              tot.Replicates = 100, eps = 2.2204e-16, 
                              stopconv = 10000, niter = 1e+06)


for(i in 1:length(tumors))
{
  # Retrieve data
  mutCounts <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_",tumors[i],".txt",sep=""),header=TRUE,check.names=FALSE)
  
  # Prepare the input object for the whole process
  samples <- colnames(mutCounts)
  colnames(mutCounts) <- NULL
  
  input <- setMutCountObject(mutCountMatrix =  as.matrix(mutCounts), 
                             mutationTypes = rownames(mutCounts), 
                             sampleNames = samples, 
                             datasetName = tumors[i])
  
  # Estimate the maximum number of signatures for that specific dataset
  max_signatures <- EstMaxNumOfSignatures(dim(mutCounts)[2])
  
  # Initialize list of results
  results <- list()
  
  # Sequentially deciphering signatures between min and max_signatures
  for(j in 2:max_signatures)
  {
    # Update the number of processes to extract
    params$num.processes.toextract <- j
    
    # Decipher mutational processes
    results[[j-1]] <- decipherMutationalProcesses_wSD(input, params)
    save(results, file=paste0("/project/devel/PCAWG/mmaqueda/Results/deNovo_signatures/Lung-AdenoCA/CDS_", tumors[i],".RData",sep=""))
  }
  
}
