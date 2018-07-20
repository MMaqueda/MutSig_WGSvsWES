
require(mutSignatures)
require(MutationalPatterns)
source("/project/devel/PCAWG/mmaqueda/Rscripts/FrobeniusNorm.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/EstMaxNum_Signatures.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/decipherMutationalProcesses_wSD.R")

mut96_mat_WGS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)  
mut96_mat_CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE) 

# General pipeline for computing deNovo signatures with mutSignatures package
# MSI tumors

load("/project/devel/PCAWG/stats_replTime2cluster2drivers.RData")
MSI_samples <- stats_replTime2cluster2drivers[which(stats_replTime2cluster2drivers[,"msi_riken"] %in% 1),c("sample_id","tumor_type")]

labels <- character()
for(j in 1:dim(MSI_samples)[1])
{labels[j] <- paste(MSI_samples$sample_id[j], MSI_samples$tumor_type[j],sep=".")}

mutCounts <- mut96_mat_CDS[, which(colnames(mut96_mat_CDS) %in% labels)]
#mutCounts <- mut96_mat_WGS[, which(colnames(mut96_mat_WGS) %in% labels)]

# Set the initial parameters for the process
# The minimum number of processes to compute is 2: we are assuming there are at least 2 underlying processes
# If 1 has to be computed, the function evaluateStability must be modified to accept one cluster

params <- setMutClusterParams(num.processes.toextract = 2, tot.iterations = 25, 
                              tot.cores = 16, remove.weak.muttypes = 0.01, 
                              remove.last.percent = 0.07, process.distance = "cosine", 
                              tot.Replicates = 100, eps = 2.2204e-16, 
                              stopconv = 10000, niter = 1e+06)


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
  save(results, file=paste0("/project/devel/PCAWG/mmaqueda/Results/deNovo_signatures/MSI/CDS_", "MSI",".RData",sep=""))
}

