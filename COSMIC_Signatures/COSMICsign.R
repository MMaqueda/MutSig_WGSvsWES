
#This script is just to obtain the COSMIC signatures and leave them prepare for further analysis

#Contribution of COSMIC signatures
#sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
#cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)


#We need to load one mut matrix to match the order to SSMs
mut96_mat <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)

# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)

# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

save.image("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/COSMICsign.RData")

