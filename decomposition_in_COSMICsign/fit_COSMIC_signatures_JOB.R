
# Identify the COSMIC signatures contribution in a SSM96 mutational profile from three different methods

.libPaths(c(.libPaths(), "/home/devel/mmaqueda/R/x86_64-pc-linux-gnu-library/3.4"))

require(MutationalPatterns)
require(deconstructSigs)
require(ggplot2)
require(reshape)
require(SomaticSignatures)
source("/project/devel/PCAWG/mmaqueda/Rscripts/select_tumour.R")
#This is for applying the correction factor. Problems to install SomaticSignatures package. Function from there 
#source("/project/devel/PCAWG/mmaqueda/Rscripts/normalizeMotifs.R")

#Functions for method3
source("/project/devel/PCAWG/mmaqueda/Rscripts/decomposeQP.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/findSigExposures.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/FrobeniusNorm.R")

#Error methods
source("/project/devel/PCAWG/mmaqueda/Rscripts/get_errors_diff_methods.R")

# Load the 30 COSMIC signatures (matrix 96 x30). Object called 'cancer_signatures'
load("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/COSMICsign.RData")

# Load the correction factor from cds2wgs and exon2wgs
load("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Correction_Factors.RData")

# Vector with the name of the tumor types

tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "ColoRect-AdenoCA")


tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS")

# In order to use the correction function we need to renamed the rownames as in the SomSign package

subs <-  c('CA','CG','CT','TA','TC','TG')
context <- c('A.A','A.C','A.G','A.T','C.A','C.C','C.G','C.T',
             'G.A','G.C','G.G','G.T','T.A','T.C','T.G','T.T')

ssm96_somsign <- paste(rep(subs,each=16), rep(context,times=4), sep=" ")

#####
# Loop for all tumor types
#####

for(i in 1:length(tumors))
{
  ##########
  # Prepare the needed data
  ##########
  
  # Select the specific mut matrix per tumour from the complete SSM96 matrix
  mut_mat_WGS <- select_tumour(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt",
                               resultsPath = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/",
                               whichtumour = tumors[i])
  
  # In case  of Exon and CDS we already have the specific mut_matrix
  mut_mat_Exon <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_",tumors[i],".txt",sep=""),
                             header=TRUE,check.names=FALSE)

  mut_mat_CDS <- read.table(paste0("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_",tumors[i],".txt",sep=""),
                            header=TRUE,check.names=FALSE)
  
  
  # Those samples that have zero mutations in all the SSM96 are removed a priori from this analysis.
  # The CDS group will be the one limiting. So let's check this one.
  
  
  ########
  # Apply Method ONE: lsqnonneg, fit_to_signatures function from MutationalPatterns 
  ########

  # This method asks for counts (not freq)

  # First, no correction is applied in Exon or CDS data
  fit_WGS <- fit_to_signatures(mut_mat_WGS, cancer_signatures)
  fit_Exon <- fit_to_signatures(mut_mat_Exon, cancer_signatures)
  fit_CDS <- fit_to_signatures(mut_mat_CDS, cancer_signatures)
  
  # Second, correction is applied to Exon or CDS region -> WGS. No manipulation to COSMIC
  ssm96_std_rownames <- rownames(mut_mat_Exon)
  rownames(mut_mat_Exon) <- ssm96_somsign
  rownames(mut_mat_CDS) <- ssm96_somsign
  
  mut_mat_Exon_corr <- normalizeMotifs(mut_mat_Exon, corr_factorExon2WGS)
  mut_mat_CDS_corr <- normalizeMotifs(mut_mat_CDS, corr_factorCDS2WGS)
  
  # Coming back to the standard nomenclature
  rownames(mut_mat_Exon) <- ssm96_std_rownames
  rownames(mut_mat_CDS) <- ssm96_std_rownames
  rownames(mut_mat_Exon_corr) <- ssm96_std_rownames
  rownames(mut_mat_CDS_corr) <- ssm96_std_rownames
  
  # Fit signatures into corrected data
  fit_Exon_corr <- fit_to_signatures(mut_mat_Exon_corr, cancer_signatures)
  fit_CDS_corr <- fit_to_signatures(mut_mat_CDS_corr, cancer_signatures)
  
  #Let's collide results
  fits_m1 <- list(WGS = fit_WGS, 
                  CDS = fit_CDS, CDS_corr = fit_CDS_corr,
                  Exon = fit_Exon, Exon_corr = fit_Exon_corr)
  
  # The following object will be valid for all three methods
  corresponding_mutmats <- list(mut_mat_WGS = mut_mat_WGS,
                                mut_mat_CDS = mut_mat_CDS, mut_mat_CDS_corr = mut_mat_CDS_corr,
                                mut_mat_Exon = mut_mat_Exon, mut_mat_Exon_corr = mut_mat_Exon_corr)
  
  
  ########
  # Apply Method TWO: heuristic approach, whichSignatures function from deconstructSigs
  ########
  
  fits_m2 <- list()
  contri <- as.data.frame(matrix(ncol=dim(corresponding_mutmats[[1]])[2],nrow=30))
  
  for(j in 1:length(corresponding_mutmats))
  {
    for(k in 1:dim(corresponding_mutmats[[j]])[2])
    {
      if(all(corresponding_mutmats[[j]][k]== 0) == TRUE)  #A column of zeros
      {
        contri[,k] <- rep(NA,30) #Zero value for all signatures
      }else 
      {
        compute <- whichSignatures(tumor.ref = as.data.frame(t(corresponding_mutmats[[j]])), #Expected to have samples in rows 
                                   sample.id = colnames(corresponding_mutmats[[j]][k]),
                                   #Package already has the cosmic signatures, but can be included (we'll use ours)
                                   signatures.ref = as.data.frame(t(cancer_signatures)), #Signatures in rows in case of inclusion
                                   associated = c(), #If we wanna narrow the list of signatures
                                   signatures.limit = NA, 
                                   signature.cutoff = 0, #Weight less than this amount (this value is the one by default)
                                   contexts.needed = TRUE, #Since the data input is in counts => This is for normalizing
                                   tri.counts.method = "default")
        contri[,k] <- t(compute$weights)}
    }
    rownames(contri) <- c(paste(rep("Signature.",30), seq(1,30,1), sep=""))
    fits_m2[[j]] <- contri
    names(fits_m2)[j] <- names(fits_m1)[j]  #Order is the same as in method1
  }
  
  ########
  # Apply Method THREE: findSigExposures function from SignatureEstimation
  ########
  # We will use QP method only (not SA)
  
  #norm == TRUE para normalizar input
  fits_m3 <- lapply(corresponding_mutmats, function(profile)
    {
    allzeros <- which(apply(profile, 2, function(sample) all(sample==0))==TRUE)
    if(length(allzeros) > 0 ) # there's some sample with all zeros
      {profile[allzeros] <- 5} #Fake the data in this column/columns to run the method
    
    decomp <- findSigExposures(M = profile, P = cancer_signatures, decomposition.method = decomposeQP, norm=TRUE)
    
    if(length(allzeros) >0) #And now fix the faked data
    {
      decomp$exposures[,allzeros] <- NA
      decomp$errors[allzeros] <- NA
    }
    return(decomp)})
  
  names(fits_m3) <- names(fits_m1) #The order of the elements is the same

  #######
  # Computing the Reconstruction Error in each case
  #######

  errors_methods <- get_errors_diff_methods(fits_lsqnonneg = fits_m1,
                                            fits_decSigs = fits_m2,
                                            fits_decQP = fits_m3,
                                            mut_mats = corresponding_mutmats, 
                                            cancer_signatures = cancer_signatures)
  
  save(fits_m1, fits_m2, fits_m3, errors_methods,
             file=paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", tumors[i],".RData",sep=""))
}


