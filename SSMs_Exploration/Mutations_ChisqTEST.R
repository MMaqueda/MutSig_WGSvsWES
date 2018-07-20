
#' This script computes a Chi-square test in order to test if the number of mutations in a specific region
#' of the genome is different from the rest (WGS-that specific region). We will compute the test with median
#' number of SSMs per tumor type and then per sample (and determine % samples with significant p-value)

# Load the needed data. CHANGE THE PATH IF FILES ARE STORED IN A DIFFERENT PLACE
load("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/SSMs_Depletion_Enrichment_RoIs/MutationDensityRatios.RData")

# Test for median values per cancer type ----------------------------------

# The table should include following data (nrow = TumourTypes):
# N samples per tumor type, median SSMs in CDS, median SSMs in Exon,
# median SSMs in WGS-CDS, median SSMs in WGS-Exon, 
# pval Chisq -  CDS (median)
# pval Chisq -  Exon (median)
# adj pval Chisq -  CDS (median) - Corrected per number of cancer types
# adj pval Chisq -  Exon (median) - Corrected per number of cancer types

Mutation_density_ratios$SSMs_CDS <- as.numeric(Mutation_density_ratios$SSMs_CDS)
Mutation_density_ratios$SSMs_Exon <- as.numeric(Mutation_density_ratios$SSMs_Exon)
Mutation_density_ratios$SSMs_WGSwoCDS  <- as.numeric(Mutation_density_ratios$SSMs_WGSwoCDS)
Mutation_density_ratios$SSMs_WGSwoExon  <- as.numeric(Mutation_density_ratios$SSMs_WGSwoExon)

stats <- function(x) {median = round(median(x),digits=0)}

Mutations_Chisq <- as.data.frame(matrix(ncol=8, 
                                       nrow=length(unique(Mutation_density_ratios$Tumour_Type))))
colnames(Mutations_Chisq) <- c("MedianSSMs_CDS","MedianSSMs_Exon","MedianSSMs_WGSwoCDS","MedianSSMs_WGSwoExon",
                              "pvalChisq_CDS", "pvalChisq_Exon","adjpvalChisq_CDS","adjpvalChisq_Exon")

for(i in 1:4)
{
  # En el df de density ratios están en el mismo orden pero empiezan en la columna 3 (en vez de 1)
  Mutations_Chisq[,i] <- as.data.frame(tapply(Mutation_density_ratios[,i+2], Mutation_density_ratios$Tumour_Type, stats))
}

Mutations_Chisq <- cbind(N=plyr::count(Mutation_density_ratios,'Tumour_Type')$freq,Mutations_Chisq)
Mutations_Chisq <- cbind(Tumour_Type = levels(Mutation_density_ratios$Tumour_Type),Mutations_Chisq)

# We are missing the median of non-mutations in each case. We should compute the values for all samples,
# and then compute the median

CDSregion <- 35194324
Exonregion <- 121892536
WGSwoCDSregion <- 2861327131 - CDSregion
WGSwoExonregion <- 2861327131 - Exonregion

NonMuts <- Mutation_density_ratios[,1:2]
NonMuts$CDSnonmuts <- rep(CDSregion,dim(Mutation_density_ratios)[1]) - Mutation_density_ratios[,3]
NonMuts$Exonnonmuts <- rep(Exonregion,dim(Mutation_density_ratios)[1]) - Mutation_density_ratios[,4]
NonMuts$WGSwoCDSnonm <- rep(WGSwoCDSregion,dim(Mutation_density_ratios)[1]) - Mutation_density_ratios[,5]
NonMuts$WGSwoExonnonm <- rep(WGSwoExonregion,dim(Mutation_density_ratios)[1]) - Mutation_density_ratios[,6]

# Median values for these non-muts

NonMuts_Median <- as.data.frame(matrix(ncol=4, nrow=length(unique(NonMuts$Tumour_Type))))

for(i in 1:4)
{NonMuts_Median[,i] <- as.data.frame(tapply(NonMuts[,i+2], NonMuts$Tumour_Type, stats))}

colnames(NonMuts_Median) <- paste0("Median",colnames(NonMuts)[3:dim(NonMuts)[2]])
NonMuts_Median <- cbind(Tumour_Type = levels(NonMuts$Tumour_Type),NonMuts_Median)

# Compute Chisq test for each tumour type based on the median SSMs values
# This is for testing association betwen muts vs non muts with region

for(i in 1:dim(Mutations_Chisq)[1])
{
  Mutations_Chisq$pvalChisq_CDS[i] <- chisq.test(matrix(c(Mutations_Chisq$MedianSSMs_CDS[i], 
                                                          Mutations_Chisq$MedianSSMs_WGSwoCDS[i], 
                                                          NonMuts_Median$MedianCDSnonmuts[i], 
                                                          NonMuts_Median$MedianWGSwoCDSnonm[i]), 
                                                        ncol = 2))$p.value
    
  Mutations_Chisq$pvalChisq_Exon[i] <- chisq.test(matrix(c(Mutations_Chisq$MedianSSMs_Exon[i], 
                                                        Mutations_Chisq$MedianSSMs_WGSwoExon[i], 
                                                        NonMuts_Median$MedianExonnonmuts[i], 
                                                        NonMuts_Median$MedianWGSwoExonnonm[i]), 
                                                      ncol = 2))$p.value
}  
  
# Get the adjusted p-values
Mutations_Chisq$adjpvalChisq_CDS  <- p.adjust(Mutations_Chisq$pvalChisq_CDS, method="fdr")
Mutations_Chisq$adjpvalChisq_Exon <- p.adjust(Mutations_Chisq$pvalChisq_Exon, method="fdr")


# Fix Table results (adj p-val) for representation ------------------------

Mutations_Chisq_xPublish <- Mutations_Chisq

fix_adjpval <- function(pval)
{
  if(pval<0.001) {
    pval <- "<.001"}else if(pval<0.05) {
      pval <- "<.05"}else if(pval<0.01) {
        pval <- "<.01"}else{pval <- "NS"}
  
  return(pval)
}

Mutations_Chisq_xPublish$adjpvalChisq_CDS <- sapply(Mutations_Chisq_xPublish$adjpvalChisq_CDS, function(pval) fix_adjpval(pval))
Mutations_Chisq_xPublish$adjpvalChisq_Exon <- sapply(Mutations_Chisq_xPublish$adjpvalChisq_Exon, function(pval) fix_adjpval(pval))


# Test for values per sample ----------------------------------------------

# The objective is to apply a Chisq test per sample and then get the % of samples with 
# adj pval <.05  JUST FOR CDS region

Mutations_Chisq_xSample <- NonMuts
colnames(Mutations_Chisq_xSample) <- c("Sample_ID","Tumour_Type",
                                       "pvalCDS","pvalExon",
                                       "adjpvalCDS","adjpvalExon")

for(i in 1:dim(Mutations_Chisq_xSample)[1])
{
  Mutations_Chisq_xSample$pvalCDS[i] <- chisq.test(matrix(c(Mutation_density_ratios$SSMs_CDS[i], 
                                                    Mutation_density_ratios$SSMs_WGSwoCDS[i],
                                                    NonMuts$CDSnonmuts[i],
                                                    NonMuts$WGSwoCDSnonm[i]), 
                                                  ncol = 2))$p.value
  
  Mutations_Chisq_xSample$pvalExon[i] <- chisq.test(matrix(c(Mutation_density_ratios$SSMs_Exon[i], 
                                                 Mutation_density_ratios$SSMs_WGSwoExon[i],
                                                 NonMuts$Exonnonmuts[i],
                                                 NonMuts$WGSwoExonnonm[i]), 
                                               ncol = 2))$p.value
}  

# Get the adjusted p-values
Mutations_Chisq_xSample$adjpvalCDS  <- p.adjust(Mutations_Chisq_xSample$pvalCDS, method="fdr")
Mutations_Chisq_xSample$adjpvalExon <- p.adjust(Mutations_Chisq_xSample$pvalExon, method="fdr")

# And now let's get the percentage of significant samples per Tumour Type

percentage <- function(x,cutoff) {perc = ((length(which(x < cutoff)))/length(x))*100}

# % < .001

#tapply(Mutations_Chisq_xSample[,5], Mutations_Chisq_xSample$Tumour_Type, 
#       function(x) round(percentage(x,0.001),digits=1))

# % < .05

#tapply(Mutations_Chisq_xSample[,5], Mutations_Chisq_xSample$Tumour_Type, 
#       function(x) round(percentage(x,0.05),digits=1))

# Let's add this info to "Mutations_Chisq_xPublish"

Mutations_Chisq_xPublish$CDS_PercSignSamples001 <- paste(tapply(Mutations_Chisq_xSample[,5], Mutations_Chisq_xSample$Tumour_Type, 
                                                             function(x) round(percentage(x,0.001),digits=1)), "%",sep="")
Mutations_Chisq_xPublish$Exon_PercSignSamples001 <- paste(tapply(Mutations_Chisq_xSample[,6], Mutations_Chisq_xSample$Tumour_Type, 
                                                            function(x) round(percentage(x,0.05),digits=1)), "%",sep="")

