
# This script is used to inspect the results from the "de novo" computation
# based on the WTSI Framework. An input is needed with the results from running
# "Compute_deNovo" script (based on mutSignatures package)
# REMARK: de novo results start with 2 signatures in any case

require(latticeExtra)
require(MutationalPatterns)
require(SomaticSignatures) 
#In case previous pacakge cannot be loaded, load the needed function from it:
#source("/project/devel/PCAWG/mmaqueda/Rscripts/normalizeMotifs.R")

# Custom scripts. Change the path if needed!!!!
source("/project/devel/PCAWG/mmaqueda/Rscripts/FrobeniusNorm.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_deNovo_signatures.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/plot_exposures_deNovo.R")

# Retrieve .RData with the results to inspect
# Example for Lung-AdenoCA in CDS region
load("/project/devel/PCAWG/mmaqueda/Results/deNovo_signatures/Lung-AdenoCA/CDS_Lung-AdenoCA_Results.RData")

# These RData only include one object (a list) named as "results"

deNovo <- results
tumor <- "Lung-AdenoCA"  
region <- "CDS"

# Plot Sign Reproducibility & Frob error vs #Signatures -------------------

silh <- unlist(lapply(deNovo, function(x) mean(x$processStabAvg)))
names(silh) <- NULL

rec_error  <- unlist(lapply(deNovo, function(x)
{
  if(length(x)==0) {"NA"}
  else {FrobeniusNorm(x$input$mutCounts, x$processes, x$exposures)}
}))

  
names(rec_error) <- NULL

toplot <- data.frame(Number_Signatures = as.factor(seq(2, length(deNovo)+1)), 
                     Reconstruction_Error = rec_error,
                     Signature_Reproducibility = silh)

# Figure with double Y axis....it can be misleading since both axis have different magnitudes

obj1 <- xyplot(Signature_Reproducibility ~ Number_Signatures, toplot, type = "b" , lwd=2,
               ylim=c(0,1.1))
obj2 <- xyplot(Reconstruction_Error ~ Number_Signatures, toplot, type = "b", lwd=2)

pdf("WGS_Silhoutte_vs_RecError.pdf")
doubleYScale(obj1, obj2, add.ylab2 = TRUE)
dev.off()

# Based on previous graph a k (number of signatures) MUST BE SELECTED
# Selection MUST BE DONE MANUALLY (BY VISUAL INSPECTION)

selected <- 3 #Here to specify the number of signatures selected or another to visualize
# To retrieve the specific results for that number of signatures we have to go to -1 element
# Reminder: the de novo computation always starts with k=2 signatures

# CDS signatures frequency correction  ------------------------------------
# This is ONLY NEEDED FOR CDS results. For WGS results, go to next section

# Load the correction factor from cds2wgs and exon2wgs
load("/project/devel/PCAWG/mmaqueda/Results/Mutations_and_Content_WGSvWES/Frequency_Corr_Factors_related/Correction_Factors.RData")

# In order to use the correction function we need to renamed the rownames as in the SomSign package

subs <-  c('CA','CG','CT','TA','TC','TG')
context <- c('A.A','A.C','A.G','A.T','C.A','C.C','C.G','C.T',
             'G.A','G.C','G.G','G.T','T.A','T.C','T.G','T.T')

ssm96_somsign <- paste(rep(subs,each=16), rep(context,times=4), sep=" ")

# Correction is applied to obtained processes
rownames(results[[selected-1]]$processes) <- ssm96_somsign
rownames(results[[selected-1]]$processesStd) <- ssm96_somsign

processes_corr <- normalizeMotifs(results[[selected-1]]$processes, corr_factorCDS2WGS)
processesStd_corr <- normalizeMotifs(results[[selected-1]]$processesStd, corr_factorCDS2WGS)


# Plot SSM96 Signatures profile -------------------------------------------

pdf("CDSK3_Signatures_Selected.pdf",height =8)
plot_deNovo_signatures(processes = deNovo[[selected-1]]$processes, 
                       processesSD = deNovo[[selected-1]]$processesStd, 
                       legend = TRUE, ymax=0.3, 
                       tumor = tumor, region = region)  
dev.off()


pdf("CDSk3_Signatures_Selected_CORR.pdf")
plot_deNovo_signatures(processes = processes_corr, 
                       processesSD = processesStd_corr, 
                       legend = TRUE, ymax=0.3, 
                       tumor = tumor, region = region)  
dev.off()


# Compare with COSMIC signatures ------------------------------------------
# Comparison is based on cosine similarity metric

# Load the 30 COSMIC signatures (matrix 96 x30). Object called 'cancer_signatures'
load("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/COSMICsign.RData")

CosSim_dNovCOSMIC <- as.data.frame(cos_sim_matrix(cancer_signatures, deNovo[[selected-1]]$processes))

# In case we have results from CDS, then use the corrected profile:
#CosSim_dNovCOSMIC <- as.data.frame(cos_sim_matrix(cancer_signatures, processes_corr))



# Plot Similarities de novo vs COSMIC -------------------------------------

# Adjust results for plotting
SIGNATURES <- paste0("Signature",LETTERS,sep="")
colnames(CosSim_dNovCOSMIC) <- SIGNATURES[1:selected]
CosSim_dNovCOSMIC$COSMIC <- rownames(CosSim_dNovCOSMIC)

toplot <- melt(CosSim_dNovCOSMIC, id.vars="COSMIC")

# Plot results
order <- CosSim_dNovCOSMIC$COSMIC
toplot$COSMIC <- factor(toplot$COSMIC, levels = order)


pdf("WGSk9_COSMIC_comparison_CosSim.pdf",width=19)

ggplot(data = toplot, aes(x = COSMIC, y = value)) + 
  geom_bar(stat = "identity", color="black", fill = "orange", size = 0.2)  + 
  facet_grid( ~ variable ) + ylab("Cosine Similiarity") + 
  coord_cartesian(ylim = c(0, 1.0)) + 
  scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
  guides(fill = FALSE) + theme_bw() + 
  theme(axis.title.y = element_text(size = 12, vjust = 1), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 6, angle = 90, 
                                   vjust = 0.4), strip.text.x = element_text(size = 9), 
        strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
        panel.spacing.x = unit(0, "lines")) +
  geom_hline(yintercept = 0.98, color="green", lty = "dashed") +
  annotate("text", x = 21, y = 1, label = "Extremely high Similarity",cex=3) +
  geom_hline(yintercept = 0.85, color="red", lty = "dashed") +
  annotate("text", x = 25, y = 0.87, label = "High Similarity",cex=3) +
  labs(title = paste0(tumor, " samples ", region, " region"))

dev.off()


# COSMIC comparison: Frob norm --------------------------------------------

# Let's compute the error (Frobenius norm) between each COSMIC signature and deNovo signatures

Error_dNovoCOSMIC <- as.data.frame(matrix(nrow = dim(cancer_signatures)[2], 
                                          ncol= selected +1))
colnames(Error_dNovoCOSMIC) <- c(SIGNATURES[1:selected],"COSMIC")
Error_dNovoCOSMIC$COSMIC <- colnames(cancer_signatures)

for(i in 1:selected)
{
  Error_dNovoCOSMIC[,i] <- apply(cancer_signatures, 2, function(COSMIC) sqrt(sum((COSMIC-deNovo[[selected-1]]$processes[,i])^2)))
  # In case we have results from CDS, then use the corrected profile:
  #Error_dNovoCOSMIC[,i] <- apply(cancer_signatures, 2, function(COSMIC) sqrt(sum((COSMIC-processes_corr[,i])^2)))
}

toplot <- melt(Error_dNovoCOSMIC, id.vars="COSMIC")

# Plot results
order <- Error_dNovoCOSMIC$COSMIC
toplot$COSMIC <- factor(toplot$COSMIC, levels = order)

pdf("WGSk9_COSMIC_comparison_FrobNorm.pdf",width=19)

ggplot(data = toplot, aes(x = COSMIC, y = value)) + 
  geom_bar(stat = "identity", color="black", fill = "blue", size = 0.2)  + 
  facet_grid( ~ variable ) + ylab("Error (Frobenius Norm)") + 
  coord_cartesian(ylim = c(0, 0.6)) + 
  scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
  guides(fill = FALSE) + theme_bw() + 
  theme(axis.title.y = element_text(size = 12, vjust = 1), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
        axis.text.x = element_text(size = 6, angle = 90, 
                                   vjust = 0.4), strip.text.x = element_text(size = 9), 
        strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
        panel.spacing.x = unit(0, "lines")) +
  geom_hline(yintercept = 0.098, color="red", lty = "dashed") +
  annotate("text", x = 21, y = 0.12, label = "Systematic Diff of 0.01",cex=3, color="red") +
  geom_hline(yintercept = 0.49, color="blue", lty = "dashed") +
  annotate("text", x = 21, y = 0.51, label = "Systematic Diff of 0.05",cex=3) +
  labs(title = paste0(tumor, " samples ", region, " region"))

dev.off()

# Be careful with prior plot since depending on the contribution values of the signature 
# a systematic error of 0.01 may be too much!!! An improvement to this graph would be to show %

# Plot Exposures in Stacked barplot ---------------------------------------------

pdf("WGSk9_Exposures2signatures.pdf")

plot_exposures_deNovo(deNovo[[selected-1]]$exposures, 
                      deNovo[[selected-1]]$input$sampleNames, 
                      coord_flip = FALSE,
                      tumor,region) 

dev.off()


# Optional: Compare signatures WGS vs CDS ---------------------------------

# This part is to compare the signatures obtained from WGS and CDS using cosine similarity metric

# The signatures obtained and selected from WGS and CDS (these CORRECTED) have to be merged 
# in the same matrix and name it as "processesALL" to obtain following figure

require(corrplot)
M <- as.data.frame(cos_sim_matrix(processesALL, processesALL))

pdf("SimMatSign.pdf")
corrplot(as.matrix(M), method="circle", is.corr=FALSE, addCoef.col = "black", 
          tl.col="black", tl.srt=45, tl.cex = 0.9,type="upper", col = brewer.pal(n = 10, name = "RdBu")) #Text label color and rotation,)
dev.off()
