
#' This script is for decomposing COSMIC signatures: We are going to decompose each of the 30 signatures 
#' using the rest 29. Decomposition thru method NNLS (MutationalPatterns) and heuristic (deconstructSigs)

require(MutationalPatterns)
require(deconstructSigs)

source("/project/devel/PCAWG/mmaqueda/Rscripts/FrobeniusNorm.R")

# In a prior step, the 30 COSMIC signatures were already obtained and prepared (matrix 96 x30).
# If wanted to be repeated, follow indications from MutationalPatterns package.

# The needed object is called 'cancer_signatures' and is stored in:
load("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/COSMICsign.RData")

# Method ONE -NNLS-: fit_to_signatures function from MutationalPatterns ----------

COSMICdecomp <- as.data.frame(matrix(nrow=30, ncol=31))
rownames(COSMICdecomp) <- colnames(cancer_signatures)
colnames(COSMICdecomp) <- c(colnames(cancer_signatures), "Rec error")

for(i in 1:dim(cancer_signatures)[2]){
  
  fit <- fit_to_signatures(as.data.frame(cancer_signatures[,i]), cancer_signatures[,-i])
  
  # Errors 
  error <- FrobeniusNorm(as.data.frame(cancer_signatures[,i]), 
                         as.matrix(cancer_signatures[,-i]),
                         fit$contribution) 
  
  #Save results in the dataframe
  COSMICdecomp[i,c(1:30)[-which(c(1:30) == i)]] <- fit$contribution
  COSMICdecomp[i,31] <- error
  
  # Plot decomposition: Barplots of the decomposition of each signatures showing the 
  # contribution from the rest. A pdf will be generated for each signature
  
  # toplot <- as.data.frame(fit$contribution)
  # colnames(toplot) <- "Signature"
  # toplot$signatures <- factor(rownames(toplot), levels=rownames(toplot))
  
  # p <- ggplot(toplot, aes(x=signatures,y=Signature),show.legend=F) + 
  #   geom_bar(stat="identity",fill="steelblue") +
  #   coord_cartesian(ylim = c(0.0, 1.0)) +
  #   scale_y_continuous(breaks = seq(0, 1, 0.1)) + guides(fill = FALSE) + 
  #   theme_bw() + 
  #   theme(axis.text.x = element_text(angle = 90, size=10,vjust=0.5),
  #         axis.text.y = element_text(size=12),
  #         plot.title = element_text(face = "bold", size = (15)),
  #         axis.title.x=element_blank()) +
  #   labs(title=paste0("Signature ", i, " Decomposition \n Reconstruction Error: ", round(error,3)),
  #        y="Sig. contribution")
  # 
  # pdf(paste0("M1_Signature",i,".pdf"),height=8, width=8)
  # print(p)
  # dev.off()
}

# Method TWO - heuristic-: whichSignatures function from deconstructSigs ----------

for(i in 1:dim(cancer_signatures)[2]){
  
  sample <- as.data.frame(t(cancer_signatures[,i]))
  rownames(sample) <- paste0("Signature",i)
  
  fit <- whichSignatures(tumor.ref = sample, #Expected to have samples in rows 
                         signatures.ref = as.data.frame(t(cancer_signatures[,-i])), #Signatures in rows in case of inclusion
                         associated = c(), #If we wanna narrow the list of signatures
                         signatures.limit = NA, 
                         signature.cutoff = 0, #Weight less than this amount (this value is the one by default)
                         contexts.needed = FALSE, #Since the data input is in counts => This is for normalizing
                         tri.counts.method = "default")
  
  # Errors 
  error <- FrobeniusNorm(as.data.frame(cancer_signatures[,i]), 
                         as.matrix(cancer_signatures[,-i]),
                         as.matrix(t(fit$weights))) 
  
  # Plot decomposition: Barplots of the decomposition of each signatures showing the 
  # contribution from the rest. A pdf will be generated for each signature
  
  # toplot <- as.data.frame(t(fit$weights))
  # colnames(toplot) <- "Signature"
  # toplot$signatures <- factor(rownames(toplot), levels=rownames(toplot))
  # 
  # p <- ggplot(toplot, aes(x=signatures,y=Signature),show.legend=F) + 
  #   geom_bar(stat="identity",fill="steelblue") +
  #   coord_cartesian(ylim = c(0.0, 1.0)) +
  #   scale_y_continuous(breaks = seq(0, 1, 0.1)) + guides(fill = FALSE) + 
  #   theme_bw() + 
  #   theme(axis.text.x = element_text(angle = 90, size=10,vjust=0.5),
  #         axis.text.y = element_text(size=12),
  #         plot.title = element_text(face = "bold", size = (15)),
  #         axis.title.x=element_blank()) +
  #   labs(title=paste0("Signature ", i, " Decomposition \n Reconstruction Error: ", round(error,3)),
  #        y="Sig. contribution")
  # 
  # pdf(paste0("M2_Signature",i,".pdf"),height=8, width=8)
  # print(p)
  # dev.off()
}

# Optionally, results could be prepared for publishing. Here, an example with those from NNLS method:

COSMICdecompM1 <- apply(COSMICdecomp, c(1,2), function(x) {if(is.na(x)) {x<-NA}else if(x<0.01){x<-"-"} else x<-round(x,2)})

