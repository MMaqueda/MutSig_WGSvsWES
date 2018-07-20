#' Plot the mean and sd values for Signatures contribution obtained from 
#' a samples data set
#'
#' @param SignContribution_summary    a dataframe with the mean and sd 
#'                                    values per signatures. Typically the output of SignContrSummary
#' @param whichtumour    Character string indicating from which tumour to include in graph title
#' @param data    character string indicating data origin (i.e. WGS) to include in graph title
#' @param norm    Boolean. Indicate if dataframe must be normalized (TRUE)  - if contribution is 
#'                in counts instead of fraction (i.e. output fit_to_signatures) or not (FALSE) 
#'                (i.e. output of whichSignatures)  
#' @return          Plot
#' 

plot_summary_SignaturesContribution <- function(SignContribution, whichtumour, norm, data){
  
  # Normalize data if it is in counts
  if(norm==TRUE)
  {SignContribution <- sweep(as.matrix(SignContribution), 2, colSums(as.matrix(SignContribution)), `/`)}
  else SignContribution <- SignContribution
  
  # Set up all the colors possible
  all.sigs <- c(paste(rep("Signature.",30), seq(1,30,1), sep=""))              
  all.colors  <- c("#023FA5","#BEC1D4","#D6BCC0","#BB7784","gold","#4A6FE3","#8595E1","#B5BBE3",
                   "#E6AFB9","#E07B91","#D33F6A","#11C638","#8DD593","#C6DEC7","#EAD3C6","#F0B98D",
                   "#EF9708","#0FCFC0","#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4","#866097",
                   "#008941","#A30059","#008080","#8B0000","#F4A460","#663399")

  # Set up values for plotting
  toplot <- melt(as.data.frame(SignContribution))
  colnames(toplot) <- c("sample","SignContribution")
  toplot$Signature <- rep(paste(rep("Signature.",30), seq(1,30,1), sep=""),dim(SignContribution)[2])
  #Number of samples where signature has any contribution
  toplot$NumSamples <- apply(SignContribution, 1, function(row) length(which(row!=0)))
  toplot$Color <- rep(all.colors, dim(SignContribution)[2])
    
  #Boxplot for Signatures Contribution
  toplot$Signature <- factor(toplot$Signature, levels = all.sigs) # to keep order
  
  pval <- sapply(all.sigs, function(x) {
    test <- t.test(toplot$SignContribution[toplot$Signature == x], mu=0.05, alternative = "greater")
    return(test$p.value)})
  
  pval <- round(pval,3)
  pval <- rep(pval,dim(SignContribution)[2])
  pval[pval>0.05]<- "ns"
  pval[pval<0.001]<- "<0.001"
  
  boxplot <- ggplot(toplot, aes(y=SignContribution, x=Signature, fill=Color),show.legend=F) + 
    geom_boxplot(show.legend=F) +
    ylim(0,1) +
    ggtitle(paste("Signature Contribution in",whichtumour, " dataset \n", data, "data (t.test > 0.05")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, size=10,vjust=0.5),
                       axis.text.y = element_text(size=14),
                       plot.title = element_text(face = "bold", size = (15)))+
    annotate("text", x = toplot$Signature, y = 1, 
             label = paste(toplot$NumSamples,"/",dim(SignContribution)[2]), angle=90)+
    annotate("text", x=toplot$Signature, y =0.85, label=paste("p=",pval), angle=90) +
    geom_hline(yintercept = 0.05, color="red", lty = "dashed") 
  
  return(boxplot)
}