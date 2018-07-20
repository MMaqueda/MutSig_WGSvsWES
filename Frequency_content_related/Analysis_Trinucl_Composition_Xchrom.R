
# This script is for analyzing the 3mer computations in WGS / Exon / CDS in the ref genome
# considering the different chromosomes (usually no chromosome differentiation is done).
# So, it will give an idea if the correction factor is missing important information

setwd("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/WGSvsWES_kmerCorrection/Analysis_kmerXchrom/")
load("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/WGSvsWES_kmerCorrection/WGSvsWES_base_kmer_counts.RData")

# Los datos que necesitamos de partida son: freq_3mer_CDS, freq_3mer_Exon and freq_3mer_WGS.
# Además se necesita roi_Exon y roi_CDS para poder converger en cromosomas estos casos

# Los nombres de freq_XXX son algo desafortunados ya que se refieren al número de veces, no
# tanto al porcentaje que al final es con lo que vamos a trabajar. 

# Converger counts_3mer_Exon y CDS a los distintos cromosomas (en vez de por regiones separadas)

converge_x_chrom <- function(counts_3mer_WES, roi_WES)
{
  data <- matrix(nrow = 24, ncol = 64)
  colnames(data) <- colnames(counts_3mer_WES)
  
  #Let's fill in this matrix with the info divided per chromosomes
  for(i in 1:dim(data)[1])
  {
    for(j in 1:dim(data)[2])
    {
      if(i==23) {k <- "X"}
      else if(i==24) k <- "Y"
      else k <- i
      data[i,j] <- sum(counts_3mer_WES[roi_WES@ranges@NAMES == k,j])
      }
  }
  return(data)
}

counts_3mer_Exon <- converge_x_chrom(freq_3mer_Exon, roi_Exon) 
counts_3mer_CDS <- converge_x_chrom(freq_3mer_CDS, roi_CDS) 

# Before saving, let's just keep the needed objects
require(gdata)
counts_3mer_WGS <- freq_3mer_WGS
keep(counts_3mer_WGS, counts_3mer_CDS, counts_3mer_Exon, sure=TRUE)
save.image("kmer_XChrom.RData")

# Transform counts matrix into relative contribution. Keeping chromosome diff
# No percentage is considered - If needed, just multiply the matrix per 100

counts2relative <- function(counts_3mer_XXX)
{
  rel_k3mer_WXX <- apply(counts_3mer_XXX,c(1,2), function(x) x/sum(as.numeric(counts_3mer_XXX)))
  return(rel_k3mer_WXX)
}

rel_3mer_WGS <- counts2relative(counts_3mer_WGS)
rel_3mer_Exon <- counts2relative(counts_3mer_Exon)
rel_3mer_CDS <- counts2relative(counts_3mer_CDS)

# Plot results
# First prepare the data to be plotted

toplot <- melt(as.data.frame(rel_3mer_WGS))
colnames(toplot) <- c("k3mer","Rel_WGS")

Rel_WESexon <- melt(as.data.frame(rel_3mer_Exon))
Rel_WEScds <- melt(as.data.frame(rel_3mer_CDS))
toplot <- cbind(toplot, Rel_Exon = Rel_WESexon$value, Rel_CDS = Rel_WEScds$value)

toplot <- melt(toplot)
colnames(toplot) <- c("k3mer","Region","Rel")

boxplot <- ggplot(toplot, aes(y=Rel, x=k3mer, fill=Region)) + 
  geom_boxplot() + 
  ggtitle("Relative counts of 3mer seq in reference genome considering WG, Exon or CDS regions") +
  stat_compare_means(aes(group = Region), label = "p.signif", hide.ns=TRUE) +  #Kruskal.Wallis test by default
  theme_bw() + theme(axis.text.x = element_text(angle=90,size=11,vjust=0.5),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) 

pdf("Comparison_kmerXchrom.pdf", width=20)
boxplot
dev.off()