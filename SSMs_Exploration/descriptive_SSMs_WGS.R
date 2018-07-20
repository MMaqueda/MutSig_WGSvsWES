
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(MutationalPatterns))

#Descriptive Statistics about the mutational profile of the samples 

#First we import the data we wanna analyze
#INPUT: a mutation matrix (SSM96) to be read from a text file. 

mut96_mat <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)  

##################
#Preparing the matrix information
##################

#Transposing the data matrix will be easier and convert it into a df
Tmut_mat <- as.data.frame(t(mut96_mat))

#Let's add a column with the tumour type (if it is a subtype, we indicate the subtype)
Tumour_Type <- sapply(rownames(Tmut_mat), function(x) unlist(strsplit(x, "[.]"))[2])
Organ <- sapply(Tumour_Type, function(x) unlist(strsplit(x, "-"))[1])

Tmut_mat$Tumour_Type <- Tumour_Type
Tmut_mat$Organ <- Organ

###################
#Table 1: Basic statistics per tumour type (Similar procedure could've done for Organ)
#Here we work with the total SSM - no base substitution differentiation
###################

TOT_ssm <- colSums(mut96_mat)

stats <- function(x) {
   c(min = min(x), max = max(x), mean = round(mean(x),digits=0),
     sd = round(sd(x),digits=0), 
     median = round(median(x),digits=0), 
     IQR = round(IQR(x),digits=0), 
     sum=sum(x))}

tab1 <- as.data.frame(
  cbind(do.call(rbind, tapply(TOT_ssm, Tmut_mat$Tumour_Type, stats))))

tab1 <- cbind(Organ = sapply(rownames(tab1), function(x) unlist(strsplit(x, "-"))[1]),
              N=plyr::count(Tmut_mat,'Tumour_Type')$freq,
              tab1)

pdf("./Exon_based/Plots/Summary_SSMs_pertumour_Exon.pdf",height=12,width=14)
grid.table(tab1)
dev.off()

###################
#Figure 0: Barplot with samples number of the whole dataset
##################

cols_pcawg <- pcawg.colour.palette(x = rownames(tab1),scheme = "tumour.subtype")
tab1$TumorType  <-  rownames(tab1)

#Let's re-order the cancer types based on the N
order <- tab1$TumorType[order(tab1$N)]
tab1$TumorType <- factor(tab1$TumorType, levels = order)

pdf("./Plots/Samples_x_Tumour.pdf")
ggplot(tab1, aes(y=N, x=TumorType, fill=TumorType)) + 
  geom_bar(stat="identity",show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = tab1$colors[order(tab1$N)]) +
  geom_text(aes(label=N), color="black", hjust=0, position=position_dodge(width=0.2)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black",size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()


###################
#Figure 1: Distribution of tumours in the dataset
##################

cols_pcawg <- pcawg.colour.palette(x = sort(unique(Tmut_mat$Tumour_Type)),scheme = "tumour.subtype")

#If we wanna do it per Organ, we should call "organ.system" scheme

slices <- tab1$N 
lbls <- sort(unique(Tmut_mat$Tumour_Type))
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 

pdf("TumourDistribution.pdf",width=14, height=14)
pie(slices,labels = lbls, col=cols_pcawg,
    main=paste0("Distribution of total ", sum(slices), " samples \n" ,"per Tumour Type"))
dev.off()

###################
#Figure 2: Distribution of SSMs per tumour type in the dataset
##################

slices <- tab1$sum #Total SSM per tumour type 
lbls <- sort(unique(Tmut_mat$Tumour_Type))
pct <- round(slices/sum(slices)*100,2)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 

#Let's reorder the slices and lbls per % of mutations
order <- order(slices,decreasing = TRUE)
slices <- slices[order]
lbls <- lbls[order]
cols_pcawg <- cols_pcawg[order]

pdf("./Plots/SSMsDistribution.pdf",width=14, height=14)
pie(slices,labels = lbls, col=cols_pcawg, clockwise = TRUE, 
    main=paste0("Distribution of total ", sum(slices), " SSMs \n" ,"per Tumour Type"))
dev.off()

##################
#Table 2: Basic statistics per tumour type (Similar procedure could've done for Organ)
#Here we work with the 6 base substitution 
##################

#First we have to collapse the 96 motif matrix into the 6 base substitutions. We will create a new df for this
#Although it is supposed that SSM96 has a specific order, let's avoid any problem due to different order

occurrences <- data.frame(row.names = rownames(Tmut_mat))
occurrences$Tumour_Type <- Tumour_Type
occurrences$`C>A` <- rowSums(Tmut_mat[,grep(pattern="C>A",x=colnames(Tmut_mat))])
occurrences$`C>G` <- rowSums(Tmut_mat[,grep(pattern="C>G",x=colnames(Tmut_mat))])
occurrences$`C>T` <- rowSums(Tmut_mat[,grep(pattern="C>T",x=colnames(Tmut_mat))])
occurrences$`T>A` <- rowSums(Tmut_mat[,grep(pattern="T>A",x=colnames(Tmut_mat))])
occurrences$`T>C` <- rowSums(Tmut_mat[,grep(pattern="T>C",x=colnames(Tmut_mat))])
occurrences$`T>G` <- rowSums(Tmut_mat[,grep(pattern="T>G",x=colnames(Tmut_mat))])

#And now we have to collapse the base substitution per tumour type
tab2 <- aggregate(. ~ Tumour_Type, data=occurrences, FUN=sum)
tab2 <- cbind(Tumour_Type = tab2$Tumour_Type,
              Organ = sapply(tab2$Tumour_Type, function(x) unlist(strsplit(x, "-"))[1]),
              tab2[,2:dim(tab2)[2]])

rownames(tab2) <- NULL

pdf("Occur_per_tumour.pdf",height=18,width=16)
grid.table(tab2)
dev.off()

##################
#Figure 3 and Figure 4:  
#Plotting mutational espectrum (6 bases substitutions) in general and per tumour type
##################

pdf("Occurrences_SSM6.pdf",width=10, height=12)
#General
plot_spectrum(occurrences[ ,2:7]) #From MutationalPatterns
#Per tumour type
plot_spectrum(occurrences[ ,2:7], by = Tumour_Type, legend = TRUE)
dev.off()

##################
#Figure 5: Plotting mutational profile (96 SSMs) as an example for some samples
##################

#Which samples to show? We can choose five sample randomly. Just to show the typical profile

samples <- sample(1:dim(Tmut_mat)[1],size=5)

pdf("Example_ProfileSSM96.pdf",width=10, height=12)
plot_96_profile(mut96_mat[,samples], condensed = TRUE)
dev.off()


##################
#Figure 8: Plot showing total SSMs per sample distributed per tumour type 
##################

#Gathering the needed data
toplot <- data.frame(Total_SSMs = TOT_ssm, Tumour_Type = Tumour_Type, stringsAsFactors = F)

order <- rownames(tab1[order(tab1$median),])
toplot$Tumour_Type <- factor(toplot$Tumour_Type, levels = order)

cols_pcawg <- pcawg.colour.palette(x = levels(toplot$Tumour_Type),scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

dotplot <- ggplot(toplot, aes(y=log10(Total_SSMs), x=Tumour_Type)) + 
  geom_jitter(width=.2, height=0.2) +
  ggtitle("SSM96 counts per tumour type \n CDS REGION") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) 

boxplot <- ggplot(toplot, aes(y=Total_SSMs, x=Tumour_Type, fill=Tumour_Type)) + 
  geom_boxplot() +
  scale_fill_manual(values= cols_pcawg) +
  scale_y_continuous(trans = "log10") + 
  ggtitle("SSM96 counts per tumour type \n CDS REGION") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12,vjust=0.5),
                     axis.text.y = element_text(size=18),
                     plot.title = element_text(face = "bold", size = (15))) 
  
pdf("CDS_BoxPlot_Total_SSMs_Tumour.pdf",width=14, height=14)
#dotplot
boxplot
dev.off()