
###################
#Compute number of bases really sequenced WG (not considering N)
#Compute the same figures for RoI: CDS, Exon, Nextera
###################

#Info from Miranda

library(Biostrings)
library(Rsamtools)

fastafileGenome <- "/project/devel/PCAWG/hsapiens_v37.unmasked.fa"

#Count all four bases, 
genome_dnastringset <- readDNAStringSet(fastafileGenome)
genome_dnastringset <- genome_dnastringset[names(genome_dnastringset) %in% c(1:22, "X", "Y")]
freq_bases_WGS <- oligonucleotideFrequency(genome_dnastringset, width=1, step=1)  #Counts per chromosome

#head(freq_bases_WGS) 
#             A        C        G        T
# [1,] 38330752 27308648 27298423 38376915
# [2,] 38307244 27236798 27268038 38317436
# [3,] 38604831 26634995 26617050 38624517
# [4,] 29336945 18412698 18414776 29425459
# [5,] 25992966 18027132 18071947 26197495
# [6,] 23620876 17247582 17228387 23597921

sum(as.numeric(freq_bases_WGS))
#[1] 2861327131

#Count frequency of trinucleotides:
freq_3mer_WGS <- trinucleotideFrequency(genome_dnastringset, step=1) #width=3 by default (tri)

#head(freq_3mer_WGS)  #Son 64 combinacions (4*4*4)
#         AAA     AAC     AAG     AAT     ACA     ACC    ACG     ACT     AGA
# [1,] 4913799 1890136 2602308 3170822 2636754 1552548 346029 2090946 2886297
# [2,] 4841424 1881914 2628274 3150674 2613457 1546804 323065 2099589 2925471
# [3,] 5043050 1906014 2610526 3253179 2609221 1521148 323795 2110273 2894755
# [4,] 3973731 1441561 1936258 2672353 1987774 1023386 223100 1565251 2112716
# [5,] 3365391 1276604 1756924 2189082 1757882 1024632 216305 1421610 1943613
# [6,] 3018731 1169386 1617890 1888125 1624775  989469 216514 1309962 1794130

#This will allow you to get a specific region of the fasta
#open fasta file: In fact, there's no need to open it
#genome_fasta <- open(FaFile(fastafileGenome))

# specific region (or regions) of interest
# specific_range_of_genome_sequence <- GRanges(seqnames = <chromosome>, ranges=IRanges(start=<starting_position>, end=<end_position>)
# In our case, the GRanges are already computed and stored in a RData
load("/project/devel/PCAWG/mmaqueda/Data/WGSvsWES_kmerCorrection/ROI_GRanges.RData")

# get the actual sequence of the region of interest
roi_CDS <- scanFa(fastafileGenome, param=CDSmerged_GRanges)
roi_Exon <- scanFa(fastafileGenome, param=Exonmerged_GRanges)
roi_Nextera <- scanFa(fastafileGenome, param=Nextera_GRanges)
roi_Agilentv5 <- scanFa(fastafileGenome, param=Agilentv5_GRanges)
                                            
# Counting the bases in the sequence you retrieved
freq_bases_CDS <- oligonucleotideFrequency(roi_CDS, width=1, step=1) 
freq_bases_Exon <- oligonucleotideFrequency(roi_Exon, width=1, step=1) 
freq_bases_Nextera <- oligonucleotideFrequency(roi_Nextera, width=1, step=1) 
freq_bases_Agilentv5 <- oligonucleotideFrequency(roi_Agilentv5, width=1, step=1) 

# Counting the 3-mer in those sequences
freq_3mer_CDS <- trinucleotideFrequency(roi_CDS, step=1)
freq_3mer_Exon <- trinucleotideFrequency(roi_Exon, step=1)
freq_3mer_Nextera <- trinucleotideFrequency(roi_Nextera, step=1)


sum(as.numeric(freq_bases_CDS))
#[1] 35334619
sum(as.numeric(freq_bases_Exon))
#[1] 121985209
sum(as.numeric(freq_bases_Nextera))
#[1] 37105383
sum(as.numeric(freq_bases_Nextera))/sum(as.numeric(freq_bases_WGS))
#[1] 0.01296789
sum(as.numeric(freq_bases_Exon))/sum(as.numeric(freq_bases_WGS))
#[1]0.04263239
sum(as.numeric(freq_bases_CDS))/sum(as.numeric(freq_bases_WGS))
#[1] 0.01234903

##############################
#Los probs instead of counts...
prob_bases_WGS <- oligonucleotideFrequency(genome_dnastringset, width=1, step=1,as.prob=TRUE) 
prob_bases_CDS <- oligonucleotideFrequency(roi_CDS, width=1, step=1,as.prob=TRUE) 
prob_bases_Exon <- oligonucleotideFrequency(roi_Exon, width=1, step=1,as.prob=TRUE) 
prob_bases_Nextera <- oligonucleotideFrequency(roi_Nextera, width=1, step=1,as.prob=TRUE) 
prob_3mer_WGS <- trinucleotideFrequency(genome_dnastringset, step=1, as.prob=TRUE) 
prob_3mer_CDS <- trinucleotideFrequency(roi_CDS, step=1, as.prob=TRUE) 
prob_3mer_Exon <- trinucleotideFrequency(roi_Exon, step=1, as.prob=TRUE) 
prob_3mer_Nextera <- trinucleotideFrequency(roi_Nextera, step=1, as.prob=TRUE) 

rm(genome_dnastringset)
save.image("./WGSvsWES_base_kmer_counts.RData")

