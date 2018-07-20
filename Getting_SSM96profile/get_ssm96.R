
Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(parallel))

#' This function reads VCF compressed files (.vcf.gz) and obtains the 96 mutational matrix
#' @param resultsPath      Directory where the results wants to be stored
#' @param VCFfilesPath     Directory where the VCF files are available
#' @param Info_VCFs.RData  RData file (full path) where the information about the VCF files is  
#' @param ref_genome       A character string specifying the reference genome
#' @param whichtumour      Tumour type to call for the mutation matrix. One tumour at a time!
#' @param region           Genome region used. Relevant since file names include a string with the region
#'                         region = c("CDS","Exon","Nextera","WGS")
#' @return                 There is no return. Result is stored in a txt file in resultsPath
#'                         Matrix colnames refer to SampleID.histology_abbreviation as in RData    
#' 
#' Intern loop is using multicore (mc) function. No cores are specified, it will automatically use
#' the cores available at that moment.
#' 
#' Example of usage:
#' get_ssm96(resultsPath = "/project/devel/PCAWG/mmaqueda/SSM96Matrix/",
#'     VCFfilesPath = "/project/devel/PCAWG/somatic_consensus_aug2016/data/filtered_calls/whitelist_12oct/ssm",
#'     InfoVCFs.RData = "/project/devel/PCAWG/somatic_consensus_aug2016/data/t_sample2tumortype_aliquot.RData", 
#'     ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     whichtumour = "Breast-DCIS",
#'     region = "CDS")
#' 
#' IMPORTANT REMARKS. It is considered that in VCFfilesPath there are .vcf.gz and .tbi files. No issue if 
#' there are no .tbi files in the directory.   
#' It is assumed that all files include ".consensus.20160830.somatic.snv_mnv" in their name.

get_ssm96 <- function(resultsPath, VCFfilesPath, InfoVCFs.RData, ref_genome,whichtumour,region)
{
  start.time <- Sys.time()
  cat("Info: Starting -get_ssm96- script at ", strftime(start.time), " ...\n", sep = "")
  
  #Load the aliquot IDs (also files id) linked to the tumor type
  infoVCF <- get(load(InfoVCFs.RData))
  vcf_files <- list.files(path=VCFfilesPath,pattern = ".vcf.gz", full.names = TRUE) #Including also .tbi
  
  #Excluding .tbi files - Only need when explored region is WGS
  if(region == "WGS")
    vcf_files <- vcf_files[-(grep(pattern=".tbi",x= vcf_files, fixed=TRUE))]
  
  #Match the order of the VCF file names with the df info about those VCF
  #First, get the file names with just the aliquot id
  file_names <- list.files(path=VCFfilesPath,pattern = ".vcf.gz", full.names = FALSE)
  
  if(region == "WGS")
    file_names <- file_names[-(grep(pattern= ".tbi", x=file_names, fixed=TRUE))]
  
  file_names <- unlist(strsplit(x = file_names, split=".vcf.gz",fixed=TRUE))
  file_names <- unlist(strsplit(x = file_names, 
                                split=".consensus.20160830.somatic.snv_mnv",fixed=TRUE))
  
  #Consider different prefix depending on the genome region (VCF source)
  if(region == "CDS")
    file_names <- sapply(file_names, function(name) unlist(strsplit(x=name,split="CDSfilt_"))[2])
  else if(region == "Exon")
    file_names <- sapply(file_names, function(name) unlist(strsplit(x=name,split="EXONfilt_"))[2])
  else if(region== "Nextera")
    file_names <- sapply(file_names, function(name) unlist(strsplit(x=name,split="NEXTERAfilt_"))[2])
  
  #Now we can match vcf files with names in the dataframe(metadata)
  new_order <- match(infoVCF$tumor_wgs_aliquot_id,file_names)  
  vcf_files <- vcf_files[new_order]
  
  sample_names <- paste0(infoVCF$t_sample_index,".",
                         infoVCF$histology_abbreviation)
  
  #Select only those corresponding to an specific tumour type
  index <- grep(x=sample_names,pattern=whichtumour)

  #From VCF to GRanges object and from there to mutational profile
  ssm.set <- mcmapply(function(vcf_file, sample_name, samples) {
    tryCatch({
      cat("Reading sample...", samples, "\n")
      cat("Info: Loading VCF file \"", basename(vcf_file), "\" ...\n", sep = "")
      
      vcfs <- read_vcfs_as_granges(vcf_files = vcf_file, 
                                   sample_names = sample_name,
                                   genome = ref_genome)
      
      mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
      cat("Info: 96 Mutational profile obtained \n")
      return(mut_mat)
    }, error = function(e) {
    stop("Failed to load this VCF file. Stop", call. = F)
    })
  },
  vcf_file = vcf_files[index], #vcf_files
  sample_name = sample_names[index],   #sample_names
  samples = index, 
  SIMPLIFY = F)
  
  #Wrap-up results: 96 mutational matrix. Save it in a txt file
  mut96_mat <- as.data.frame(ssm.set, check.names=FALSE)

  write.table(x = mut96_mat,
              file = paste0(resultsPath,"mut96_mat_",
                            whichtumour,".txt",sep=""),
              sep = "\t")
}




