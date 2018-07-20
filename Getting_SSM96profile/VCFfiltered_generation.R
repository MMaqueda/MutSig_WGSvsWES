
#Los archivos que queremos filtrar:
VCFfilesPath = "/project/devel/PCAWG/somatic_consensus_aug2016/data/filtered_calls/whitelist_12oct/ssm"
vcf_files <- list.files(path=VCFfilesPath,pattern = ".vcf.gz", full.names = TRUE) #Including also .tbi
#Excluding .tbi files
vcf_files <- vcf_files[-(grep(pattern=".tbi", x= vcf_files, fixed=TRUE))]

#Los archivos generados una vez filtrados:
file_names <- list.files(path=VCFfilesPath,pattern = ".vcf.gz", full.names = FALSE)
file_names <- file_names[-(grep(pattern= ".tbi", x=file_names, fixed=TRUE))]

setwd("/project/devel/PCAWG/mmaqueda/WES_VCFfiles/VCF_Filtered/") #Definimos el directorio donde guardar los archivos filtrados

file_namesCDS <- paste0("CDSfilt_",file_names,sep="")
file_namesExon <- paste0("EXONfilt_",file_names,sep="")
file_namesNext <- paste0("NEXTERAfilt_",file_names,sep="")

#Aquí consideramos que test.bed incluye los rangos a incluir (i.e. Exoma)
for (i in 2201:length(vcf_files))
{
  #Esto lo tenemos que hacer por el error de la entrada de INFO - temporalmente creamos un archivo sin esa info
  commandRemINFO <- paste0('zgrep -v "INFO=<ID=dbsnp_somatic,Number=.,Type=Flag" ',vcf_files[i], " | gzip -c > tmp.vcf.gz")
  system(commandRemINFO)

  commandFILT_CDS <- paste0("vcftools --gzvcf tmp.vcf.gz --bed /project/devel/PCAWG/mmaqueda/WES_VCFfiles/BED/gencode.v19.CDS_Merged_nochr.bed --recode --stdout | gzip -c > ./CDS/",
                    file_namesCDS[i], sep="")
  system(commandFILT_CDS)
  
  commandFILT_EXON <- paste0("vcftools --gzvcf tmp.vcf.gz --bed /project/devel/PCAWG/mmaqueda/WES_VCFfiles/BED/gencode.v19.Exon_Merged_nochr.bed --recode --stdout | gzip -c > ./EXON/",
                            file_namesExon[i], sep="")
  system(commandFILT_EXON)
  
  commandFILT_Nextera <- paste0("vcftools --gzvcf tmp.vcf.gz --bed /project/devel/PCAWG/mmaqueda/WES_VCFfiles/BED/Nexteramerged_nochrString.bed --recode --stdout | gzip -c > ./NEXTERAplatform/",
                             file_namesNext[i], sep="")
  system(commandFILT_Nextera)
  
  #Borramos el archivo que habíamos creado temporalmente
  system("rm tmp.vcf.gz")
}




# Comments / Issues -------------------------------------------------------
#It is not possible to apply VCFtools directly since following error emerges:
#Error: Flag Type must have 0 entries: ID=dbsnp_somatic,Number=.,Type=Flag,Description="Known-somatic dbSNP variant">

#VCFtools is expecting Number=0 since it is Type=Flag