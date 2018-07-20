

# Paper WES vs WGS analysis Belkadi PNAS ----------------------------------

#Paper Belkadi A., PNAS 2015
https://github.com/HGID/WES_vs_WGS/blob/master/Analysis_WES_WGS.R


# VCF generation with WES data --------------------------------------------
#Para generar los nuevos archivos VCF con la info WGS
#http://vcftools.sourceforge.net/man_latest.html

#Los archivos vcf.gz los puedo leer del directorio donde estan
#Y guardar los nuevos en el directorio de trabajo
#El comando es el siguiente, se tiene que parametrizar para generarlos en loop

#Antes de abrir R, hay que cargar  module load VCFtools

#Aquí consideramos que test.bed incluye los rangos a incluir (i.e. Exoma)
system("vcftools --gzvcf ./WES_VCFfiles/testsample.vcf.gz --bed test.bed --recode --stdout | gzip -c > filtered.vcf.gz")

#Si en vez de --bed utilizamos --exclude-bed, haremos lo contrario WGS -  exoma- Para el random/control!!
--exclude-bed



# Definicion del BED file de filtrado -------------------------------------
#Acerca del BED file WES. Hay diferentes fuentes de información. Algunos ejemplos:

# BED - from UCSC genome browser ------------------------------------------
# Esta opción NO la utilizaremos
#Info desde UCSC genome browser - database annotation (Downloads > Human....)
wget --timestamping 
'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/wgEncodeGencodeExonSupportV19.txt.gz ' -O WES.txt.gz

# BED - from Illumina kits ------------------------------------------------
#Los archivos BED (tienen el mismo número de entries, el de Nextera está OK) - TBD si se utilizará
#Info desde Library Kits de Illumina
https://support.illumina.com/sequencing/sequencing_kits/nextera-rapid-capture-exome-kit/downloads.html
https://support.illumina.com/sequencing/sequencing_kits/truseq-exome-kit/downloads.html

# BED - from GENCODE annot ------------------------------------------------
#Info de GENCODE directamente. Esta es la opción SELECCIONADA

#Servidor desde donde bajar el archivo de anotaciones de Gencode v19 en formato .gtf
http://ftp.sanger.ac.uk/pub/gencode/GencodeHuman/release_19/gencode.v19.annotation.gtf
tail -n +X gencode.v19.annotation.gtf > gencode.v19.annotation_noH.gtf
#where X has to be a number indicating the line from where to start

setwd("/Users/mariamaqueda/Documents/TFM/WESfromWGS/GENCODE_annot/")

# CDS regions -------------------------------------------------------------
#HAVANA and ENSEMBL annotations are both kept 
#Take the CDS regions. Filtramos en command line para para reducir tamaño antes de importar
system("awk '{ if ($3 == "CDS") {print $0}}' gencode.v19.annotation_noH.gtf > gencode.v19.CDS.gtf")

#Evaluemos qué tipos de gene_type y trancript type hay en esta info
#La columna 8 es el frame
gencodeCDS <- read.table(file = "gencode.v19.CDS.gtf", header= FALSE, sep="\t",stringsAsFactors = FALSE) #Info columna 9
info_entry <- as.data.frame(gencodeCDS$V9,stringAsFactors=FALSE)

gene_types <- apply(info_entry,1, function(k) unlist(strsplit(unlist(strsplit(k,split=";"))[3], split=" "))[3])
transcript_type <- apply(info_entry,1, function(k) unlist(strsplit(unlist(strsplit(k,split=";"))[6], split=" "))[3])
#Las definiciones de los diferentes tipos encontrados https://www.gencodegenes.org/gencode_biotypes.html

#Ahora creamos el BED file en command line: solo tomaremos el chr, start y end. El resto de info no la mantenemos
awk 'BEGIN{OFS="\t"}{print $1,$4,$5}' gencode.v19.CDS.gtf > tmp.gtf && mv tmp.gtf gencode.v19.CDS.gtf
head gencode.v19.CDS.gtf
# chr1	69091	70005
# chr1	138533	139309
# chr1	367659	368594
# chr1	621099	622034
# chr1	739121	739137

#Faltan eliminar las anotaciones de chrM (esto habría que hacerlo al inicio)
awk 'BEGIN{OFS="\t"}{if ($1 != "chrM") {print $0}}' gencode.v19.CDS.gtf > tmp.gtf && mv tmp.gtf gencode.v19.CDS.gtf

#Realizamos un merge de coordenadas del archivo .bed para contar número de nucleótidos como CDS
#Mediante BEDtools. Sin distinguir entre strands, mínimo un 1bp de overlap: flattened interval
sort -k1,1 -k2,2n gencode.v19.CDS.gtf | uniq > gencode.v19.CDSsorted_uniq.gtf
#Primero se ordenan datos (requisito para mergeBed) - Hay muchos duplicados!!

#Convertimos los .gtf (1-based) to .bed (0-based) with awk (ver notas al final)
awk 'BEGIN{OFS="\t"} {print $1,$2-1,$3}' gencode.v19.CDSsorted_uniq.gtf > gencode.v19.CDSsorted_uniq.bed

#Ahora ya podemos hacer el merge en el cluster (BEDtools there)
#Info de cómo hacer el merge http://bedtools.readthedocs.io/en/latest/content/tools/merge.html
#BEDTools/2.16.2(default)
module load intel/16.3.067
module load BEDTools
  

#Por último, generamos un nuevo archivo donde no aparezca "chr", a VCFtools no le gusta...
awk 'BEGIN{OFS="\t"}{gsub(/chr/,"",$1); print $0}' gencode.v19.CDSsorted_uniq_merged.bed > gencode.v19.CDSnochr.bed

#Conversión de BED a GRanges. Tomamos los datos del Merge (será lo mismo que el no merge para este caso)
#Tiene que ser sin chr al indicar el cromosoma
require(GenomicRanges)
CDSmerged_bed <- read.table("./Final_BED_files/gencode.v19.CDS_Merged_nochr.bed",
                        header=FALSE, sep="\t")

# Start+1 ya que BED es 0-based
CDSmerged_GRanges <- GRanges(seqnames = CDSmerged_bed$V1, 
                             ranges = IRanges(start = CDSmerged_bed$V2 +1, end = CDSmerged_bed$V3))

# Exon regions ------------------------------------------------------------
#THIS CODE IS EXACTLY THE SAME AS FOR CDS REGION - it's copy paste basically

#HAVANA and ENSEMBL annotations are both kept 
#Take the EXON regions. Filtramos en command line para para reducir tamaño antes de importar
system("awk '{ if ($3 == "exon") {print $0}}' gencode.v19.annotation_noH.gtf > gencode.v19.exon.gtf")

#Evaluemos qué tipos de gene_type y trancript type hay en esta info
#La columna 8 es el frame
gencodeExon <- read.table(file = "gencode.v19.exon.gtf", header= FALSE, sep="\t",stringsAsFactors = FALSE) #Info columna 9
info_entry <- as.data.frame(gencodeExon$V9,stringAsFactors=FALSE)

gene_types <- apply(info_entry,1, function(k) unlist(strsplit(unlist(strsplit(k,split=";"))[3], split=" "))[3])
transcript_type <- apply(info_entry,1, function(k) unlist(strsplit(unlist(strsplit(k,split=";"))[6], split=" "))[3])
#Las definiciones de los diferentes tipos encontrados https://www.gencodegenes.org/gencode_biotypes.html

#Ahora creamos el BED file en command line
awk 'BEGIN{OFS="\t"}{print $1,$4,$5}' gencode.v19.exon.gtf > tmp.gtf && mv tmp.gtf gencode.v19.exon.gtf
head gencode.v19.exon.gtf
# chr1	11869	12227
# chr1	12613	12721
# chr1	13221	14409
# chr1	11872	12227
# chr1	12613	12721
# chr1	13225	14412

#Faltan eliminar las anotaciones de chrM (esto habría que hacerlo al inicio)
awk 'BEGIN{OFS="\t"}{if ($1 != "chrM") {print $0}}' gencode.v19.exon.gtf > tmp.gtf && mv tmp.gtf gencode.v19.exon.gtf

#Realizamos un merge de coordenadas del archivo .bed para contar número de nucleótidos como CDS
#Mediante BEDtools. Sin distinguir entre strands, mínimo un 1bp de overlap: flattened interval
sort -k1,1 -k2,2n gencode.v19.exon.gtf | uniq > gencode.v19.Exonsorted_uniq.gtf
#Primero se ordenan datos (requisito para mergeBed) - Se eliminan duplicados

#Convertimos los .gtf (1-based) to .bed (0-based) with awk (ver notas al final)
awk 'BEGIN{OFS="\t"} {print $1,$2-1,$3}' gencode.v19.Exonsorted_uniq.gtf > gencode.v19.Exonsorted_uniq.bed

#Ahora ya podemos hacer el merge
bedtools merge -i gencode.v19.Exonsorted_uniq.bed > gencode.v19.Exon_Merged.bed  

#Por último, generamos un nuevo archivo donde no aparezca "chr"
awk 'BEGIN{OFS="\t"}{gsub(/chr/,"",$1); print $0}' gencode.v19.Exon_Merged.bed > gencode.v19.Exon_Merged_nochr.bed

#Conversión de BED a GRanges
require(GenomicRanges)
Exonmerged_bed <- read.table("./Final_BED_files/gencode.v19.Exon_Merged_nochr.bed",
                            header=FALSE, sep="\t")

Exonmerged_GRanges <- GRanges(seqnames = Exonmerged_bed$V1, 
                             ranges = IRanges(start = Exonmerged_bed$V2 +1 , end = Exonmerged_bed$V3))


# Other: Nextera annotations ----------------------------------------------
Nextera_bed <- read.table("/Users/mariamaqueda/Documents/TFM/WESfromWGS/Annotations_Info/WES_platforms/nexterarapidcapture_exome_targetedregions_hg19.bed",
                             header=FALSE, sep="\t")
#Eliminamos el chr en la primera columna para generar el GRanges
Nextera_bed$V1 <- sapply(as.character(Nextera_bed$V1), function(chr) unlist(strsplit(x=chr,split="chr"))[2])

Nextera_GRanges <- GRanges(seqnames = Nextera_bed$V1, 
                           ranges = IRanges(start = Nextera_bed$V2 +1, end = Nextera_bed$V3))

# Other: Agilent SureSelect v5 annotations --------------------------------

Agilentv5_bed <- read.table("/project/devel/PCAWG/mmaqueda/WES_VCFfiles/BED/Agilent_SureSelect.v5.bed",
                          header=FALSE, sep="\t")

# There's no 'chr' to be removed from first column in this case
Agilentv5_GRanges <- GRanges(seqnames = Agilentv5_bed$V1, 
                           ranges = IRanges(start = Agilentv5_bed$V2 +1, end = Agilentv5_bed$V3))

# BED and GTF formats -----------------------------------------------------

#GTF files are 1-based and inclusive on both sides of the interval
#BED is 0-based and non-inclusive (or exclusive) on the right
