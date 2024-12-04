#### Install ggRibo and its required packages:
(1) Install required packages.
```
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install Bioconductor packages
bioconductor_packages <- c(
    "GenomicRanges", 
    "GenomicFeatures", 
    "GenomicAlignments", 
    "Rsamtools", 
    "IRanges", 
    "BiocParallel", 
    "txdbmaker"
)

BiocManager::install(bioconductor_packages)

# Install CRAN packages
cran_packages <- c(
    "ggplot2", 
    "cowplot", 
    "grid", 
    "dplyr", 
    "R6"
)

install.packages(setdiff(cran_packages, rownames(installed.packages())))
```
(2) Install RiboPlotR.
```
#Install RiboPlotR
library(devtools)
install_github("hsinyenwu/ggRibo")
```
### Basic usage of ggRibo
#### Load RNA-seq, Ribo-seq and annotation files 
```
#Path for example data from ggRibo package
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE) #Annotation
ugtf <- system.file("extdata", "AT3G02468.gtf", package = "ggRibo", mustWork = TRUE) #uORF gtf
Root_RNA <- system.file("extdata", "Root_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNA <- system.file("extdata", "Shoot_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
Root_RNAse <- system.file("extdata", "Root_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNAse <- system.file("extdata", "Shoot_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
Root_Ribo <- system.file("extdata", "riboRoot.bed", package = "ggRibo", mustWork = TRUE) #Root Ribo-seq data
Shoot_Ribo <- system.file("extdata", "riboShoot.bed", package = "ggRibo", mustWork = TRUE) #Shoot Ribo-seq data

#Define sample names  
Samples=c("Root","Shoot")
#Load Ribo-seq data
RiboseqData=Ribo_data(c(Root_Ribo,Shoot_Ribo),SampleNames=Samples)
#Make list for paths of RNA-seq datasets 
RNAseqData=c(Root_RNA,Shoot_RNA)
#RNA-seq is paired-end or single-end?
RNAseqBamPairorSingle=c("paired","paired")

# check single-end data
# RNAseqData=c(Root_RNAse,Shoot_RNAse)
# RNAseqBamPairorSingle=c("single","single")

#Load example transcriptome annotation file
gtf_import(annotation=agtf,format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```
#### Plot different isoforms 
```
ggRibo(gene_id="AT4G21910",tx_id="AT4G21910.1",
       Y_scale="each",Extend=c(400,50),
       NAME = "MATE efflux family protein")
```
![image](https://github.com/user-attachments/assets/3aa258cb-718e-4a99-96da-359998f43c03)
```
ggRibo(gene_id="AT4G21910",tx_id="AT4G21910.2",
       Y_scale="each",Extend=c(400,50),
       NAME = "MATE efflux family protein")
```
![image](https://github.com/user-attachments/assets/c217a5ef-d2ff-4069-bdf7-a54c29ab7f22)

#### Plot a uORF
```
#Load CPuORF gtf
eORF_import(annotation=ugtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
ggRibo(gene_id="AT3G02470",tx_id="AT3G02470.1",
       eORF.tx_id = "AT3G02468.1",
       Y_scale="each",Extend=50,
       NAME = "SAMDC, CPuORF")
```
![image](https://github.com/user-attachments/assets/7dba3356-8459-493e-afcb-6dfd5cc69c0a)




```
