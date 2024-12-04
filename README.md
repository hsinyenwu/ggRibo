Install ggRibo and its required packages:

Install required packages.
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
Install RiboPlotR.
```
#Install RiboPlotR
library(devtools)
install_github("hsinyenwu/ggRibo")
```

```
# Load example datasets
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE) #Annotation
ugtf <- system.file("extdata", "AT3G02468.gtf", package = "ggRibo", mustWork = TRUE) #uORF annotation
RRNA <- system.file("extdata", "Root_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
SRNA <- system.file("extdata", "Shoot_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
RRNAse <- system.file("extdata", "Root_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
SRNAse <- system.file("extdata", "Shoot_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
RRibo <- system.file("extdata", "riboRoot.bed", package = "ggRibo", mustWork = TRUE) #Root Ribo-seq data
SRibo <- system.file("extdata", "riboShoot.bed", package = "ggRibo", mustWork = TRUE) #Shoot Ribo-seq data

FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
Samples=c("Root","Shoot")
RiboseqData=Ribo_data(c(RRibo,SRibo),SampleNames=Samples)
RNAseqData=c(RRNA,SRNA)
RNAseqBamPairorSingle=c("paired","paired")

# check single-end data
# RNAseqData=c(RRNAse,SRNAse)
# RNAseqBamPairorSingle=c("single","single")

gtf_import(annotation=agtf,format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
eORF_import(annotation=ugtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```
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


```
ggRibo(gene_id="AT3G02470",tx_id="AT3G02470.1",
       eORF.tx_id = "AT3G02468.1",
       Y_scale="each",Extend=50,
       NAME = "SAMDC, CPuORF")
```
![image](https://github.com/user-attachments/assets/7dba3356-8459-493e-afcb-6dfd5cc69c0a)




```
