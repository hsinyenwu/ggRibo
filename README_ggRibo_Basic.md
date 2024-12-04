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
(2) Install ggRibo.
```
#Install ggRibo
library(devtools)
install_github("hsinyenwu/ggRibo")
```
### Basic usage of ggRibo
#### Load RNA-seq, Ribo-seq and annotation files 
```
#Path for example data from ggRibo package
#You can create paths for your own data files
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE) #Annotation
ugtf <- system.file("extdata", "AT3G02468.gtf", package = "ggRibo", mustWork = TRUE) #uORF gtf
Root_RNA <- system.file("extdata", "Root_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNA <- system.file("extdata", "Shoot_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
Root_RNAse <- system.file("extdata", "Root_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNAse <- system.file("extdata", "Shoot_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
Root_Ribo <- system.file("extdata", "riboRoot.bed", package = "ggRibo", mustWork = TRUE) #Root Ribo-seq data
Shoot_Ribo <- system.file("extdata", "riboShoot.bed", package = "ggRibo", mustWork = TRUE) #Shoot Ribo-seq data
```
**Note the riboseq files is a table with the following  organization:**
(1) No column names  
(2) Four columns for "counts", "chromosome number", "chromosome position", "strand" from left to right 
*The riboseq file contain the number and distribution of the 1st position of the P-site of ribosome footprint. Alternatively, you could create files from other sequencing data with single-nucleotide resolution (SNR). SNR data include, but not limit to, PARE-seq, CAGE-seq, or TI-seq, which defines the 5' nucleotide of mRNA degredation intermediates, the 5' CAP positions or the transcription start sites, or transation initiation sites, respectively.*
```
1   1  1000000      +
3   1 10000007      +
3   1 10000010      +
3   1 10000016      +
1   1 10000018      +
4   1 10000019      +
```

**Setup variables for the ggRibo function:**  
```
#Define sample names  
Samples=c("Root","Shoot")
#Load Ribo-seq data
RiboseqData=Ribo_data(c(Root_Ribo,Shoot_Ribo),SampleNames=Samples)
#Make list for paths of RNA-seq datasets 
RNAseqData=c(Root_RNA,Shoot_RNA)
#RNA-seq is paired-end or single-end?
RNAseqBamPairorSingle=c("paired","paired")

# Please do not change the names of the following variables "Samples", "RiboseqData", "RNAseqData", "RNAseqBamPairorSingle" or you have to change the input variable names for the ggRibo function.

# check single-end data
# RNAseqData=c(Root_RNAse,Shoot_RNAse)
# RNAseqBamPairorSingle=c("single","single")
```
**Load transcriptome annotation:**  
```
#Load example transcriptome annotation file
gtf_import(annotation=agtf,format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```
**Please do not change the names of the following variables "Samples", "RiboseqData", "RNAseqData", "RNAseqBamPairorSingle" or you have to change the corresponding input parameters for the ggRibo function.**

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

#### Check sequences for the uORF
Download GTF (Araport11+CTRL_20181206.gtf) and FASTA (TAIR10_chr_all_2.fas) from [here](https://data.mendeley.com/datasets/89j7snbm2r/2)

#### Key parameters for ggRibo
(1) Extend (integer or a two integer vector): extend the plot range for both side of the plot. You can either use one number, which means same extension for both side, or use a vector with two values to extend left and right sides differently.  
(2) Y_scale (Boolean): the y-axis scale for each sample for the gene of interest. It could be "each", means each sample scale by itself to its max. The alternative is "all", means all samples are scaled together (same max Y-axis scale).  
(3) fExtend (integer): entend the 5' side of annotated CDS and also extend the frame of the annotated CDS. This is designed for visualizing non-AUG start.  
(4) tExtend (integer): entend the 3' side of annotated CDS and also extend the frame of the annotated CDS.  This is designed for visualizing stop codon readthrough. 

