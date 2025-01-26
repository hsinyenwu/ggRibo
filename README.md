## *ggRibo*: a ggplot-based single-gene viewer for visualizing Ribo-seq and related omics datasets
### Introduction
Ribo-seq (ribosome profiling) is a powerful technique for studying mRNA translation by deep sequencing ribosome-protected footprints. A key feature of Ribo-seq data is 3-nucleotide periodicity, which reflects the ribosome’s codon-by-codon progression during translation. This 3-nucleotide periodicity facilitates the discovery of unannotated translation events and provides insights into translational regulation. Here, we present ggRibo, an R package designed for visualizing 3-nucleotide periodicity within a genomic context. ggRibo enables visual confirmation of translated and unannotated isoforms, as well as additional translation events, including upstream open reading frames (ORFs), downstream ORFs, stop codon readthrough, and correction of misannotated ORFs due to genome sequencing errors.  

Additionally, ggRibo allows for the comparison of Ribo-seq data with other sequencing methods that provide single nucleotide resolution (SNR), such as Translation Initiation sequencing (TI-seq), degradome sequencing (PARE-seq, Parallel Analysis of RNA Ends), GMUCT (Genome-Wide Mapping of Uncapped Transcripts), and Cap Analysis of Gene Expression sequencing (CAGE-seq). Some epitranscriptomic sequencing methods that detect the exact position of mRNA modifications, such as m6A-SAC-seq (N⁶-methyladenosine-Selective Alkylation Cleavage sequencing, which detects m6A sites) and BID-seq (bisulfite-induced deletion sequencing, which detects pseudouridine (Ψ) sites), are also SNR data. In SNR data, only one nucleotide position within each sequencing read carries the entire biological meaning for that read. For example, the P-site nucleotide of TI-seq reads indicates translation initiation sites on mRNAs, the first nucleotide of degradome-seq reads marks the 5’ end of RNA degradation intermediates, and the first nucleotide of CAGE-seq reads denotes transcription start sites.  

By integrating these diverse datasets, ggRibo enables researchers to identify factors that influence translation or are associated with translational processes, thereby facilitating the generation of hypotheses about the mechanisms governing diverse steps of gene expression and mRNA translation.  

### Plotting Ribo-seq reads
Each Ribo-seq read is represented with its first nucleotide aligned to the P-site (Figure 1A). The offset indicates the distance from the first nucleotide of the Ribo-seq read to the P-site of the ribosome. The offset can be obtained from metagene analysis of Ribo-seq reads using RiboTaper, Ribo-seQC, or other Ribo-seq analysis software. The cumulative P-site counts from all reads within the selected gene range were plotted (e.g., Figure 1B). Note the P-site offsets could vary in different organisms and organelles (see panels D-E).   
  
<img width="675" alt="image" src="https://github.com/user-attachments/assets/b1b16e9a-2a0d-45bd-b55e-77a4c5c68aad" />

### Gene-context plot vs single transcript plot for presenting Ribo-seq plots
Here we show one example gene with 3 isoforms (Figure 2A). Using the single transcript style plot, it is impossible to check which transcript(s) is translated (Figure 2B). The isoform 3 is not transcribed in the sample and leads to a confusing plot (bottom panel of Figure 2B). In gene-context plot, we can clear see the first and second isoforms are transcribed and translated (Figure 2C) even though only isoform 1 is colored for periodicity. Therefore, ***gene-context Ribo-seq plot provides a bird’s-eye view of the translation for all isoforms.*** 
  
<img width="675" alt="image" src="https://github.com/user-attachments/assets/7cbcacb4-a42d-45ab-bbcb-cd45bc1923a6" />

### How to read ggRibo plots  
The gene-context plot shown in Figure 2C is a ggRibo plot, where Ribo-seq reads are color-coded to demonstrate the 3-nucleotide periodicity: red for the first (expected/annotated) reading frame, blue for the second, and green for the third (Figures 2C). Reads outside the ORF range are displayed in gray. RNA-seq coverage is represented with a light yellow background (Figures 2C).   

### Steps and examples for the basic usage of ggRibo  

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

#### Load RNA-seq, Ribo-seq and annotation files  
1. For preparing your Ribo-seq reads for ggRibo, see [Here](https://github.com/hsinyenwu/ggRibo/blob/v2025.1.25/README_a0_preparing_Ribo-seq_for_ggRibo.md).  
2. RNA-seq files are just the bam files after Ribo-seq reads aligned with STAR or HISAT2.  
3. The FASTA and gtf/gff files are just the files you used to map the RNA-seq and Ribo-seq reads.  
4. Other Single Nucleotide Resolution data such as PARE-seq or TSS-seq (both 1st nt position and counts for each reads) could also be loaded for ggRibo plotting.  

**Here are the files for ggRibo plotting (Figure 3):**

<img width="675" alt="image" src="https://github.com/user-attachments/assets/efc31d7d-7d9f-4b01-a95e-f541183ddde6" />

**Load example files in the ggRibo package.**
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
(1) No column names.  
(2) Four columns for "counts", "chromosome number", "chromosome position", "strand" from left to right.  
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
It is clear that the root an shoot in Arabidopsis express different transcripts. And the 4th isoforms is likely not transcribed and translated.
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
# eORF means extra ORF. the eORF_import could be used to import gtf/gff3 for uORF, overlapping uORF, nested ORF, overlapping dORF and dORF.  

eORF_import(annotation=ugtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
ggRibo(gene_id="AT3G02470",tx_id="AT3G02470.1",
       eORF.tx_id = "AT3G02468.1",
       Y_scale="each",Extend=50,
       NAME = "SAMDC, CPuORF")
```
![image](https://github.com/user-attachments/assets/7dba3356-8459-493e-afcb-6dfd5cc69c0a)

#### Check sequences for the uORF
Download:  
(1) GTF (Araport11+CTRL_20181206.gtf)  
(2) FASTA (TAIR10_chr_all_2.fas)  
(3) RNA bam file (RNA_CTRL_merged.bam)   
(4) Ribo file (CTRL_expressed_P_sites_sort_count) from [here]  (https://data.mendeley.com/datasets/89j7snbm2r/2)  
Load data and import gtfs:
```
library(ggRibo)
CTRL_RNA="/path/to/RNA_CTRL_merged.bam"
CTRL_ribo="/path/to/CTRL_expressed_P_sites_sort_count"
FA <- FaFile("/path/to/TAIR10_chr_all_2.fas")
Samples = c("Seedlings")
RiboseqData = Ribo_data(c(CTRL_ribo),SampleNames=Samples)
RNAseqData = CTRL_RNA
RNAseqBamPairorSingle="paired"
gtf_import(annotation="/path/to/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
CiPS_TuORFs_gff3 <- system.file("extdata", "CiPS_TuORFs_Sep5d_2024.gff3", package = "ggRibo", mustWork = TRUE) #Load uORFs
eORF_import(annotation=CiPS_TuORFs_gff3, format="gff3",dataSource="Araport",organism="Arabidopsis thaliana")
```
Make the simple plot.
```
# plot the entire gene
ggRibo(
  gene_id = "AT3G50500",
  tx_id = "AT3G50500.1",
  NAME="SnRK2.2",
  Extend=50)
```
![image](https://github.com/user-attachments/assets/3b69990f-8e8e-4ef2-9324-689fde936fc4)
There is a strong peak in the 5'UTR suggesting the presence of a translated uORF.  
Show DNA sequence and focus on the uORF.  
**Need plot_range, show_seq = TRUE, FASTA**
```
# show minimum uORF
ggRibo(
  gene_id = "AT3G50500",
  tx_id = "AT3G50500.1",
  eORF.tx_id = "AT3G50500.1_227_232",
  NAME="SnRK2.2",
  plot_range = c(18743960,18743920),
  show_seq = TRUE,FASTA = FA,
  Extend=50)
```
![image](https://github.com/user-attachments/assets/3363a6b5-1447-470a-b9ce-d484782ca9ff)


#### Key parameters for ggRibo
(1) Extend (integer or a two integer vector): extend the plot range for both side of the plot. You can either use one number, which means same extension for both side, or use a vector with two values to extend left and right sides differently.  
(2) Y_scale (Boolean): the y-axis scale for each sample for the gene of interest. It could be "each", means each sample scale by itself to its max. The alternative is "all", means all samples are scaled together (same max Y-axis scale).  
(3) fExtend (integer): entend the 5' side of annotated CDS and also extend the frame of the annotated CDS. This is designed for visualizing non-AUG start.  
(4) tExtend (integer): entend the 3' side of annotated CDS and also extend the frame of the annotated CDS.  This is designed for visualizing stop codon readthrough.   
(5) eORF.tx_id (text): input the transcript id for extra ORFs. Remember the eORF gtf should be input with the eORF_import function and the transcript id for extra ORFs is included in the eORF gtf.   
(6) plot_genomic_direction (Boolean): plot the direction of the gene on the genome browser on top right side of the top plot.  
(7) sample_color (text vector): the color of the reads in each sample (from top to bottom). If you want the reads in the plot are color according to the 3 frames, use "color". Otherwise just give a single color. For example, if we provide: sample_color=c("color","purple"), reads in the first plot will be colored according to their frames, but all reads in the second plot will be colored purple. The default for all plots are "color".   
(8) frame_colors (text vector): colors for the 3 frames, default is c("0"="#FF0000", "1"="#3366FF", "2"="#009900"), you can choose the color you like.  

