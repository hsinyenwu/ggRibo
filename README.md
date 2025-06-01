## *ggRibo*: a ggplot-based single-gene viewer for visualizing Ribo-seq and related omics datasets
[ggRibo Basic](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README.md)  
[Prepare_Ribo-seq data](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_0A_preparing_Ribo-seq_for_ggRibo.md)  
[Prepare bedGraph and bigWig for RNA-seq and Ribo-seq](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_0B_bedGraph_bigwig.md)  
#[Visualizing human data](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_human_ggRibo_examples.md)  
[Visualizing noncoding RNA](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_1_noncodingRNA.md)  
[Visualizing overlapping ORFs](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_2_overlapping_ORFs.md)  
[Visualizing downstream ORFs](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_3_dORF_translation.md)  
[Examples for ggRNA](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_ggRNA.md)  
[Examples for ggRibo_tx](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_ggRibo_tx.md)  
**[Multiomics visualization with ggRibo](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README_multiomics.md)**

### Introduction
Ribo-seq (ribosome profiling) is a powerful technique for studying mRNA translation by deep sequencing ribosome-protected footprints. A key feature of Ribo-seq data is 3-nucleotide periodicity, which reflects the ribosome’s codon-by-codon progression during translation. This 3-nucleotide periodicity facilitates the discovery of unannotated translation events and provides insights into translational regulation. Here, we present ggRibo, an R package designed for visualizing 3-nucleotide periodicity within a genomic context. ggRibo enables visual confirmation of translated and unannotated isoforms, as well as additional translation events, including upstream open reading frames (ORFs), downstream ORFs, stop codon readthrough, and correction of misannotated ORFs due to genome sequencing errors.  

Additionally, ggRibo allows for the comparison of Ribo-seq data with other sequencing methods that provide single nucleotide resolution (SNR), such as Translation Initiation sequencing (TI-seq), degradome sequencing (PARE-seq, Parallel Analysis of RNA Ends), GMUCT (Genome-Wide Mapping of Uncapped Transcripts), and Cap Analysis of Gene Expression sequencing (CAGE-seq). Some epitranscriptomic sequencing methods that detect the exact position of mRNA modifications, such as m6A-SAC-seq (N⁶-methyladenosine-Selective Alkylation Cleavage sequencing, which detects m6A sites) and BID-seq (bisulfite-induced deletion sequencing, which detects pseudouridine (Ψ) sites), are also SNR data. In SNR data, only one nucleotide position within each sequencing read carries the entire biological meaning for that read. For example, the P-site nucleotide of TI-seq reads indicates translation initiation sites on mRNAs, the first nucleotide of degradome-seq reads marks the 5’ end of RNA degradation intermediates, and the first nucleotide of CAGE-seq reads denotes transcription start sites.  

By integrating these diverse datasets, ggRibo enables researchers to identify factors that influence translation or are associated with translational processes, thereby facilitating the generation of hypotheses about the mechanisms governing diverse steps of gene expression and mRNA translation.  

### Plotting Ribo-seq reads
Each Ribo-seq read is represented with its first nucleotide aligned to the P-site (Figure 1A). The offset indicates the distance from the first nucleotide of the Ribo-seq read to the P-site of the ribosome. The offset can be obtained from metagene analysis of Ribo-seq reads using RiboTaper, Ribo-seQC, or other Ribo-seq analysis software. The cumulative P-site counts from all reads within the selected gene range were plotted (e.g., Figure 1B). Note the P-site offsets could vary in different organisms and organelles (see panels D-E).   
  
<img width="675" alt="image" src="https://github.com/user-attachments/assets/b1b16e9a-2a0d-45bd-b55e-77a4c5c68aad" />

### Gene-context plot vs single transcript plot for presenting Ribo-seq plots
Here we show one example gene with 3 isoforms (Figure 2A). Using the single transcript style plot, it is impossible to check which transcript(s) is translated (Figure 2B). The isoform 3 is not transcribed in the sample and leads to a confusing plot (bottom panel of Figure 2B). In gene-context plot, we can clear see the first and second isoforms are transcribed and translated (Figure 2C) even though only isoform 1 is colored for periodicity. Therefore, ***gene-context Ribo-seq plot provides a bird’s-eye view of the translation for all isoforms.*** ***However, single transcript plot could still be helpful for genes have long introns. In that case, you can check the Gene-context plot first to check expressed isoforms then use single transcript plot for presentation.***
  
<img width="675" alt="image" src="https://github.com/user-attachments/assets/7cbcacb4-a42d-45ab-bbcb-cd45bc1923a6" />

### How to understand a ggRibo plot   
The gene-context plot shown in Figure 2C is a ggRibo plot, where Ribo-seq reads are color-coded to demonstrate the 3-nucleotide periodicity: red for the first (expected/annotated) reading frame, blue for the second, and green for the third (Figures 2C). Reads outside the ORF range are displayed in gray. RNA-seq coverage is represented with a light yellow background (Figures 2C).   

### Steps and examples for the basic usage of ggRibo  

#### Install ggRibo and its required packages:  

(1) Install required packages.
```
## 1. Make sure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

## 2. Bioconductor packages
bioc_pkgs <- c(
    "GenomicRanges",
    "GenomicFeatures",
    "GenomicAlignments",
    "Rsamtools",
    "IRanges",
    "txdbmaker",
    "rtracklayer",     
    "Biostrings",      
    "GenomeInfoDb",   
    "BSgenome"       
)

BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())))

## 3. CRAN packages
cran_pkgs <- c("ggplot2", "cowplot", "dplyr", "R6")

install.packages(setdiff(cran_pkgs, rownames(installed.packages())))
```
(2) Install ggRibo.
```
#Install ggRibo
library(devtools)
install_github("hsinyenwu/ggRibo")
```

#### Load RNA-seq, Ribo-seq and annotation files  
1. Ribo-seq input could be a tabular format with 4 columns for (1) read counts, (2) chromosome, (3) position of the 1st nucleotide of P-site and (4) strand. Alternatively, you can also input bedGraph or bigWig format files. For preparing files for ggRibo, see [Here](https://github.com/hsinyenwu/ggRibo/blob/v2025.1.25/README_a0_preparing_Ribo-seq_for_ggRibo.md).  
2. RNA-seq files could be the bam files from RNA-seq reads aligned with STAR or HISAT2. You can also convert your data to bedGraph or bigWig formats, see [Here](https://github.com/hsinyenwu/ggRibo/blob/v2025.3.30/README_0B_bedGraph_bigwig.md). 
3. The FASTA (or a BSGenome object) and gtf/gff files for visualizing DNA and amino acid sequences.  
4. Other Single Nucleotide Resolution data such as PARE-seq or TSS-seq could also be loaded for ggRibo plotting. Similar to Ribo-seq, the SNR data could be the tabular format, bedGraph or bigWig.  

**Here are the files for ggRibo plotting (Figure 3):**

<img width="675" alt="image" src="https://github.com/user-attachments/assets/efc31d7d-7d9f-4b01-a95e-f541183ddde6" />

**Load package and example files in the ggRibo package.**
```
library(ggRibo)
#Path for example data from ggRibo package
#path to annotated gtf
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE) #Annotated gtf
#path to the gtf files for uORFs
ugtf <- system.file("extdata", "AT3G02468.gtf", package = "ggRibo", mustWork = TRUE) #uORF gtf
#path to paired-end RNA-seq reads
Root_RNA <- system.file("extdata", "Root_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNA <- system.file("extdata", "Shoot_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
#path to single end RNA-seq reads
Root_RNAse <- system.file("extdata", "Root_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNAse <- system.file("extdata", "Shoot_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
#path to Ribo-seq reads
Root_Ribo <- system.file("extdata", "riboRoot.bed", package = "ggRibo", mustWork = TRUE) #Root Ribo-seq data
Shoot_Ribo <- system.file("extdata", "riboShoot.bed", package = "ggRibo", mustWork = TRUE) #Shoot Ribo-seq data
#You need to create paths for your own data files
```

**Load transcriptome annotation:**  
```
#Load example transcriptome annotation file
gtf_import(annotation=agtf,format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```

**Setup variables for the ggRibo function:**  
```
Samples=c("Root","Shoot")
# input single-end data
inputs_full <- create_seq_input(
  rna_files = c(Root_RNAse,Shoot_RNAse),
  ribo_files = c(Root_Ribo,Shoot_Ribo),
  sample_names = Samples,
  rna_paired = c("single","single") # The default for "rna_paired" is "paired" for all samples.
)

# Plot with ggRibo (Example)
ggRibo(
  gene_id = "AT4G21910",
  tx_id = "AT4G21910.1"
)
#Input paired-end data
inputs_full <- create_seq_input(
  rna_files = c(Root_RNA,Shoot_RNA),
  ribo_files = c(Root_Ribo,Shoot_Ribo),
  sample_names = Samples
) # The default for "rna_paired" is "paired" for all samples so we do not need to change it here.

# Plot with ggRibo (Example)
ggRibo(
  gene_id = "AT4G21910",
  tx_id = "AT4G21910.1"
)
```

#### Plot different isoforms 
The result below shows that the root an shoot in Arabidopsis express different transcripts. And the 4th isoforms is likely not transcribed and translated.
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
![image](https://github.com/user-attachments/assets/52b62a05-2d93-4f34-8b63-b2485a6555b5)

#### Check sequences for the uORF
Download annotation and data files from [here](https://data.mendeley.com/datasets/89j7snbm2r/2):  
(1) GTF (Araport11+CTRL_20181206.gtf)  
(2) FASTA (TAIR10_chr_all_2.fas) #you can also use a BSGenome object 
(3) RNA bam file (RNA_CTRL_merged.bam)   
(4) Ribo file (CTRL_expressed_P_sites_sort_count)




Load data and import gtfs:  
```
library(ggRibo)
CTRL_RNA="/path/to/RNA_CTRL_merged.bam"
CTRL_ribo="/path/to/CTRL_expressed_P_sites_sort_count"
FA <- FaFile("/path/to/TAIR10_chr_all_2.fas")
RNA_files <- list(CTRL_RNA)
Ribo_files <- list(CTRL_ribo)

# Prepare coverage descriptors
inputs_full <- create_seq_input(
  rna_files = RNA_files,
  ribo_files = Ribo_files,
  sample_names="CTRL")

gtf_import(annotation="/path/to/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
CiPS_TuORFs_gff3 <- system.file("extdata", "CiPS_TuORFs_Sep5d_2024.gff3", package = "ggRibo", mustWork = TRUE) #Load uORFs
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
#Input minimum uORF gtf
eORF_import(annotation=CiPS_TuORFs_gff3, format="gff3",dataSource="Araport",organism="Arabidopsis thaliana")
# show minimum uORF
ggRibo(
  gene_id = "AT3G50500",
  tx_id = "AT3G50500.1",
  eORF.tx_id = "AT3G50500.1_227_232",
  NAME="SnRK2.2",
  plot_range = c(18743960,18743920), #select plotting range
  show_seq = TRUE,FASTA = FA, #'show_seq = TRUE' means you want to see the sequence, than you need to input FASTA using 'FASTA=FA'
  Extend=50)
#You could also use a BSGenome object for FASTA input, For example, FASTA=BSGenome.Athaliana 
```
![image](https://github.com/user-attachments/assets/0bc4ee83-633c-4082-97af-e97cfa5f4a5c)

**Use nucleotide_color_scheme="colorblind" for an alternative coloring scheme**
```
ggRibo(
    gene_id = "AT3G50500",
    tx_id = "AT3G50500.1",
    eORF.tx_id = "AT3G50500.1_227_232",
    NAME="SnRK2.2",
    plot_range = c(18743960,18743920), #select plotting range
    show_seq = TRUE,FASTA = FA, #'show_seq = TRUE' means you want to see the sequence, than you need to input FASTA using 'FASTA=FA'
    Extend=50, nucleotide_color_scheme="colorblind")
```
![image](https://github.com/user-attachments/assets/1fe171f2-4036-4a9e-8fee-2f8a466ce8e5)

#### Key parameters for ggRibo
(1) Extend (integer or a two integer vector): extend the plot range for both side of the plot. You can either use one number, which means same extension for both side, or use a vector with two values to extend left and right sides differently.  
(2) Y_scale (Boolean): the y-axis scale for each sample for the gene of interest. It could be "each", means each sample scale by itself to its max. The alternative is "all", means all samples are scaled together (same max Y-axis scale).  
(3) fExtend (integer): entend the 5' side of annotated CDS and also extend the frame of the annotated CDS. This is designed for visualizing non-AUG start.  
(4) tExtend (integer): entend the 3' side of annotated CDS and also extend the frame of the annotated CDS.  This is designed for visualizing stop codon readthrough.   
(5) eORF.tx_id (text): input the transcript id for extra ORFs. Remember the eORF gtf should be input with the eORF_import function and the transcript id for extra ORFs is included in the eORF gtf.   
(6) plot_genomic_direction (Boolean): plot the direction of the gene on the genome browser on top right side of the top plot.  
(7) sample_color (text vector): the color of the reads in each sample (from top to bottom). If you want the reads in the plot are color according to the 3 frames, use "color". Otherwise just give a single color. For example, if we provide: sample_color=c("color","purple"), reads in the first plot will be colored according to their frames, but all reads in the second plot will be colored purple. The default for all plots are "color".   
(8) frame_colors (text vector): colors for the 3 frames, default is c("0"="#FF0000", "1"="#3366FF", "2"="#009900"), you can choose the color you like.  
(9) selected_isoforms (text vector): you can select certain isoforms to plot. 
(10) data_types (text vector): This parameter is for the right Y-axis labels. Default for all data is Ribo-seq.
(11) dna_aa_height_ratio (numeric): change is height of DNA/AA plot.
(12) gene_model_height_ratio (numeric): change is height of transcript model plot.
(13) show_seq (Boolean): show DNA (when plot range <=201 nucleotides) and AA sequences 
(14) FASTA (a genome FASTA file or a BSGenome object): contain genomic sequences, required when show_seq=T
(15) plot_range (2 inegter numeric vector): defined the range for the plot.
(16) oORF_coloring:	Character string specifying coloring method for overlapping ORFs ("oORF_colors" or "extend_mORF"). "oORF_colors" means only show the 3-nt periodicity for the oORF. "extend_mORF" means extend the mORF frames to cover the extra ORFs.

## Citation: [ggRibo: a ggplot-based single-gene viewer for visualizing Ribo-seq and related omics datasets](https://www.biorxiv.org/content/10.1101/2025.01.30.635743v1)


