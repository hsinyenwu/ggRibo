### Multiomics Visualization
Here we show the alternative transcription start sites for the Arabidopsis BCA4 gene:  
The RNA-seq examples are from Arabidopsis Root, Shoot and Seedling. The single-nucleotide resolution data are Arabidopsis Root Ribo-seq, Arabidopsis Shoot Ribo-seq, and Arabidopsis seedlings CAGE-seq.  
For Ribo-seq, the nucleotide position selected for each read count is the 1st nucleotide of the inferred P-site. See the [Introduction](https://github.com/hsinyenwu/ggRibo/blob/v2025.5.30/README.md) for detail.    
For CAGE-seq, the nucleotide position selected for each read count is the 1st position of the CAGE-seq read (i.e., transcription start site).   

***For visualizing single-nucleotide resolution (SNR) multiomics data:***  
***1. Identify/process the SNR data as the tabular/bigWig/bedGraph file as the Ribo-seq data.***   
***2. Include the SNR data as one of the Ribo-seq data for the "create_seq_input" function (see below)***  
```
#path to annotated gtf
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE)

#path to RNA-seq datasets
Root_RNA <- system.file("extdata", "Root_test_PE.bam", package = "ggRibo", mustWork = TRUE)
Shoot_RNA <- system.file("extdata", "Shoot_test_PE.bam", package = "ggRibo", mustWork = TRUE)
Seedling_RNA <- system.file("extdata", "RNA_CTRL_merged_sub34.bam", package = "ggRibo", mustWork = TRUE)

#path to Ribo-seq datasets
Root_Ribo <- system.file("extdata", "riboRoot.bed", package = "ggRibo", mustWork = TRUE) #Root Ribo-seq data
Shoot_Ribo <- system.file("extdata", "riboShoot.bed", package = "ggRibo", mustWork = TRUE) #Shoot Ribo-seq data
CAGE_seq <- system.file("extdata", "wt_R123_chr34.txt", package = "ggRibo", mustWork = TRUE) #Seedling CAGE-seq data
```
Load annotation with **gtf_import**.
```
#Load example transcriptome annotation file
gtf_import(annotation=agtf,format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```
Input data files with the **create_seq_input** function.  
```
#include CAGE_seq file in the "ribo_files" vector
inputs_full <- create_seq_input(
    rna_files = c(Root_RNA,Shoot_RNA,Seedling_RNA),
    ribo_files = c(Root_Ribo,Shoot_Ribo,CAGE_seq),
    sample_names = c("Root","Shoot","CAGE")
)
```

Example plot.
```
#data_types changes the right y-axis labels for your datasets. Default is "Ribo-seq"
ggRibo(tx_id="AT4G21910.1",
       Y_scale="each",Extend=c(400,50),
       data_types=c("Ribo-seq","Ribo-seq","CAGE-seq"),
       NAME = "MATE efflux family protein")
```
![image](https://github.com/user-attachments/assets/6656964a-574c-4eed-a1fc-5f0065ada8d5)  

Show CAGE-seq reads in blue.
```
#sample_color could change the nucleotide resolution read colors
#sample_color default is "color", which shows 3-nucleotide periodicity 
#Here use sample_color="blue" for CAGE-seq
ggRibo(tx_id="AT4G21910.1",
       Y_scale="each",Extend=c(400,50),
       data_types=c("Ribo-seq","Ribo-seq","CAGE-seq"),
       sample_color=c("color","color","blue"),
       NAME = "MATE efflux family protein")
```
![image](https://github.com/user-attachments/assets/6583a0bf-45ef-43b2-b556-74cd2f5a7795)
Since seedling sample is the mix of root and shoot samples, we see both isoforms expressed and translated in this sample.

### Visualizing other omics data independent of ggRibo
```
#path to annotated gtf
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE)
#path to RNA-seq datasets
Seedling_RNA <- system.file("extdata", "RNA_CTRL_merged_sub34.bam", package = "ggRibo", mustWork = TRUE)
#path to Ribo-seq datasets
CAGE_seq <- system.file("extdata", "wt_R123_chr34.txt", package = "ggRibo", mustWork = TRUE) #Seedling CAGE-seq data
```
Load annotation with **gtf_import**.  
```
#Load example transcriptome annotation file
gtf_import(annotation=agtf,format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```
Input CAGE data files with the **create_seq_input** function as above.  
```
#include CAGE_seq file in the "ribo_files" vector
inputs_full <- create_seq_input(
    rna_files = c(Seedling_RNA),
    ribo_files = c(CAGE_seq),
    sample_names = "Seedling CAGE"
)
```
Load CAGE-seq data. Note that there is no Ribo-seq data here.  
```
#Here use sample_color="blue" for CAGE-seq
ggRibo(tx_id="AT4G21910.1",
       Y_scale="each",Extend=c(400,50),
       data_types=c("CAGE-seq"),
       sample_color=c("blue"),
       NAME = "MATE efflux family protein")
```
![image](https://github.com/user-attachments/assets/666a53f6-2b91-449f-ad7f-6e973bfa1202)

### Note for other omics data.
For other multiomics data types, it is crucial to identify the nucleotide position(s) for quantification and plotting. Here are some examples:  
1. TI-seq for translation initiation sites: 1st position of the p-site.  
2. TSS-seq for transcription start sites: as CAGE-seq, uses 1st position of the sequencing read.  
3. m6A-SAC-seq and ac4C-seq for mRNA modifications: nucleotide positions contain the modifications.  
