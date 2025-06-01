
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

```
inputs_full <- create_seq_input(
    rna_files = c(Root_RNA,Shoot_RNA,Seedling_RNA),
    ribo_files = c(Root_Ribo,Shoot_Ribo,CAGE_seq),
    sample_names = Samples,
    rna_types = rep("bam", 3),
    ribo_types = rep("tabular", 3)
)
```

Example plot.
```
#data_types changes the right y-axis labels for your datasets. Default is "Ribo-seq"
ggRibo(gene_id="AT4G21910",tx_id="AT4G21910.1",
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
ggRibo(gene_id="AT4G21910",tx_id="AT4G21910.1",
       Y_scale="each",Extend=c(400,50),
       data_types=c("Ribo-seq","Ribo-seq","CAGE-seq"),
       sample_color=c("color","color","blue"),
       NAME = "MATE efflux family protein")
```
![image](https://github.com/user-attachments/assets/6583a0bf-45ef-43b2-b556-74cd2f5a7795)


