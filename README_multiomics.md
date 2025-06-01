
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

```
#data_types changes the right y-axis labels for your datasets. Default is "Ribo-seq"
ggRibo(gene_id="AT4G21910",tx_id="AT4G21910.1",
       Y_scale="each",Extend=c(400,50),
       data_types=c("Ribo-seq","Ribo-seq","CAGE-seq"),
       NAME = "MATE efflux family protein")
```
![image](https://github.com/user-attachments/assets/6656964a-574c-4eed-a1fc-5f0065ada8d5)
