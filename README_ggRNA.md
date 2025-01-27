### Plot RNA-seq reads with ggRNA
Here we use the example files included in the ggRibo package:
```
#Path for example data from ggRibo package
#path to annotated gtf
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE) #Annotated gtf
#path to paired-end RNA-seq reads
Root_RNA <- system.file("extdata", "Root_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNA <- system.file("extdata", "Shoot_test_PE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data
#path to single-end RNA-seq reads
Root_RNAse <- system.file("extdata", "Root_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Root RNA-seq data
Shoot_RNAse <- system.file("extdata", "Shoot_test_SE.bam", package = "ggRibo", mustWork = TRUE) #Shoot RNA-seq data

#Define sample names
Samples=c("Root_PE","Root_SE","Shoot_PE","Shoot_SE")
#Make list for paths of RNA-seq datasets
RNAseqData=c(Root_RNA,Root_RNAse,Shoot_RNA,Shoot_RNAse)
#RNA-seq is paired-end or single-end?
RNAseqBamPairorSingle=c("paired","single","paired","single")

#Load example transcriptome annotation file
gtf_import(annotation=agtf,format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```

**Plot with the y-axis from each sample.**
```
#Plot AT4G21910
ggRNA(gene_id="AT4G21910",tx_id="AT4G21910.1",Y_scale = "each",NAME="MATE efflux family protein",Extend=c(700,0))
```
![image](https://github.com/user-attachments/assets/5f933434-9973-4e54-9d6e-1406cff5f3c3)

**Plot with the max y-axis from all samples.**
```
#Plot AT4G21910
ggRNA(gene_id="AT4G21910",tx_id="AT4G21910.1",Y_scale = "all",NAME="MATE efflux family protein",Extend=c(700,0))
```
![image](https://github.com/user-attachments/assets/200450a8-6aab-4f2c-8d1d-00c159c02b25)

**Plot with different backfround**
```
#Plot AT4G21910
ggRNA(gene_id="AT4G21910",tx_id="AT4G21910.1",Y_scale = "all",NAME="MATE efflux family protein",Extend=c(700,0),RNAbackground = c("lightblue","lightgreen","yellow","red"))
```
![image](https://github.com/user-attachments/assets/cfc45a0b-8f6b-4b26-8651-4a011e1ba27b)


