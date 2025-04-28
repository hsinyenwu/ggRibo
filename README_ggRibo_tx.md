ggRibo_tx provided a tool to visualize the ribo-seq reads (colored for periodicity) and RNA-seq coverage for a single transcript. While we prefer using single gene view with ggRibo, some genes have long introns and difficult to visualize with ggRibo. In that case, ggRibo_tx will be useful. However, visualizing with ggRibo first for checking expressed isoforms is still recommended.    
Here is an example:
```
#Install the new version and load test data (Arabidopsis) if you have not done so.
library(devtools)
install_github("hsinyenwu/ggRibo@v2025.3.30", dependencies = TRUE, force = TRUE)
library(ggRibo)
#path to annotated gtf
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE)
#path to RNA-seq bigwig coverage files
#Remember two files per sample (one for plus strand and one for minus strand)
Root_RNA_p <- system.file("extdata", "Root_plus_strand.bw", package = "ggRibo", mustWork = TRUE)
Root_RNA_m <- system.file("extdata", "Root_minus_strand.bw", package = "ggRibo", mustWork = TRUE)
Shoot_RNA_p <- system.file("extdata", "Shoot_plus_strand.bw", package = "ggRibo", mustWork = TRUE)
Shoot_RNA_m <- system.file("extdata", "Shoot_minus_strand.bw", package = "ggRibo", mustWork = TRUE)
#path to Ribo-seq bigwig coverage files
#Remember two files per sample (one for plus strand and one for minus strand)
Root_Ribo_p <- system.file("extdata", "riboRoot_plus.bw", package = "ggRibo", mustWork = TRUE)
Root_Ribo_m <- system.file("extdata", "riboRoot_minus.bw", package = "ggRibo", mustWork = TRUE)
Shoot_Ribo_p <- system.file("extdata", "riboShoot_plus.bw", package = "ggRibo", mustWork = TRUE)
Shoot_Ribo_m <- system.file("extdata", "riboShoot_minus.bw", package = "ggRibo", mustWork = TRUE)
# Load gtf
gtf_import(annotation=agtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

#path to the gtf files for uORFs
ugtf <- system.file("extdata", "AT3G02468.gtf", package = "ggRibo", mustWork = TRUE) #uORF gtf
eORF_import(annotation=ugtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

# Load necessary data (assuming Txome_Range and other objects are prepared)
RNA_files <- list(
  list(plus = Root_RNA_p, minus = Root_RNA_m),
  list(plus = Shoot_RNA_p,  minus = Shoot_RNA_m)
)

Ribo_files <- list(
  list(plus = Root_Ribo_p, minus = Root_Ribo_m),
  list(plus = Shoot_Ribo_p,  minus = Shoot_Ribo_m)
)

# Prepare coverage descriptors
inputs_full <- create_seq_input(
  rna_files = RNA_files,
  ribo_files = Ribo_files,
  sample_names = c("Root", "Shoot")
)
```
#### Just ggRibo (isoform 3 is expressed) for the entire gene.
```
ggRibo(gene_id="AT3G02470",tx_id="AT3G02470.3",
       eORF.tx_id = "AT3G02468.1",
       plot_ORF_ranges=T,
       NAME = "SAMDC, CPuORF")
```
<img src=https://github.com/user-attachments/assets/cd27d8ae-98e5-4659-90c5-ff19f9ee4b43 width="300" height="200">

#### ggRibo_tx: Plot isoform 1 (not expressed, see a strong intron in the plot)
Since this is a single transcript plot (e.g. plotting for cDNA), we should not see a strong intron.
```
ggRibo_tx(gene_id="AT3G02470",tx_id="AT3G02470.1",
       eORF.tx_id = "AT3G02468.1",
       plot_ORF_ranges=T,
       gene_model_height_ratio=1.9,
       NAME = "SAMDC, CPuORF")
```
![image](https://github.com/user-attachments/assets/a1d65c5f-bab0-4099-adf8-5f940932c77c)

#### Plot isoform 3 (expressed, no intron showen in the plot)
```
ggRibo_tx(gene_id="AT3G02470",tx_id="AT3G02470.3",
       eORF.tx_id = "AT3G02468.1",
       plot_ORF_ranges=T,
       gene_model_height_ratio=1.9,
       NAME = "SAMDC, CPuORF")
```
![image](https://github.com/user-attachments/assets/0d05810c-292f-497f-a14d-1819b2e7772a)

#### Plot isoform 4 (not expressed, miss RNA-seq coverage in the first exon)
```
ggRibo_tx(gene_id="AT3G02470",tx_id="AT3G02470.3",
       eORF.tx_id = "AT3G02468.1",
       plot_ORF_ranges=T,
       gene_model_height_ratio=1.9,
       NAME = "SAMDC, CPuORF")
```
![image](https://github.com/user-attachments/assets/3f3271af-2c2a-4744-83fa-5c7137f25e5d)

#### Zoom in view
```
ggRibo_tx(gene_id="AT3G02470",tx_id="AT3G02470.3",
          eORF.tx_id = "AT3G02468.1",
          plot_ORF_ranges=T,
          plot_range=c(210,380),
          gene_model_height_ratio=1.9,
          NAME = "SAMDC, CPuORF")
```
![image](https://github.com/user-attachments/assets/a6335df4-2c7f-4a48-bea0-43e32d8c525b)


