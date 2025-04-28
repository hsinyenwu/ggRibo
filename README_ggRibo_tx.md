

```
#Install the new version and load test data (Arabidopsis) if you have not done so.
library(devtools)
install_github("hsinyenwu/ggRibo@v2025.3.30", dependencies = TRUE, force = TRUE)
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
#### Just ggRibo (isoform 3 is expressed)
```
ggRibo(gene_id="AT3G02470",tx_id="AT3G02470.3",
       eORF.tx_id = "AT3G02468.1",
       plot_ORF_ranges=T,
       NAME = "SAMDC, CPuORF”)
```

#### ggRibo_tx: Plot isoform 1 (not expressed, see an intron)
```
ggRibo_tx(gene_id="AT3G02470",tx_id="AT3G02470.1",
       eORF.tx_id = "AT3G02468.1",
       plot_ORF_ranges=T,
       NAME = "SAMDC, CPuORF")
```

#### Plot isoform 3 (expressed, no intron)
```
ggRibo_tx(gene_id="AT3G02470",tx_id="AT3G02470.3",
       eORF.tx_id = "AT3G02468.1",
       plot_ORF_ranges=T,
       NAME = "SAMDC, CPuORF")
```
#### Zoom in view
```
ggRibo_tx(gene_id="AT3G02470",tx_id="AT3G02470.3",
          eORF.tx_id = "AT3G02468.1",
          plot_ORF_ranges=T,
          plot_range=c(210,380),
          NAME = "SAMDC, CPuORF”)
```
