# Plotting human Ribo-seq data with ggRibo
Ribo-seq and RNA-seq data for embryonic stem cells (ESC) and brain tissue from [Chothani et al., 2022](https://doi.org/10.1016/j.molcel.2022.06.023) acquired from NCBI BioProjects [PRJNA756018](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA756018) and [PRJNA756023](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA756023).

Ribo-seq and RNA-seq data from induced pluripotent stem cells (iPSC) and cardiomyocytes from [Chen et al., 2020](https://www.science.org/doi/10.1126/science.aay0262) acquired from NCBI BioProject [PRJNA544411](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA544411)

Aligned to the GRCh38 human genome acquired from [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)

## Example of uORF in humans
The uORF in MRPL11 has been shown to inhibit translation of its downstream ORF ([Calvo, Pagliarini, & Mootha, 2009](https://www.pnas.org/doi/10.1073/pnas.0810916106)).

```
riboseq <- "ESC/mergedBamFile.bam_P_sites.ggRibo"
rnaseq <- "ESC/small_merged_rna_subset.bam"
rna_type <- "paired"

Samples <- "Human ESC"
RiboseqData <- Ribo_data(riboseq,SampleNames=Samples)
RNAseqData <- rnaseq
RNAseqBamPairorSingle <- rna_type

gtf_path <- "Homo_sapiens.GRCh38.112.chr.gff3"
gtf_import(annotation = gtf_path, format = "gff3", organism = "Homo sapiens")

uorf_path <- "ENST00000310999_uORF.gtf"
eORF_import(annotation = uorf_path, format = "gtf", organism = "Homo sapiens")

ggRibo(gene_id = "ENSG00000174547", tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999",
       NAME = "MRPL11")
```
![uORF1_MRPL11_Human_ESC_full_gene](https://github.com/user-attachments/assets/b1c7c941-ba05-40b7-b734-88721ddd21a1)

## Add genome/peptide sequence
```
fasta <- FaFile("Homo_sapiens.GRCh38.dna.primary_assembly.fa")

ggRibo(gene_id = "ENSG00000174547", tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999",
       NAME = "MRPL11", FASTA = fasta, show_seq = T)
```
![uORF2_MRPL11_Human_ESC_full_with_seq](https://github.com/user-attachments/assets/55688691-b0cf-4988-bc9a-14ba22781771)

## Closer look at uORF (adjusting plot range)
```
range <- c(66438710, 66438850)
ggRibo(gene_id = "ENSG00000174547", tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999",
       NAME = "MRPL11", FASTA = fasta, show_seq = T,
       plot_range = range, Extend = 0)
```
![uORF3_MRPL11_Human_ESC_zoom_with_seq](https://github.com/user-attachments/assets/5f10f2ae-864e-4cfc-8b56-28d6caf8aaae)

## Multiple samples
```
rna_type <- c("single", "single")
riboseq <- c("IPSC/mergedBamFile.bam_P_sites.ggRibo",
             "Cardio/mergedBamFile.bam_P_sites.ggRibo")
rnaseq <- c("IPSC/small_merged_rna_subset.bam",
            "Cardio/small_merged_rna_subset.bam")

Samples <- c("iPSCs", "Cardiomyocytes")
RiboseqData <- Ribo_data(riboseq,SampleNames=Samples)
RNAseqData <- rnaseq
RNAseqBamPairorSingle <- rna_type

gtf_path <- "Homo_sapiens.GRCh38.112.chr.gff3"
gtf_import(annotation = gtf_path, format = "gff3", organism = "Homo sapiens")

ggRibo(gene_id = "ENSG00000106211", tx_id = "ENST00000248553", NAME = "HSPB1")
ggRibo(gene_id = "ENSG00000133112", tx_id = "ENST00000530705", NAME = "TPT1")
```
<ins>Higher translation in Cardiomyocytes:</ins>
![TwoSample1_HSPB1_Human](https://github.com/user-attachments/assets/846a182b-3325-43a9-8390-109a0f3dceee)

<ins>Higher translation in iPSCs:</ins>
![TwoSample2_TPT1_Human](https://github.com/user-attachments/assets/48489cdb-f1b8-4a11-b57f-a3be5ce74d6e)

## Plotting subset of transcripts
Many genes are annotated with a large number of isoforms.
This can be a nuisance when plotting transcripts as the unexpressed isoforms can take up space,
and expand the x-axis beyond the ideal range if there are distant exons.

While it is good practice to view all isoforms to see which are the most important, it may be desirable to ultimately subset the number of transcripts plotted:

```
riboseq <- "Brain/mergedBamFile.bam_P_sites.ggRibo"
rnaseq <- "Brain/small_merged_rna_subset.bam"
rna_type <- "paired"

Samples <- "Human brain"
RiboseqData <- Ribo_data(riboseq,SampleNames=Samples)
RNAseqData <- rnaseq
RNAseqBamPairorSingle <- rna_type

gtf_path <- "Homo_sapiens.GRCh38.112.chr.gff3"
gtf_import(annotation = gtf_path, format = "gff3", organism = "Homo sapiens")

### Plot all 27 isoforms of IFITM3
ggRibo(gene_id = "ENSG00000142089", tx_id = "ENST00000399808", NAME = "IFITM3")

### Plot subset of isoforms
tx_subset <- c("ENST00000399808", "ENST00000680209", "ENST00000681198",
               "ENST00000681304", "ENST00000681840")  # list of desired isoforms

ggRibo(gene_id = "ENSG00000142089", tx_id = "ENST00000399808", NAME = "IFITM3",
       selected_isoforms = tx_subset)  # use subset of isoforms

```
<ins>All 27 isoforms:</ins>
![IsoSubset_IFITM3_human_allIsoforms](https://github.com/user-attachments/assets/eb5d290a-f92d-4bb8-a6de-fb226b4699d6)
<ins>Subset of 5 isoforms:</ins>
![IsoSubset_IFITM3_human_subset](https://github.com/user-attachments/assets/7dbcaf44-d36f-4276-be1f-18b5613b16da)
