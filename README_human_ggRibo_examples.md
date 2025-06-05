# Plotting human Ribo-seq data with ggRibo
Ribo-seq and RNA-seq data for embryonic stem cells (ESC) and brain tissue from [Chothani et al., 2022](https://doi.org/10.1016/j.molcel.2022.06.023) acquired from NCBI BioProjects [PRJNA756018](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA756018) and [PRJNA756023](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA756023).

Ribo-seq and RNA-seq data from induced pluripotent stem cells (iPSC) and cardiomyocytes from [Chen et al., 2020](https://www.science.org/doi/10.1126/science.aay0262) acquired from NCBI BioProject [PRJNA544411](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA544411)

Aligned to the GRCh38 human genome acquired from [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index)

***You can find processed human data for testing here:[Human Data](https://data.mendeley.com/datasets/m3t293k4wr/1)***

## Example of uORF in humans
The uORF in MRPL11 has been shown to inhibit translation of its downstream ORF ([Calvo, Pagliarini, & Mootha, 2009](https://www.pnas.org/doi/10.1073/pnas.0810916106)).

```
gtf_path <- "Homo_sapiens.GRCh38.112.chr.gff3"
gtf_import(annotation = gtf_path, format = "gff3", organism = "Homo sapiens")  # Load annotation

uorf_path <- "ENST00000310999_uORF.gtf"
eORF_import(annotation = uorf_path, format = "gtf", organism = "Homo sapiens")  # Load uORF coordinates

inputs_full <- create_seq_input(
  rna_files = "ESC/small_merged_rna_subset.bam",
  ribo_files = "ESC/mergedBamFile.bam_P_sites.ggRibo",
  sample_names = "Human ESC",
  rna_paired = "paired"
)  # load sequencing data

# Use ggRibo to plot data in gene structure context
ggRibo(tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999", NAME = "MRPL11")
```
![uORF1_MRPL11_Human_ESC_full_gene](https://github.com/user-attachments/assets/0c68b41c-f145-4988-8697-cb6ddbfaa575)

## Add genome/peptide sequence
```
fasta <- FaFile("Homo_sapiens.GRCh38.dna.primary_assembly.fa")  # set path to genome FASTA

ggRibo(tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999", NAME = "MRPL11",
       FASTA = fasta, show_seq = T)
```
![uORF2_MRPL11_Human_ESC_full_with_seq](https://github.com/user-attachments/assets/72fc96fd-6e8e-489f-b6a9-3c0f4cbc6607)

## Closer look at uORF (adjusting plot range)
```
range <- c(66438710, 66438850)  # define a range on the chromosome
ggRibo(tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999", NAME = "MRPL11",
       FASTA = fasta, show_seq = T,
       plot_range = range, Extend = 0)
```
![uORF3_MRPL11_Human_ESC_zoom_with_seq](https://github.com/user-attachments/assets/5f151cda-291d-4d1d-a281-32b82b0bd378)

## Transcript view
Using transcript view, we can take a close look at specific isoforms.
```
ggRibo_tx(tx_id="ENST00000310999", eORF.tx_id = "ENST00000310999", NAME = "MRPL11",
          show_seq = T, FASTA = fasta)
ggRibo_tx(tx_id="ENST00000329819", eORF.tx_id = "ENST00000310999", NAME = "MRPL11",
          show_seq = T, FASTA = fasta)
```
<ins>ENST00000310999 (most abundant isoform):</ins>  
![uORF4_MRPL11_Human_ESC_tx1](https://github.com/user-attachments/assets/2bb86858-2923-4ff1-a5af-c1f17a726648)

<ins>ENST00000329819 (lowly translated isoform):</ins>  
![uORF5_MRPL11_Human_ESC_tx2_zoom](https://github.com/user-attachments/assets/870bef04-f62f-48af-a36d-06d9ec8a844f)

Transcript view is also useful for genes in which the introns are very long relative to the exons
```
ggRibo(tx_id = "ENST00000311672", NAME = "UQCRH")  # gene view
ggRibo_tx(tx_id = "ENST00000311672", NAME = "UQCRH")  # transcript view
```
<ins>Gene view with large introns:</ins>  
![Intron_UQCRH_Human_ESC_gene](https://github.com/user-attachments/assets/1f5fb7d5-e501-407a-be52-76d85f6034c7)
<ins>Transcript view, introns removed:</ins>  
![Intron_UQCRH_Human_ESC_tx](https://github.com/user-attachments/assets/6876114e-346a-4f1f-8c7b-b348c07ec2e6)

## Multiple samples
```
inputs_full <- create_seq_input(
  rna_files = c("IPSC/small_merged_rna_subset.bam",
                "Cardio/small_merged_rna_subset.bam"),
  ribo_files = c("IPSC/mergedBamFile.bam_P_sites.ggRibo",
                 "Cardio/mergedBamFile.bam_P_sites.ggRibo"),
  sample_names = c("iPSCs", "Cardiomyocytes"),
  rna_paired = c("single", "single")
)

gtf_path <- "Homo_sapiens.GRCh38.112.chr.gff3"
gtf_import(annotation = gtf_path, format = "gff3", organism = "Homo sapiens")

ggRibo(tx_id = "ENST00000248553", NAME = "HSPB1")
ggRibo(tx_id = "ENST00000530705", NAME = "TPT1")
```
<ins>Higher translation in Cardiomyocytes:</ins>  
![TwoSample1_HSPB1_Human](https://github.com/user-attachments/assets/29f93c10-f36d-4492-b638-d57b535ebf52)

<ins>Higher translation in iPSCs:</ins>  
![TwoSample2_TPT1_Human](https://github.com/user-attachments/assets/1ccfea33-ddd8-423d-a9b9-a1d5dfb530e6)

## Plotting subset of transcripts
Many genes are annotated with a large number of isoforms.
This can be a nuisance when plotting transcripts as the unexpressed isoforms can take up space,
and expand the x-axis beyond the ideal range if there are distant exons.

While it is good practice to view all isoforms to see which are the most important, it may be desirable to ultimately subset the number of transcripts plotted:

```
inputs_full <- create_seq_input(
  rna_files = "Brain/small_merged_rna_subset.bam",
  ribo_files = "Brain/mergedBamFile.bam_P_sites.ggRibo",
  sample_names = "Human ESC",
  rna_paired = "paired"
)

gtf_path <- "Homo_sapiens.GRCh38.112.chr.gff3"
gtf_import(annotation = gtf_path, format = "gff3", organism = "Homo sapiens")

### Plot all 27 isoforms of IFITM3
ggRibo(tx_id = "ENST00000399808", NAME = "IFITM3")

### Plot subset of isoforms
tx_subset <- c("ENST00000399808", "ENST00000680209", "ENST00000681198",
               "ENST00000681304", "ENST00000681840")  # list of desired isoforms

ggRibo(tx_id = "ENST00000399808", NAME = "IFITM3",
       selected_isoforms = tx_subset)  # use subset of isoforms

```
<ins>All 27 isoforms:</ins>  
![IsoSubset_IFITM3_human_allIsoforms](https://github.com/user-attachments/assets/3bf15ecd-2935-4794-8ffa-dbdae0cf73b6)

<ins>Subset of 5 isoforms:</ins>  
![IsoSubset_IFITM3_human_subset](https://github.com/user-attachments/assets/538eff8c-fc20-46e4-9cff-cad0e6bef376)

## Using BSgenome
BSgenome packages can also be used to plot DNA/peptide sequences. BSgenomes are available for the UCSC and NCBI human genome/annotation, so if you plan to use BSgenome, it is best to align to UCSC or NCBI genome versions.
```
library(BSgenome.Hsapiens.UCSC.hg38)  # Load UCSC BSgenome package

gtf_path <- "hg38.knownGene_ucsc.gtf"
gtf_import(annotation = gtf_path, format = "gtf", organism = "Homo sapiens")  # Load UCSC annotation

uorf_path <- "ENST00000310999_uORF_ucsc.gtf"
eORF_import(annotation = uorf_path, format = "gtf", organism = "Homo sapiens")  # Load uORF coordinates

inputs_full <- create_seq_input(
  rna_files = list(list(plus = "RNAseq_ucsc_plus.bedgraph",
                        minus = "RNAseq_ucsc_minus.bedgraph")),
  ribo_files = list(list(plus = "Riboseq_ucsc_plus.bedgraph",
                         minus = "Riboseq_ucsc_minus.bedgraph")),
  sample_names = "Human ESC",
  rna_paired = c("paired")
)

# To use BSgenome, provide name of BSgenome library in place of FASTA file
ggRibo(tx_id="ENST00000310999.11", eORF.tx_id = "ENST00000310999", NAME = "MRPL11",
       FASTA = BSgenome.Hsapiens.UCSC.hg38, show_seq = T)  # zoomed out (starts/stops)

ggRibo(tx_id = "ENST00000310999.11", eORF.tx_id = "ENST00000310999", NAME = "MRPL11",
                      FASTA = BSgenome.Hsapiens.UCSC.hg38, show_seq = T,
                      plot_range = c(66438710, 66438850), Extend = 0)  # zoomed in (full sequence)

ggRibo_tx(tx_id="ENST00000310999.11", eORF.tx_id = "ENST00000310999", NAME = "MRPL11",
          FASTA = BSgenome.Hsapiens.UCSC.hg38, show_seq = T)  # transcript view
```
