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
_Insert uORF1 here_

## Add genome/peptide sequence
```
fasta <- FaFile("Homo_sapiens.GRCh38.dna.primary_assembly.fa")

ggRibo(gene_id = "ENSG00000174547", tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999",
       NAME = "MRPL11", FASTA = fasta, show_seq = T)
```
_Insert uORF2 here_

## Closer look at uORF (adjusting plot range)
```
range <- c(66438710, 66438850)
ggRibo(gene_id = "ENSG00000174547", tx_id = "ENST00000310999", eORF.tx_id = "ENST00000310999",
       NAME = "MRPL11", FASTA = fasta, show_seq = T,
       plot_range = range, Extend = 0)
```
_Insert uORF3 here_

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

# Higher translation in Cardiomyocytes
ggRibo(gene_id = "ENSG00000106211", tx_id = "ENST00000248553", NAME = "HSPB1")

# Higher translation in iPSCs
ggRibo(gene_id = "ENSG00000133112", tx_id = "ENST00000530705", NAME = "TPT1")
```
_Insert TwoS1 here_

_Insert TwoS2 here_

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

gene <- "ENSG00000142089"
tx <- "ENST00000399808"
name <- "IFITM3"

### Plot all 27 isoforms of IFITM3
ggRibo(gene_id = gene, tx_id = tx, NAME = name)

### Plot subset
sub_Txome <- Txome_Range$clone()  # make copy of Txome_Range
tx_subset <- c("ENST00000399808", "ENST00000680209", "ENST00000681198",
               "ENST00000681304", "ENST00000681840")  # define list of desired isoforms
sub_Txome$txByGene[[gene]] <- sub_Txome$txByGene[[gene]][sub_Txome$txByGene[[gene]]$tx_name %in% tx_subset]

ggRibo(gene_id = gene, tx_id = tx, NAME = name, GRangeInfo = sub_Txome)  # use subsetted Txome object

```
_Insert Iso1 here_

_Insert Iso2 here_
