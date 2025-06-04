## Advanced ggRibo part 1: plot Ribo-seq reads on annotated noncoding RNAs

### Reading frame coloring for ncRNA:
For noncoding gene or noncoding isoforms, their ribo-seq coloring for reading frame is different from coding ORFs.  
* **Frame 0 for coding ORFs:** start from the first nucleotide of the annotated ORF/CDS. So the translated ORFs should be red.   
* **Frame 0 for noncoding RNA:** start from the first nucleotide of the transcript. So the translated ORFs could be red, blue or green. The enrichment of a color in an area suggest a potential coding ORF/CDS.  

### Other notes:
* No ORF/CDS border dash lines is plotted for ncORFs.  
* The 3 frames of one ORF are colored consistently even cross intron (e.g., if frame 0 is blue in exon 1, it will still be blue in other exons).    

### Plot an annotated noncoding gene TAS3
```
#Load data first
CTRL_RNA="~/path/to/RNA_CTRL_merged.bam"
CTRL_ribo="~/path/to/CTRL_expressed_P_sites_sort_count"

RNA_files <- list(CTRL_RNA)
Ribo_files <- list(CTRL_ribo)
Samples <- c("Seedlings")

# Prepare coverage descriptors
inputs_full <- create_seq_input(
  rna_files = RNA_files,
  ribo_files = Ribo_files,
  sample_names = Samples,
)

#Load annotated transcript gtf
gtf_import(annotation="~/path/to/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

ggRibo(
  gene_id = "AT3G17185",
  tx_id = "AT3G17185.1",
  NAME="TAS3",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/d452f80e-7703-4295-85e4-b64860ca1e0c)

Check the reads in 3 frames with ggRibo_decom (decom for decomposition).
```
ggRibo_decom(gene_id = "AT3G17185",
             tx_id = "AT3G17185.1",
             NAME="TAS3",
             plot_genomic_direction = TRUE,
             Extend=50)
```
![image](https://github.com/user-attachments/assets/b5be24a2-4289-4d32-992f-704c9d5a5752)

### Plot another annotated noncoding gene sORF1 (from Hsu et al., PNAS 2016)
As mentioned above, for non-coding RNAs, ggRibo assigns the reading frame from the first nucleotide of the annotated RNA sequence, rather than from the start of a CDS/ORF, as it does for coding RNAs. As a result, a translated ORF in a non-coding RNA may enrich one of the reading frames in red, blue, or green. For sORF1 (below), the main translated ORF is colored green. However, you can still provide a gtf with annotated ORF ranges for visualizing this sORF. 

```
ggRibo(gene_id = "AT1G10682",
       tx_id = "AT1G10682.1",
       NAME="sORF1",
       Extend=50)
```
![image](https://github.com/user-attachments/assets/a17af73c-0b76-4454-8255-7ac0e0cfc8e8)

Ribo-seq reads decomposition for frame enrichment:
```
ggRibo_decom(gene_id = "AT1G10682",
             tx_id = "AT1G10682.1",
             NAME="sORF1",
             Extend=50)
```
![image](https://github.com/user-attachments/assets/c3658d18-e6bf-44b3-a4e3-59eb0868525f)

### Plot a gene with both coding or noncoding isoforms
The situation: **a gene with both coding and noncoding isoforms** occurs a lot in later versions of animal and plant annotations. Here we show that ggRibo can still plot those genes.  
To make an example: I artifically removed the CDSs for both *AT1G01060.5* and *AT1G01060.7* transcripts for the *AT1G01060 (LHY)* gene.  
Plot the noncoding isoform *AT1G01060.7*  
```
tgtf <- system.file("extdata", "AT1G01060_test.gtf", package = "ggRibo", mustWork = TRUE)
gtf_import(annotation=tgtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
ggRibo(
  gene_id = "AT1G01060",
  tx_id = "AT1G01060.7",
  NAME="",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/082cc0f1-4ffc-48ab-aa14-e291e59e1ee4)

Plot the noncoding isoform *AT1G01060.7* with ggRibo_decom for frames enriched. 
```
ggRibo_decom(
  gene_id = "AT1G01060",
  tx_id = "AT1G01060.7",
  NAME="",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/7b837b75-bbcc-4861-90b0-85490b7b5670)

Plot a coding isoform for LHY.
```
ggRibo(
  gene_id = "AT1G01060",
  tx_id = "AT1G01060.4",
  NAME="",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/4a148431-3ede-4125-9b3e-247c62a62bc7)


