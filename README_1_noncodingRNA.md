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
As mentioned above, for non-coding RNAs, ggRibo assigns the reading frame from the first nucleotide of the annotated RNA sequence, rather than from the start of a CDS/ORF, as it does for coding RNAs. As a result, a translated ORF in a non-coding RNA may enrich one of the reading frames in red, blue, or green. For sORF1 (below), the main translated ORF is colored green.

```
ggRibo(gene_id = "AT1G10682",
       tx_id = "AT1G10682.1",
       NAME="sORF1",
       Extend=50)
```
![image](https://github.com/user-attachments/assets/a17af73c-0b76-4454-8255-7ac0e0cfc8e8)

Ribo-seq reads decomposition:
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
ggRibo(
  gene_id = "AT1G01060",
  tx_id = "AT1G01060.5",
  NAME="",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/d62096ae-79bf-421f-a15f-42dc6de4667f)
Plot the noncoding isoform *AT1G01060.5*  
```
ggRibo(
  gene_id = "AT1G01060",
  tx_id = "AT1G01060.7",
  NAME="",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/fc99ea78-e5f9-4b64-9cd4-cd4f60b81ef4)  

Plot a coding isoform for LHY.
```
ggRibo(
  gene_id = "AT1G01060",
  tx_id = "AT1G01060.4",
  NAME="",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/c3b7a1a4-0520-42b9-b130-1795735385c7)

