## Advanced ggRibo part 1: plot Ribo-seq reads on annotated noncoding RNAs

### Reading frame coloring for ncRNA:
For noncoding gene or noncoding isoforms, their ribo-seq coloring for reading frame is different from coding ORFs.  
**Frame 0 for coding ORFs:** start from the first nucleotide of the annotated ORF/CDS. So the translated ORFs should be red.   
**Frame 0 for noncoding RNA:** start from the first nucleotide of the transcript. So the translated ORFs could be red, blue or green. The enrichment of a color in an area suggest a potential coding ORF/CDS.  

### Other notes:
* No ORF/CDS border dash lines is plotted for ncORFs.  
* The 3 frames of one ORF are colored consistently even cross intron (e.g., if frame 0 is blue in exon 1, it will still be blue in other exons).    

### Plot a noncoding gene
```
ggRibo(
  gene_id = "AT3G17185",
  tx_id = "AT3G17185.1",
  NAME="TAS3",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/d452f80e-7703-4295-85e4-b64860ca1e0c)

### Plot a gene with both coding or noncoding isoforms
To make an example: I artifically removed the CDS for both *AT1G01060.5* and *AT1G01060.7* transcript for the *AT1G01060 (LHY)* gene.  
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


