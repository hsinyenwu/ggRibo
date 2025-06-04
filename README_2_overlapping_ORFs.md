### Advanced ggRibo part 2: show Ribo-seq reads on overlapped translation events 
#### *(Also deconvolute the reading frames for overlapping translation events)* 

Some uORFs could overlap with their main ORFs. Furthermore, ORFs inside the annotated ORFs (nested ORFs) could also be translated. Recently, studies has found downstream ORFs (dPRFs) and dORFs that overlap with the main ORFs (odORFs).       
<img width="535" alt="image" src="https://github.com/user-attachments/assets/2b900fee-888a-40cc-803f-772898a3ed37">

Here we show an example for a gene with uORF and overlapping uORF for AT3G57170 (N-acetylglucosaminyl transferase component family protein). The overlapping uORF on AT3G57170 is conserved in plants.

Here is the gtf for the uORF and overlapping uORF, you can copy this file and save as AT3G57170_uORFs.gtf:  
```
3	Araport11	gene	21159373	21163322	.	-	.	gene_id AT3G57170; gene_biotype protein_coding;
3	Araport11	mRNA	21159373	21163322	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21159373	21159731	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21159824	21159917	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21160026	21160435	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21160688	21161182	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21161420	21161641	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21161731	21161837	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21161924	21162064	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21162158	21162524	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21162619	21162798	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	exon	21163121	21163322	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.1; gene_biotype protein_coding;
3	Araport11	CDS	21162619	21162789	.	-	.	gene_id "AT3G57170"; transcript_id "AT3G57170.1"; gene_biotype "protein_coding";
3	Araport11	CDS	21162158	21162524	.	-	.	gene_id "AT3G57170"; transcript_id "AT3G57170.1"; gene_biotype "protein_coding";
3	Araport11	CDS	21162051	21162064	.	-	.	gene_id "AT3G57170"; transcript_id "AT3G57170.1"; gene_biotype "protein_coding";
3	Araport11	mRNA	21159373	21163322	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21159373	21159731	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21159824	21159917	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21160026	21160435	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21160688	21161182	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21161420	21161641	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21161731	21161837	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21161924	21162064	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21162158	21162524	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21162619	21162798	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	exon	21163121	21163322	.	-	.	gene_id AT3G57170; transcript_id AT3G57170.2; gene_biotype protein_coding;
3	Araport11	CDS	21163146	21163187	.	-	.	gene_id "AT3G57170"; transcript_id "AT3G57170.2"; gene_biotype "protein_coding";
```
The reason that we need to use two transcript ID (AT3G57170.1 and AT3G57170.2) for the two extra ORFs is because gtf/gff files only allow one CDS per transcript. Remember that rule when you create your own gtf files for eORFs (extra ORFs other than annotated main ORFs). That is also the reason that we do not prefer to create one gtf for main ORF gtf and uORF gtf.  

You can also load the above file with:  
```
ugtf <- system.file("extdata", "AT3G57170_uORFs.gtf", package = "ggRibo", mustWork = TRUE) #uORF gtf
eORF_import(annotation=ugtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```
**Code to load RNA-seq and Ribo-seq files, sample names:**  
```
CTRL_RNA="~/path/to/RNA_CTRL_merged.bam"
CTRL_ribo="~/path/to/CTRL_expressed_P_sites_sort_count"
FA <- FaFile("~/path/to/TAIR10_chr_all_2.fas")

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
#Load uORF ouORF gtf
eORF_import(annotation="/path/to/AT3G57170_uORFs.gtf", format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
```
**Plot the overlapping uORF:**  
```
#Run the plotting function
ggRibo(
  tx_id = "AT3G57170.1",
  eORF.tx_id = "AT3G57170.1",
  NAME="Gpi1 family protein",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/78087fd8-b588-446a-9514-c650d4b49c82)

**Plot 2 extra ORFs (one uORF and one ouORF):**  
```
ggRibo(
  tx_id = "AT3G57170.1",
  eORF.tx_id = c("AT3G57170.1","AT3G57170.2"),
  NAME="Gpi1 family protein",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/7f07e921-b761-4900-92a3-5442fd700266)

**Focus on the first uORF:**  
```
ggRibo(
  tx_id = "AT3G57170.1",
  eORF.tx_id = "AT3G57170.2",
  NAME="Gpi1 family protein",
  plot_range = c(21163130,21163250))
```
![image](https://github.com/user-attachments/assets/88d6dcb3-0089-43f4-ab19-db18d1949b78)

**Focus on the first uORF and add sequence viewer:**  
```
ggRibo(
  gene_id = "AT3G57170",
  tx_id = "AT3G57170.1",
  eORF.tx_id = "AT3G57170.2",
  NAME="Gpi1 family protein",
  show_seq = TRUE,FASTA = FA,
  plot_range = c(21163130,21163250))
```
![image](https://github.com/user-attachments/assets/959f8d9b-1824-4e98-b5fb-ee484263314b)

**Show the ouORF frame relative to the main ORF**  
```
ggRibo(
  tx_id = "AT3G57170.1",
  eORF.tx_id = "AT3G57170.1",
  NAME="Gpi1 family protein",
  oORF_coloring = "extend_mORF",
  Extend=200)
#oORF_coloring = "extend_mORF" does not do anything if the eORF is not overlapping with the main ORF.
```
![image](https://github.com/user-attachments/assets/455b64ab-e5ad-43dd-8060-b3ad96f76f05)
**Show the ouORF frame alone**  
```
ggRibo(
  tx_id = "AT3G57170.1",
  eORF.tx_id = "AT3G57170.1",
  NAME="Gpi1 family protein",
  oORF_coloring = "oORF_colors",
  gene_model_height_ratio = 0.3,
  Extend=200)
```
![image](https://github.com/user-attachments/assets/5f90a8c1-47ab-45e8-bc88-497f0a3e5d11)

### Decomposition of Ribo-seq reads shows the translated frames of overlapping uORF and main ORF  
ggRibo_decon takes only one Ribo-seq and RNA-seq samples and plot the 3 frames separately.   
**If we only assign frame colors to the annotated CDS**    
```
ggRibo_decom(
  tx_id = "AT3G57170.1",
  NAME="Gpi1 family protein",
  oORF_coloring = "extend_mORF",
  frame_logic ="CDS_start",
  gene_model_height_ratio = 0.7,
  Extend=200)
```
![image](https://github.com/user-attachments/assets/7d8abbfb-b1d3-4e1b-b400-cb5ceb0dbcd8)

**The above "frame_logic" parameter has 3 options:**
1. frame_logic="tx_start", the frame 0 starts from the beginning of the transcript. This is the default for ncRNA.  
2. frame_logic="CDS_start",  the frame 0 starts from the beginning of the annotated CDS. This is the default for coding transcript.   
3. frame_logic="CDS_extend", the frame 0 starts from the beginning of the annotated CDS and extend to the two ends of the transcript.  

**We can try frame_logic = "CDS_extend")**   
Now the eORF ranges are shown, but do not guide the coloring of eORF ribo-seq reads.
```
ggRibo_decom(
    gene_id = "AT3G57170",
    tx_id = "AT3G57170.1",
    eORF.tx_id = c("AT3G57170.1","AT3G57170.2"),
    NAME="Gpi1 family protein",gene_model_height_ratio = 0.7,
    Extend=200,frame_logic = "CDS_extend")
```
![image](https://github.com/user-attachments/assets/99895619-e136-4847-8b46-0231aea45504)

### Single transcript view of oORF
From the above data, we can see isoform 1 is expressed for the above gene, so we can use ggRibo_tx to see the uORF and ouORFs.
```
ggRibo_tx(
    tx_id = "AT3G57170.1",
    eORF.tx_id = "AT3G57170.1",
    NAME="Gpi1 family protein",
    oORF_coloring = "extend_mORF", #default
    gene_model_height_ratio =1.8)
```
![image](https://github.com/user-attachments/assets/13ff7de6-0412-4c8d-9317-2bd2e66a2c0e)
To see only uORF and ouORF (their frame coloring is based on their own frames):
```
ggRibo_tx(
    tx_id = "AT3G57170.1",
    eORF.tx_id = c("AT3G57170.1","AT3G57170.2"),
    NAME="Gpi1 family protein",
    gene_model_height_ratio =1.2,
    oORF_coloring = "oORF_colors")
```
![image](https://github.com/user-attachments/assets/8ffc7082-01b2-4abc-bbdf-2abcbf8ed7b9)
To see ouORF (its frame coloring scheme is extended from the main ORF):
```
ggRibo_tx(
    tx_id = "AT3G57170.1",
    eORF.tx_id = c("AT3G57170.1","AT3G57170.2"),
    NAME="Gpi1 family protein",
    gene_model_height_ratio =1.2,
    oORF_coloring = "extend_mORF")
```
![image](https://github.com/user-attachments/assets/1a69f318-bb31-4d78-a34e-b08d28458bd9)





