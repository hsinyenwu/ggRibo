### Advanced ggRibo part 2: show Ribo-seq reads on overlapping translation events 

Some uORFs could overlap with their main ORFs. Furthermore, ORFs inside the annotated ORFs (nested ORFs) could also be translated. Recently, studies has found downstream ORFs (dPRFs) and dORFs that overlap with the main ORFs (odORFs).       
<img width="535" alt="image" src="https://github.com/user-attachments/assets/2b900fee-888a-40cc-803f-772898a3ed37">

Here we show an example for a gene with uORF and overlapping uORF for AT3G57170 (N-acetylglucosaminyl transferase component family protein). The overlapping uORF on AT3G57170 is conserved in plants:  
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
Code to load RNA-seq and Ribo-seq files, sample names:  
```
CTRL_RNA="~/path/to/RNA_CTRL_merged.bam"
CTRL_ribo="~/path/to/CTRL_expressed_P_sites_sort_count"
FA <- FaFile("~/path/to/TAIR10_chr_all_2.fas")
Samples = c("Seedlings")
RiboseqData = Ribo_data(c(CTRL_ribo),SampleNames=Samples)
RNAseqData = CTRL_RNA
RNAseqBamPairorSingle="paired"
#Load annotated transcript gtf
gtf_import(annotation="~/path/to/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
#Load uORF ouORF gtf
eORF_import(annotation="/path/to/AT3G57170_uORFs.gtf", format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

#Run the plotting function
ggRibo(
  gene_id = "AT3G57170",
  tx_id = "AT3G57170.1",
  eORF.tx_id = "AT3G57170.1",
  NAME="Gpi1 family protein",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/78087fd8-b588-446a-9514-c650d4b49c82)

Plot 2 extra ORFs (one uORF and one ouORF):  
```
ggRibo(
  gene_id = "AT3G57170",
  tx_id = "AT3G57170.1",
  eORF.tx_id = c("AT3G57170.1","AT3G57170.2"),
  NAME="Gpi1 family protein",
  Extend=200)
```
![image](https://github.com/user-attachments/assets/7f07e921-b761-4900-92a3-5442fd700266)

Focus on the first uORF:  
```
ggRibo(
  gene_id = "AT3G57170",
  tx_id = "AT3G57170.1",
  eORF.tx_id = "AT3G57170.2",
  NAME="Gpi1 family protein",
  plot_range = c(21163130,21163250))
```
![image](https://github.com/user-attachments/assets/88d6dcb3-0089-43f4-ab19-db18d1949b78)

Focus on the first uORF and add sequence viewer:  
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
