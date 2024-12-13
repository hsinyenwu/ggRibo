## Advanced ggRibo part 3: visualizing dORF translation
Some ORFs are identified downstream of the main ORF. Here we visualize a case of dORF in Arabidopsis.  

**gtf for dORF identified from PMID: 38000896**
You can copy and paste in your computer to make the *AT1G27950_dORF.gtf* file.   
```
1	Araport11	gene	9740508	9742159	.	+	.	gene_id AT1G27950; gene_biotype protein_coding;
1	Araport11	mRNA	9740508	9742159	.	+	.	gene_id AT1G27950; transcript_id AT1G27950.1; gene_biotype protein_coding;
1	Araport11	exon	9740508	9741091	.	+	.	gene_id AT1G27950; transcript_id AT1G27950.1; gene_biotype protein_coding;
1	Araport11	exon	9741762	9742159	.	+	.	gene_id AT1G27950; transcript_id AT1G27950.1; gene_biotype protein_coding;
1	Araport11	CDS	9742001	9742114	.	+	0	gene_id AT1G27950; transcript_id AT1G27950.1; gene_biotype protein_coding;
```

Read in the dORF gtf and plot with ggRibo.   
```
eORF_import(annotation="~/Desktop/AT1G27950_dORF.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
ggRibo(gene_id="AT1G27950",
       tx_id = "AT1G27950.1",
       eORF.tx_id="AT1G27950.1",
       Extend = 300)
```
Almost cannot see it. The dORF must be very poorly translated comparing to the main ORF.    
![image](https://github.com/user-attachments/assets/7e51defc-f4d1-4340-812d-e9e03c2352ca)

But you can set the height (Y-axis value) for the Ribo-seq reads in the plots with ***Ribo_fix_height***. Here I set the height of Ribo-seq to 8, which means the maximum of Y-axis for Ribo-seq is 8 reads.  
```
ggRibo(gene_id="AT1G27950",
       tx_id = "AT1G27950.1",
       eORF.tx_id="AT1G27950.1",
       Ribo_fix_height = 12,
       Extend = 300)
```
![image](https://github.com/user-attachments/assets/6eccf63b-4ef2-4794-b3a4-e3236a350e5f)

We can futher zoom in the dORF and see it more clearly.  
```
ggRibo(gene_id="AT1G27950",
       tx_id = "AT1G27950.1",
       eORF.tx_id="AT1G27950.1",
       plot_range = c(9741980,9742125),
       show_seq = TRUE,FASTA = FA,
       Ribo_fix_height = 8,
       Extend = 300)
```
![image](https://github.com/user-attachments/assets/9b3bce91-3ff3-4121-8ccc-d252d475caee)



