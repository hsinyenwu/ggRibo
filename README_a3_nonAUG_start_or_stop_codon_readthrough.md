## Advanced ggRibo part 3: visualizing dORFs and translation readthrough

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


