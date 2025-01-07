# Using Ribo-seQC output for ggRibo
[Ribo-seQC](https://github.com/lcalviell/Ribo-seQC) is an R package that performs comprehensive quality analysis of Ribo-seq data.

Among the outputs of Ribo-seQC is two BedGraph files of the plus and minus strand P-sites. The files for the plus and minus strand P-sites match the name of the inputted BAM appended with "_P_sites_plus.bedgraph" and "_P_sites_minus.bedgraph". The BedGraph format contains all the data necessary for ggRibo, but needs to be reformatted.

### Reformatting Ribo-seQC P-site output
The two P-site files can be converted into the format compatible with ggRibo using the following code in R:
```
### Define file paths and load the plus and minus strand BedGraph files
in_plus <- "riboseq_aligned.bam_P_sites_plus.bedgraph"
in_minus <- "riboseq_aligned.bam_P_sites_minus.bedgraph"
out_file <- "riboseq_Psites_combined.ggRibo"

plus <- read.table(in_plus, sep="\t")
minus <- read.table(in_minus, sep="\t")

### Reformat dataframe and output table
plus[,ncol(plus)+1] <- "+"  # assign strands
minus[,ncol(minus)+1] <- "-"

comb <- rbind(plus, minus)  # combine dataframes
comb <- comb[,c(4, 1, 3, 5)]  # reorder columns to match ggRibo format

write.table(comb, file=out_file, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep="\t")
```