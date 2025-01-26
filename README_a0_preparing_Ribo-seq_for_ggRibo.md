## Two ways to generate P-site file (psf) format for ggRibo

**Note: the Ribo-seq P-site coordinate file should look like this:**
The first to forth columns are 
(1) total counts  
(2) chromosome number   
(3) P-site chromsome coordinates  
(4) strand (+ or -)  
**No header in the P-site file please.**
```
1   1  1000000      +
3   1 10000007      +
3   1 10000010      +
3   1 10000016      +
1   1 10000018      +
4   1 10000019      +
```

### 1. Using RiboTaper output
In RiboTaper output files, there is one called *P_sites_all*. You can further process it to the psf format.  
For how to run RiboTaper: see [Here](https://github.com/hsinyenwu/Riboseq_pipeline/blob/main/README_pipeline.md#9-orf-visualization-ggribo).  

```
## Format and copy files needed for ggRibo
# Run this code in linux command line.
# ${taper_out} contains the path of output P-site file.
# $p_site_file contains the path and file name of output P-site file.
cut -f 1,3,6 ${taper_out}/P_sites_all | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' > $p_site_file 
```

### 2. Using Ribo-seQC output for ggRibo
Ribo-seQC is an R package that performs comprehensive quality analysis of Ribo-seq data. For R version >4.0, please download a updated version of Ribo-seQC [Here](https://github.com/hsinyenwu/RiboseQC_R4.2.1). For older version of R, download Ribo-seQC [here](https://github.com/lcalviell/Ribo-seQC).

Among the outputs of Ribo-seQC is two BedGraph files of the plus and minus strand P-sites. The files for the plus and minus strand P-sites match the name of the inputted BAM appended with "_P_sites_plus.bedgraph" and "_P_sites_minus.bedgraph". The BedGraph format contains all the data necessary for ggRibo, but needs to be reformatted.

**Reformatting Ribo-seQC P-site output**
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
