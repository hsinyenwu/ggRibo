## BedGraph and BigWig for ggRibo

We added the support for bedGraph and BigWig formats in the ggRibo package to load both RNA-seq reads and Ribo-seq P-sites. Users can now leverage these coverage-focused file formats for their analyses. Unlike BAM files, which store detailed alignments for RNA-seq reads, bedGraph and BigWig formats concentrate on representing genome-wide coverage or signal intensity data in a concise and user-friendly manner. **For strand-specific datasets, files are typically separated into two: one for reads aligned to the forward strand and another for the reverse strand.** A particular advantage of BigWig over bedGraph is that it is compressed and indexed, enabling rapid retrieval of specific regions and **reducing both file size and RAM usage for visualization.**  

### 1. RNA-seq bam file to bedGraph and BigWig
Here is an example for how to convert RNA-seq bam files to bedGraph and BigWig formats. In this example there are two samples. You will need samtools, bedtools and bedGraphToBigWig to run the code in linux. Also remember to change the paths for your system.  
***Please modify your path and file names accordingly***
```bash
# Essential: Load samtools and bedtools first
# Define the base path for input and output files (same as the provided code)
PATH1=/path/to/data
cd $PATH1

# Define the list of samples to process
SAMPLES=("Sample1" "Sample2")

# Define CHROM_SIZES file path for the bedGraphToBigWig function.
CHROM_SIZES="/path/to/TAIR10.fas.fai" # CHROM_SIZES, just a fai (FASTA index) file


# Loop over each sample
for sample in "${SAMPLES[@]}"; do
    echo "${sample}"

    BAM="${PATH1}/${sample}_test_PE.bam"

    # ------------------------------------------------------------------
    # 1.  Strand‑specific BAMs that contain **both** reads of each pair
    #     (reverse‑stranded library: Read1– / Read2+  ⇒  plus‑strand mRNA)
    # ------------------------------------------------------------------
    PLUS_BAM="RNA_${sample}_plus_reads.bam"
    MINUS_BAM="RNA_${sample}_minus_reads.bam"

    #  83 (0x53)  = read1 reverse,  read2 forward,  proper pair
    # 163 (0xA3) = read2 reverse,  read1 forward,  proper pair
    samtools view -b -f 83  "$BAM" >  r1rev.bam     # part of plus‑strand transcript
    samtools view -b -f 163 "$BAM" >  r2rev.bam
    samtools merge -f "$PLUS_BAM"  r1rev.bam r2rev.bam
    rm r1rev.bam r2rev.bam

    #  99 (0x63)  = read1 forward, read2 reverse  → minus‑strand transcript
    # 147 (0x93) = read2 forward, read1 reverse
    samtools view -b -f 99  "$BAM" >  r1fwd.bam
    samtools view -b -f 147 "$BAM" >  r2fwd.bam
    samtools merge -f "$MINUS_BAM" r1fwd.bam r2fwd.bam
    rm r1fwd.bam r2fwd.bam

    # -----------------------------------------------------------
    # 2.  Per‑base BedGraph, no scaling, introns removed (‑split)
    # -----------------------------------------------------------
    PLUS_BG="RNA_${sample}_plus.bedgraph"
    MINUS_BG="RNA_${sample}_minus.bedgraph"

    genomeCoverageBed -ibam "$PLUS_BAM"  -bg -split \
      | LC_ALL=C sort -k1,1 -k2,2n > "$PLUS_BG"

    genomeCoverageBed -ibam "$MINUS_BAM" -bg -split \
      | LC_ALL=C sort -k1,1 -k2,2n > "$MINUS_BG"

    # ------------------------------
    # 3.  Convert to strand BigWigs
    # ------------------------------
    /path/to/bedGraphToBigWig "$PLUS_BG"  "$CHROM_SIZES"  "RNA_${sample}_plus.bw"
    /path/to/bedGraphToBigWig "$MINUS_BG" "$CHROM_SIZES"  "RNA_${sample}_minus.bw"

    rm "$PLUS_BAM" "$MINUS_BAM"   # optional cleanup
done
```

### 2. Ribo-seq P-sites from RiboTaper (P_sites_all) to bedGraph and BigWig
Here is just an example with RiboTaper, you can also convert the P-site files from other software to bedGraph and BigWig formats.
```
#Load bedtools bedGraphToBigWig first
# (Optional) If bedGraphToBigWig is available in another module, load it here:
# module load UHTS/Analysis/bedGraphToBigWig/<version>
SAMPLE1=/path/to/P_sites_all_sample1
SAMPLE2=/path/to/P_sites_all_sample2

# Path to your chromosome size file
CHROM_SIZES="/path/to/TAIR10.fas.fai" # CHROM_SIZES, just a fai (FASTA index) file

cd /path/to/data

# Loop over the two samples
for SAMPLE in SAMPLE1 SAMPLE2
do
    FILE=${!SAMPLE}
    ############################################################################
    # 1) Prepare a clean 6-column BED (chrom, start, end, name, score=0, strand)
    ############################################################################
    # If your input file has extra columns, we can cut the first 6 columns to keep it “clean.”
    # Adjust this if your file has a different name or layout.
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,"0",$6}' "$FILE" > ${SAMPLE}.6col.bed

    ############################################################################
    # 2) Split by strand (+ vs -) and sort
    ############################################################################
    awk '$6=="+"' ${SAMPLE}.6col.bed | sort -k1,1 -k2,2n > ${SAMPLE}.plus.bed
    awk '$6=="-"' ${SAMPLE}.6col.bed | sort -k1,1 -k2,2n > ${SAMPLE}.minus.bed

    ############################################################################
    # 3) Convert each strand’s reads to a bedGraph of coverage
    #    (-bg sums the per-position coverage automatically)
    ############################################################################
    genomeCoverageBed \
        -strand + \
        -bg \
        -i ${SAMPLE}.plus.bed \
        -g ${CHROM_SIZES} \
        > ${SAMPLE}.plus.bedGraph

    genomeCoverageBed \
        -strand - \
        -bg \
        -i ${SAMPLE}.minus.bed \
        -g ${CHROM_SIZES} \
        > ${SAMPLE}.minus.bedGraph

    ############################################################################
    # 4) Convert each bedGraph to BigWig
    #    (Requires bedGraphToBigWig and a proper chromosome size file)
    ############################################################################
    bedGraphToBigWig ${SAMPLE}.plus.bedGraph  ${CHROM_SIZES}  ${SAMPLE}.plus.bw
    bedGraphToBigWig ${SAMPLE}.minus.bedGraph ${CHROM_SIZES}  ${SAMPLE}.minus.bw

    # Done with this sample!
    echo "Finished processing $SAMPLE"
    #remove intermediate files
    rm ${SAMPLE}.6col.bed ${SAMPLE}.plus.bed ${SAMPLE}.minus.bed
done
```


### 3. Run the test code for example files
```
#path to annotated gtf
agtf <- system.file("extdata", "TAIR10.29_part.gtf", package = "ggRibo", mustWork = TRUE)
#path to RNA-seq bigwig coverage files
#Remember two files per sample (one for plus strand and one for minus strand)
Root_RNA_p <- system.file("extdata", "Root_plus_strand.bw", package = "ggRibo", mustWork = TRUE)
Root_RNA_m <- system.file("extdata", "Root_minus_strand.bw", package = "ggRibo", mustWork = TRUE)
Shoot_RNA_p <- system.file("extdata", "Shoot_plus_strand.bw", package = "ggRibo", mustWork = TRUE)
Shoot_RNA_m <- system.file("extdata", "Shoot_minus_strand.bw", package = "ggRibo", mustWork = TRUE)
#path to Ribo-seq bigwig coverage files
#Remember two files per sample (one for plus strand and one for minus strand)
Root_Ribo_p <- system.file("extdata", "riboRoot_plus.bw", package = "ggRibo", mustWork = TRUE)
Root_Ribo_m <- system.file("extdata", "riboRoot_minus.bw", package = "ggRibo", mustWork = TRUE)
Shoot_Ribo_p <- system.file("extdata", "riboShoot_plus.bw", package = "ggRibo", mustWork = TRUE)
Shoot_Ribo_m <- system.file("extdata", "riboShoot_minus.bw", package = "ggRibo", mustWork = TRUE)

gtf_import(annotation=agtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
#gtf_import(annotation="~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

#path to the gtf files for uORFs
ugtf <- system.file("extdata", "AT3G02468.gtf", package = "ggRibo", mustWork = TRUE) #uORF gtf
eORF_import(annotation=ugtf, format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

# Load necessary data (assuming Txome_Range and other objects are prepared)
RNA_files <- list(
  list(plus = Root_RNA_p, minus = Root_RNA_m),
  list(plus = Shoot_RNA_p,  minus = Shoot_RNA_m)
)

Ribo_files <- list(
  list(plus = Root_Ribo_p, minus = Root_Ribo_m),
  list(plus = Shoot_Ribo_p,  minus = Shoot_Ribo_m)
)

# Prepare coverage descriptors

inputs_full <- create_seq_input(
  rna_files = RNA_files,
  ribo_files = Ribo_files,
  sample_names = c("Root", "Shoot")
)

ggRibo(gene_id="AT3G02470",tx_id="AT3G02470.3",
          eORF.tx_id = "AT3G02468.1",
          Y_scale="each",Extend=50,
          gene_model_height_ratio=0.8,
          plot_ORF_ranges=T,
          NAME = "SAMDC, CPuORF")

ggRibo(gene_id="AT4G21910",tx_id="AT4G21910.1",
       Y_scale="each",Extend=c(400,50),
       NAME = "MATE efflux family protein")

ggRibo(gene_id="AT4G21910",tx_id="AT4G21910.2",
       Y_scale="each",Extend=c(400,50),
       NAME = "MATE efflux family protein")
```




