# ---------------------------
# Load Necessary Packages
# ---------------------------

# Suppress warnings and package startup messages when loading libraries
suppressWarnings({
  suppressPackageStartupMessages({
    # GenomicRanges: Representation and manipulation of genomic intervals
    library(GenomicRanges)
    # GenomicFeatures: Tools for working with genomic annotation data
    library(GenomicFeatures)
    # GenomicAlignments: Efficient manipulation of large-scale sequencing data
    library(GenomicAlignments)
    # Rsamtools: Interface to the SAM/BAM sequence alignment format
    library(Rsamtools)
    # ggplot2: Data visualization package
    library(ggplot2)
    # cowplot: For assembling multiple ggplot2 plots into one
    library(cowplot)
    # grid: Low-level graphics functions
    library(grid)
    # IRanges: Representation and manipulation of integer ranges
    library(IRanges)
    # dplyr: Data manipulation package
    library(dplyr)
    # txdbmaker: Tools for creating transcript databases
    library(txdbmaker)
    # R6: Implementation of classes with reference semantics
    library(R6)
    # BiocParallel: Support for parallel computation in Bioconductor packages
    library(BiocParallel)
  })
})

# ---------------------------
# Class Definitions (R6 Classes)
# ---------------------------

# R6 Class to store range information from the annotation
#' @title Range_info Class
#' @description
#' An R6 class to store range information from the genomic annotation.
#'
#' @field exonsByTx Genomic ranges of exons by transcript.
#' @field txByGene Transcripts grouped by gene.
#' @field cdsByTx Coding sequences grouped by transcript.
#' @field fiveUTR Five prime untranslated regions grouped by transcript.
#' @field threeUTR Three prime untranslated regions grouped by transcript.
Range_info <- R6::R6Class("Range_info",
  public = list(
    exonsByTx = NULL,      # Exons grouped by transcript
    txByGene = NULL,       # Transcripts grouped by gene
    cdsByTx = NULL,        # Coding sequences grouped by transcript
    fiveUTR = NULL,        # 5' UTR regions grouped by transcript
    threeUTR = NULL,       # 3' UTR regions grouped by transcript
    
    # Initialization method for Range_info
    initialize = function(exonsByTx, txByGene, cdsByTx, fiveUTR, threeUTR) {
      # Assign provided genomic ranges to class fields
      self$exonsByTx <- exonsByTx
      self$txByGene <- txByGene
      self$cdsByTx <- cdsByTx
      self$fiveUTR <- fiveUTR
      self$threeUTR <- threeUTR
    }
  )
)

# R6 Class to store eORF (extra Open Reading Frame) range information
#' @title eORF_Range_info Class
#' @description
#' An R6 class to store eORF (extra Open Reading Frame) range information.
#'
#' @field eORFByTx Genomic ranges of eORFs by transcript.
eORF_Range_info <- R6::R6Class("eORF_Range_info",
  public = list(
    eORFByTx = NULL,      # eORFs grouped by transcript
    
    # Initialization method for eORF_Range_info
    initialize = function(eORFByTx) {
      # Assign provided eORF ranges to class field
      self$eORFByTx <- eORFByTx
    }
  )
)

# R6 Class to store gene-specific information
#' @title Gene_info Class
#' @description
#' An R6 class to store gene-specific information.
#'
#' @field gene_id Gene identifier.
#' @field tx_id Transcript identifier.
#' @field txByGene Transcripts grouped by gene.
#' @field cdsByYFGtx Coding sequences of the gene's transcripts.
#' @field chr Chromosome information.
#' @field generanges Genomic ranges for the gene.
#' @field generangesplus Extended genomic ranges for the gene.
#' @field range_left Left boundary of the gene range.
#' @field range_right Right boundary of the gene range.
#' @field num_isoforms Number of isoforms.
#' @field tx_names Names of transcripts.
#' @field isoforms.w.3UTR Isoforms with 3' UTR.
#' @field isoforms.w.5UTR Isoforms with 5' UTR.
#' @field threeUTRByYFGtx 3' UTRs by transcript.
#' @field fiveUTRByYFGtx 5' UTRs by transcript.
#' @field exonByYFGtx Exons by transcript.
#' @field Extend Extension length.
#' @field strand Strand information.
#' @field xlimCds CDS ranges for plotting.
#' @field Riboseq_list List of Ribo-seq data.
#' @field cds_left Left boundary of CDS.
#' @field cds_right Right boundary of CDS.
Gene_info <- R6::R6Class("Gene_info",
  public = list(
    gene_id = NULL,               # Identifier for the gene
    tx_id = NULL,                 # Main transcript ID
    txByGene = NULL,              # Transcripts grouped by gene
    cdsByYFGtx = NULL,            # Coding sequences grouped by transcript
    chr = NULL,                   # Chromosome name
    generanges = NULL,            # Genomic ranges for the gene
    generangesplus = NULL,        # Genomic ranges with additional information
    range_left = NULL,            # Left boundary of the gene range
    range_right = NULL,           # Right boundary of the gene range
    num_isoforms = NULL,          # Number of transcript isoforms
    tx_names = NULL,              # Names of transcripts
    isoforms.w.3UTR = NULL,       # Isoforms with 3' UTR
    isoforms.w.5UTR = NULL,       # Isoforms with 5' UTR
    threeUTRByYFGtx = NULL,       # 3' UTR regions grouped by transcript
    fiveUTRByYFGtx = NULL,        # 5' UTR regions grouped by transcript
    exonByYFGtx = NULL,           # Exons grouped by transcript
    Extend = NULL,                # Extension parameter for plotting
    strand = NULL,                # Strand information ('+' or '-')
    xlimCds = NULL,               # CDS limits for plotting
    Riboseq_list = NULL,          # List of Ribo-seq data
    cds_left = NULL,              # Left boundary of CDS
    cds_right = NULL,             # Right boundary of CDS
    
    # Initialization method for Gene_info
    initialize = function(gene_id, tx_id, txByGene, cdsByYFGtx, chr, generanges, generangesplus,
                          range_left, range_right, num_isoforms, tx_names, isoforms.w.3UTR,
                          isoforms.w.5UTR,  # Added this parameter
                          threeUTRByYFGtx, fiveUTRByYFGtx, exonByYFGtx, Extend, strand,
                          xlimCds, Riboseq_list, cds_left, cds_right) {
      # Assign provided gene-specific data to class fields
      self$gene_id <- gene_id
      self$tx_id <- tx_id
      self$txByGene <- txByGene
      self$cdsByYFGtx <- cdsByYFGtx
      self$chr <- chr
      self$generanges <- generanges
      self$generangesplus <- generangesplus
      self$range_left <- range_left
      self$range_right <- range_right
      self$num_isoforms <- num_isoforms
      self$tx_names <- tx_names
      self$isoforms.w.3UTR <- isoforms.w.3UTR
      self$isoforms.w.5UTR <- isoforms.w.5UTR  # Assign the additional parameter
      self$threeUTRByYFGtx <- threeUTRByYFGtx
      self$fiveUTRByYFGtx <- fiveUTRByYFGtx
      self$exonByYFGtx <- exonByYFGtx
      self$Extend <- Extend
      self$strand <- strand
      self$xlimCds <- xlimCds
      self$Riboseq_list <- Riboseq_list
      self$cds_left <- cds_left
      self$cds_right <- cds_right
    }
  )
)

# R6 Class to store eORF-specific information
#' @title eORF_info Class
#' @description
#' An R6 class to store eORF-specific information.
#'
#' @field eORF.tx_id eORF transcript identifiers.
#' @field eORF_Riboseq_list List of Ribo-seq data for eORFs.
#' @field xlim.eORF eORF ranges for plotting.
#' @field eORF_left Left boundaries of eORFs.
#' @field eORF_right Right boundaries of eORFs.
eORF_info <- R6::R6Class("eORF_info",
  public = list(
    eORF.tx_id = NULL,           # eORF transcript IDs
    eORF_Riboseq_list = NULL,    # List of Ribo-seq data for eORFs
    xlim.eORF = NULL,            # eORF genomic range limits
    eORF_left = NULL,            # Left boundary of eORF
    eORF_right = NULL,           # Right boundary of eORF
    
    # Initialization method for eORF_info
    initialize = function(eORF.tx_id, eORF_Riboseq_list, xlim.eORF, eORF_left, eORF_right) {
      # Assign provided eORF-specific data to class fields
      self$eORF.tx_id <- eORF.tx_id
      self$eORF_Riboseq_list <- eORF_Riboseq_list
      self$xlim.eORF <- xlim.eORF
      self$eORF_left <- eORF_left
      self$eORF_right <- eORF_right
    }
  )
)
