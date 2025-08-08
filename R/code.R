# ---------------------------
# Function Definitions
# ---------------------------

# The following functions are used for importing annotation data, reading Ribo-seq data,
# and assigning reading frames to sequencing reads. They collectively help in preparing
# genomic annotations and Ribo-seq datasets for downstream analyses and plotting.

# Function to import GTF/GFF annotation and create a Range_info object
#' @title Import GTF/GFF Annotation
#' @description
#' Imports GTF/GFF annotation files and creates a \code{Range_info} object.
#'
#' @param annotation Path to the annotation file.
#' @param format Format of the annotation file ("gtf" or "gff").
#' @param dataSource Optional data source description.
#' @param organism Optional organism name.
#' @return A \code{Range_info} object stored in the global environment as \code{Txome_Range}.
gtf_import <- function(annotation, format = "gtf", dataSource = "", organism = "") {
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(annotation, format = format, dataSource = dataSource, organism = organism))
  exonsByTx <- exonsBy(txdb, by = 'tx', use.names = TRUE)
  txByGene <- transcriptsBy(txdb, by = 'gene')
  cdsByTx <- cdsBy(txdb, by = "tx", use.names = TRUE)
  fiveUTR <- fiveUTRsByTranscript(txdb, use.names = TRUE)
  threeUTR <- threeUTRsByTranscript(txdb, use.names = TRUE)
  
  # Create a transcript-to-gene lookup table efficiently
  tx_to_gene <- AnnotationDbi::select(txdb, keys = keys(txdb, keytype = "TXNAME"), 
                       columns = c("TXNAME", "GENEID"), 
                       keytype = "TXNAME")
  colnames(tx_to_gene) <- c("tx_id", "gene_id")
  
  # Sort the lookup table by tx_id (smallest to largest)
  tx_to_gene <- tx_to_gene[order(tx_to_gene$tx_id), ]
  
  # Create Range_info object
  Txome_Range <- Range_info$new(
    exonsByTx = exonsByTx,
    txByGene = txByGene,
    cdsByTx = cdsByTx,
    fiveUTR = fiveUTR,
    threeUTR = threeUTR,
    tx_to_gene = tx_to_gene
  )
  assign("Txome_Range", Txome_Range, envir = .GlobalEnv)
}

# Function to import eORF annotation and create an eORF_Range_info object
#' @title Import eORF Annotation
#' @description
#' Imports eORF annotation files and creates an \code{eORF_Range_info} object.
#'
#' @param annotation Path to the eORF annotation file.
#' @param format Format of the annotation file ("gtf" or "gff").
#' @param dataSource Optional data source description.
#' @param organism Optional organism name.
#' @return An \code{eORF_Range_info} object stored in the global environment as \code{eORF_Range}.
eORF_import <- function(annotation, format = "gtf", dataSource = "", organism = "") {
  # Create a TxDb object from the eORF annotation file
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(annotation, format = format, dataSource = dataSource, organism = organism))

  # Extract CDS ranges by transcript, which correspond to eORFs
  cdsByTx <- cdsBy(txdb, by = "tx", use.names = TRUE)

  # Create the eORF_Range_info object containing extracted eORF ranges
  eORF_Range <- eORF_Range_info$new(
    eORFByTx = cdsByTx
  )

  # Assign the eORF_Range_info object to the global environment
  assign("eORF_Range", eORF_Range, envir = .GlobalEnv)
}

# Function to read in Ribo-seq data files
#' @title Read Ribo-seq Data Files
#' @description
#' Reads in Ribo-seq data files and returns a list of data frames.
#'
#' @param RiboseqData Vector of file paths to Ribo-seq data files.
#' @param SampleNames Vector of sample names corresponding to the data files.
#' @return A named list of data frames containing Ribo-seq data.
Ribo_data <- function(RiboseqData, SampleNames) {
  # Read each Ribo-seq data file into a list
  Ribo_data_list <- lapply(RiboseqData, function(file) {
    # Read the file as tab-delimited, without headers
    Ribo1 <- read.delim(file = file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    # Assign meaningful column names: count, chr, position, strand
    colnames(Ribo1) <- c("count", "chr", "position", "strand")
    # Return the processed data frame
    Ribo1
  })
  # Name each list element according to the sample name for easy identification
  names(Ribo_data_list) <- SampleNames
  return(Ribo_data_list)
}

# Helper Function: Assign frames with extended CDS ranges
#' @title Assign Frames with Extended CDS Ranges
#' @description
#' Assigns reading frames to Ribo-seq reads based on extended CDS ranges, including overlapping eORFs.
#'
#' @param Ribo_data Data frame containing Ribo-seq reads.
#' @param extended_cds_ranges Extended CDS ranges as a \code{GRanges} object.
#' @param strand Strand information ("+" or "-").
#' @param main_cds_ranges Main CDS ranges as a \code{GRanges} object.
#' @return Data frame with an added \code{frame} column indicating the reading frame.
assign_frames_extended <- function(Ribo_data, extended_cds_ranges, strand, main_cds_ranges) {
  # Check if extended CDS ranges are present
  if (length(extended_cds_ranges) > 0) {
    # Sort exons depending on strand to ensure proper ordering
    if (strand == "+") {
      exons <- sort(extended_cds_ranges, decreasing = FALSE)
    } else {
      exons <- sort(extended_cds_ranges, decreasing = TRUE)
    }

    # Initialize vectors to store genomic positions and corresponding transcript positions
    positions <- integer(0)
    tx_positions <- integer(0)
    cum_len <- 0

    # Map each exon’s genomic positions to transcript-based positions
    for (idx in seq_along(exons)) {
      exon <- exons[idx]
      pos <- seq(start(exon), end(exon))
      # Reverse if on negative strand
      if (strand == "-") {
        pos <- rev(pos)
      }
      len <- length(pos)
      # Assign consecutive transcript positions
      tx_pos <- seq_len(len) + cum_len
      # Accumulate positions
      positions <- c(positions, pos)
      tx_positions <- c(tx_positions, tx_pos)
      cum_len <- cum_len + len
    }

    # Create a data frame linking genomic positions to transcript positions
    position_df <- data.frame(position = positions, tx_pos = tx_positions)

    # Determine main ORF start position in genomic coordinates
    if (strand == "+") {
      main_orf_start <- min(start(main_cds_ranges))
    } else {
      main_orf_start <- max(end(main_cds_ranges))
    }

    # Find transcript position of the main ORF start
    main_orf_tx_pos <- position_df$tx_pos[position_df$position == main_orf_start][1]

    # Compute frames relative to the main ORF start
    position_df$frame <- (position_df$tx_pos - main_orf_tx_pos) %% 3
    position_df$frame <- factor(position_df$frame, levels = c(0,1,2))

    # Merge the frame information into the Ribo_data
    Ribo_data <- merge(Ribo_data, position_df[, c("position", "frame")], by = "position", all.x = TRUE)
  } else {
    # If no extended CDS ranges, assign NA frames
    Ribo_data$frame <- factor(NA, levels = c(0,1,2))
  }

  return(Ribo_data)
}

# Helper function to assign frames to Ribo-seq reads
#' @title Assign Frames to Ribo-seq Reads
#' @description
#' Assigns reading frames to Ribo-seq reads based on CDS or eORF ranges.
#'
#' @param Ribo_data Data frame containing Ribo-seq reads.
#' @param ranges Genomic ranges (\code{GRanges} object) for CDS or eORF.
#' @param strand Strand information ("+" or "-").
#' @return Data frame with an added \code{frame} column indicating the reading frame.
assign_frames <- function(Ribo_data, ranges, strand) {
  # Check if ranges are provided
  if (length(ranges) > 0) {
    # Sort exons by strand orientation
    if (strand == "+") {
      exons <- sort(ranges, decreasing = FALSE)
    } else {
      exons <- sort(ranges, decreasing = TRUE)
    }

    # Initialize vectors for genomic and transcript positions
    positions <- integer(0)
    tx_positions <- integer(0)
    cum_len <- 0

    # Map genomic positions to transcript positions
    for (idx in seq_along(exons)) {
      exon <- exons[idx]
      pos <- seq(start(exon), end(exon))
      if (strand == "-") {
        pos <- rev(pos)
      }
      len <- length(pos)
      tx_pos <- seq_len(len) + cum_len
      positions <- c(positions, pos)
      tx_positions <- c(tx_positions, tx_pos)
      cum_len <- cum_len + len
    }

    # Create a data frame for position to transcript position mapping
    position_df <- data.frame(position = positions, tx_pos = tx_positions)

    # Compute frame based on transcript positions
    position_df$frame <- (position_df$tx_pos - 1) %% 3
    position_df$frame <- factor(position_df$frame, levels = c(0,1,2))

    # Merge frame information into Ribo_data
    Ribo_data <- merge(Ribo_data, position_df[, c("position", "frame")], by = "position", all.x = TRUE)
  } else {
    # No ranges, assign NA frames
    Ribo_data$frame <- factor(NA, levels = c(0,1,2))
  }

  return(Ribo_data)
}

# Corrected helper function to assign frames with extension into UTRs
#' @title Assign Frames with Extension into UTRs
#' @description
#' Assigns reading frames to Ribo-seq reads with extension into UTRs. This is useful for analyzing
#' alternative ORFs (including those in UTRs) by extending the coding region.
#'
#' @param Ribo_data Data frame containing Ribo-seq reads.
#' @param cds_ranges CDS ranges as a \code{GRanges} object.
#' @param exons Exon ranges as a \code{GRanges} object.
#' @param fExtend Number of nucleotides to extend into the 5' UTR.
#' @param tExtend Number of nucleotides to extend into the 3' UTR.
#' @param strand Strand information ("+" or "-").
#' @return Data frame with an added \code{frame} column indicating the reading frame.
assign_frames_with_extension <- function(Ribo_data, cds_ranges, exons, fExtend, tExtend, strand) {
  # Check that exons and CDS ranges are provided
  if (length(exons) > 0 && length(cds_ranges) > 0) {
    # Sort exons by strand
    if (strand == "+") {
      exons <- sort(exons, decreasing = FALSE)
    } else {
      exons <- sort(exons, decreasing = TRUE)
    }

    # Build transcript coordinate mapping from exons
    positions <- integer(0)
    tx_positions <- integer(0)
    cum_len <- 0
    for (idx in seq_along(exons)) {
      exon <- exons[idx]
      pos <- seq(start(exon), end(exon))
      if (strand == "-") {
        pos <- rev(pos)
      }
      len <- length(pos)
      tx_pos <- seq_len(len) + cum_len
      positions <- c(positions, pos)
      tx_positions <- c(tx_positions, tx_pos)
      cum_len <- cum_len + len
    }
    position_df <- data.frame(position = positions, tx_pos = tx_positions)

    # Map CDS genomic positions to transcript positions
    cds_positions <- unlist(lapply(seq_along(cds_ranges), function(idx) {
      seq(start(cds_ranges[idx]), end(cds_ranges[idx]))
    }))
    if (strand == "-") {
      cds_positions <- rev(cds_positions)
    }
    cds_tx_pos <- position_df$tx_pos[position_df$position %in% cds_positions]

    # Identify start and end of CDS in transcript coordinates
    cds_start_tx <- min(cds_tx_pos)
    cds_end_tx <- max(cds_tx_pos)

    # Extend the CDS in transcript coordinates by fExtend and tExtend
    extended_start_tx <- max(1, cds_start_tx - fExtend)
    extended_end_tx <- min(max(position_df$tx_pos), cds_end_tx + tExtend)
    extended_tx_pos <- seq(extended_start_tx, extended_end_tx)

    # Map extended transcript positions back to genomic positions
    extended_positions <- position_df$position[position_df$tx_pos %in% extended_tx_pos]

    # Compute frames relative to the original CDS start
    frames <- ((extended_tx_pos - cds_start_tx) %% 3)
    frames <- factor(frames, levels = c(0,1,2))

    # Build a data frame linking positions to frames
    frame_df <- data.frame(tx_pos = extended_tx_pos, position = extended_positions, frame = frames)

    # Remove duplicate rows (if any) so that each position maps to a unique frame
    frame_df <- frame_df[!duplicated(frame_df$position), ]

    # Merge frame information into Ribo_data
    Ribo_data <- merge(Ribo_data, frame_df[, c("position", "frame")], by = "position", all.x = TRUE)
  } else {
    # If no exons or CDS ranges, assign NA frames
    Ribo_data$frame <- factor(NA, levels = c(0,1,2))
  }

  return(Ribo_data)
}


# Function to exclude Ribo-seq reads that overlap eORF regions
# Helper function to exclude eORF reads from main Ribo-seq data
#' @title Exclude eORF Reads from Main Ribo-seq Data
#' @description
#' Excludes Ribo-seq reads overlapping eORF regions from the main Ribo-seq data, ensuring the main ORF analysis is not confounded by eORF reads.
#'
#' @param Ribo_data Data frame containing Ribo-seq reads.
#' @param eORFTxInfo An \code{eORF_info} object containing eORF information.
#' @param strand Strand information ("+" or "-").
#' @return Filtered Ribo-seq data frame without eORF-overlapping reads.
exclude_eORF_reads <- function(Ribo_data, eORFTxInfo, strand) {
  # If no eORF info, just return the original Ribo_data
  if (is.null(eORFTxInfo) || length(eORFTxInfo$xlim.eORF) == 0) {
    return(Ribo_data)
  }

  # Create a GRanges object from Ribo_data
  Ribo_gr <- GRanges(
    seqnames = Ribo_data$chr,
    ranges = IRanges(Ribo_data$position, Ribo_data$position),
    strand = Ribo_data$strand
  )

  # Create a GRangesList of eORF regions
  eORF_grl <- GRangesList(eORFTxInfo$xlim.eORF)

  # Find overlaps between Ribo-seq reads and eORF regions
  overlaps <- findOverlaps(Ribo_gr, unlist(eORF_grl))

  # Remove reads overlapping with eORFs
  if (length(overlaps) > 0) {
    Ribo_data <- Ribo_data[-queryHits(overlaps), ]
  }

  return(Ribo_data)
}



#' Get RNA-seq coverage for a sample
#'
#' Processes RNA-seq input to return a coverage vector over a genomic range,
#' safely handling cases where the requested region extends beyond the data.
#'
#' @param RNAseq_sample List specifying type ("bam", "bigwig", "bedgraph") and data.
#' @param gene_range GRanges object specifying the region of interest.
#' @param strand_info Character, "+" or "-", indicating the gene's strand.
#' @return Numeric vector of coverage values (zeros outside the available range).
get_RNAseq_coverage <- function(RNAseq_sample, gene_range, strand_info) {
  # Determine desired genomic window
  global_start <- start(gene_range)
  global_end   <- end(gene_range)
  chr          <- as.character(seqnames(gene_range))
  window_len   <- global_end - global_start + 1
  coverage_vec <- numeric(window_len)  # initialize all zeros

  if (RNAseq_sample$type == "bam") {
    # Compute coverage from BAM
    param <- ScanBamParam(which=gene_range, what=c("rname","strand","pos","qwidth"))
    if (RNAseq_sample$paired == "paired") {
      readPairs <- readGAlignmentPairs(RNAseq_sample$file, param=param, strandMode=2)
      readPairs <- readPairs[strand(readPairs)==strand_info]
      cvg <- coverage(readPairs)
    } else {
      alignments <- readGAlignments(RNAseq_sample$file, param=param)
      alignments <- alignments[strand(alignments)==strand_info]
      cvg <- coverage(alignments)
    }

    # Safe extraction: only fill positions that actually exist
    if (!is.null(cvg[[chr]])) {
      chr_cvg <- as.numeric(cvg[[chr]])
      chr_len <- length(chr_cvg)
      # compute overlap of requested window with [1, chr_len]
      start_idx <- max(global_start, 1)
      end_idx   <- min(global_end,   chr_len)
      if (start_idx <= end_idx) {
        out_idx   <- (start_idx: end_idx) - global_start + 1
        coverage_vec[out_idx] <- chr_cvg[start_idx: end_idx]
      }
    }
    return(coverage_vec)

  } else if (RNAseq_sample$type %in% c("bigwig", "bedgraph")) {
    # Import coverage from bigwig or bedgraph
    strand_file <- if (strand_info == "+") RNAseq_sample$plus else RNAseq_sample$minus
    cvg_gr <- rtracklayer::import(strand_file, which=gene_range, format=RNAseq_sample$type)
    cvg_gr <- cvg_gr[strand(cvg_gr) == strand_info | strand(cvg_gr) == "*"]
    if (length(cvg_gr) > 0) {
      # Compute coverage
      cvg_full <- coverage(cvg_gr, weight=score(cvg_gr))
      if (!is.null(cvg_full[[chr]])) {
        chr_cvg <- as.numeric(cvg_full[[chr]])
        chr_len <- length(chr_cvg)
        start_idx <- max(global_start, 1)
        end_idx   <- min(global_end,   chr_len)
        if (start_idx <= end_idx) {
          out_idx   <- (start_idx: end_idx) - global_start + 1
          coverage_vec[out_idx] <- chr_cvg[start_idx: end_idx]
        }
      }
    }
    return(coverage_vec)

  } else {
    stop("Invalid type for RNAseq_sample: must be 'bam', 'bigwig', or 'bedgraph'")
  }
}



#' Get Ribo-seq data for a sample
#'
#' Processes Ribo-seq input to return a data frame with position and count.
#'
#' @param Riboseq_sample List specifying type ("tabular", "bigwig", "bedgraph") and data.
#' @param gene_range GRanges object specifying the region of interest.
#' @param strand_info Character, "+" or "-", indicating the gene's strand.
#' @return Data frame with columns "position", "count", "strand", "chr".
get_Riboseq_data <- function(Riboseq_sample, gene_range, strand_info) {
  if (Riboseq_sample$type == "tabular") {
    # Use provided data frame, filter to gene range and strand
    df <- Riboseq_sample$data
    df <- df[df$chr == as.character(seqnames(gene_range)) &
             df$position >= start(gene_range) &
             df$position <= end(gene_range) &
             df$strand == strand_info, ]
    return(df)
  } else if (Riboseq_sample$type %in% c("bigwig", "bedgraph")) {
    # Import from bigwig or bedgraph
    strand_file <- if (strand_info == "+") Riboseq_sample$plus else Riboseq_sample$minus
    gr <- rtracklayer::import(strand_file, which=gene_range, format=Riboseq_sample$type)
    if (length(gr) == 0) {
      return(data.frame(position=integer(0), count=numeric(0), strand=character(0), chr=character(0)))
    }
    # Expand ranges to individual positions (1bp resolution for Ribo-seq)
    positions <- unlist(lapply(seq_along(gr), function(i) {
      seq(start(gr[i]), end(gr[i]))
    }))
    counts <- rep(score(gr), width(gr))
    df <- data.frame(
      position = positions,
      count = counts,
      strand = strand_info,
      chr = as.character(seqnames(gene_range))
    )
    # Aggregate counts if multiple entries per position
    df <- aggregate(count ~ position + strand + chr, data=df, sum)
    return(df)
  } else {
    stop("Invalid type for Riboseq_sample: must be 'tabular', 'bigwig', or 'bedgraph'")
  }
}


#' Create input lists for ggRNA and ggRibo
#'
#' A wrapper to construct RNAseq and Riboseq input lists from file paths with automatic file format detection.
#'
#' @param rna_files List of RNA-seq file paths (BAM) or named list with plus/minus (bigwig/bedgraph).
#' @param ribo_files List of Ribo-seq file paths (tabular) or named list with plus/minus (bigwig/bedgraph).
#' @param sample_names Character vector of sample names.
#' @param rna_types Character vector of file types for RNA-seq ("bam", "bigwig", "bedgraph"). If NULL, detected automatically.
#' @param ribo_types Character vector of file types for Ribo-seq ("tabular", "bigwig", "bedgraph"). If NULL, detected automatically.
#' @param include_rna Logical, whether to include RNA-seq data.
#' @param rna_paired Character vector, "paired" or "single" for BAM files.
#' @return List with RNAseq and Riboseq inputs for ggRNA/ggRibo.
create_seq_input <- function(rna_files = NULL, ribo_files = NULL, sample_names,
                            rna_types = NULL, ribo_types = NULL, include_rna = TRUE,
                            rna_paired = rep("paired", length(sample_names))) {
  if (length(sample_names) != length(rna_files) && include_rna) {
    stop("Number of RNA-seq files must match sample_names when include_rna is TRUE")
  }
  if (length(sample_names) != length(ribo_files) && !is.null(ribo_files)) {
    stop("Number of Ribo-seq files must match sample_names")
  }

  RNAseq <- list()
  Riboseq <- list()

  # Detect file type based on extension
  detect_file_type <- function(file) {
    if (is.list(file) && all(c("plus", "minus") %in% names(file))) {
      ext <- tolower(tools::file_ext(file$plus))
      if (ext %in% c("bigwig", "bw")) return("bigwig")
      if (ext %in% c("bedgraph", "bedGraph", "BedGraph", "bdg")) return("bedgraph")
    } else {
      ext <- tolower(tools::file_ext(file))
      if (ext == "bam") return("bam")
      if (ext %in% c("bigwig", "bw")) return("bigwig")
      if (ext %in% c("bedgraph", "bedGraph", "BedGraph", "bdg")) return("bedgraph")
    }
    return("tabular") # Default for ribo_files if not bam/bigwig/bedgraph
  }

  # Process RNA-seq inputs
  if (include_rna && !is.null(rna_files)) {
    if (is.null(rna_types)) {
      rna_types <- sapply(rna_files, detect_file_type)
    }
    if (length(rna_types) != length(rna_files)) stop("rna_types must match rna_files length")
    if (any(!rna_types %in% c("bam", "bigwig", "bedgraph"))) {
      stop("rna_types must be 'bam', 'bigwig', or 'bedgraph'")
    }
    for (i in seq_along(rna_files)) {
      if (rna_types[i] == "bam") {
        if (is.null(rna_paired)) stop("rna_paired must be provided for BAM files")
        RNAseq[[i]] <- list(type = "bam", file = rna_files[[i]], paired = rna_paired[i])
      } else {
        RNAseq[[i]] <- list(type = rna_types[i], plus = rna_files[[i]]$plus, minus = rna_files[[i]]$minus)
      }
    }
  }

  # Process Ribo-seq inputs
  if (!is.null(ribo_files)) {
    if (is.null(ribo_types)) {
      ribo_types <- sapply(ribo_files, detect_file_type)
    }
    if (length(ribo_types) != length(ribo_files)) stop("ribo_types must match ribo_files length")
    if (any(!ribo_types %in% c("tabular", "bigwig", "bedgraph"))) {
      stop("ribo_types must be 'tabular', 'bigwig', or 'bedgraph'")
    }
    for (i in seq_along(ribo_files)) {
      if (ribo_types[i] == "tabular") {
        df <- read.delim(ribo_files[[i]], header = FALSE, stringsAsFactors = FALSE, sep = "\t")
        colnames(df) <- c("count", "chr", "position", "strand")
        Riboseq[[i]] <- list(type = "tabular", data = df)
      } else {
        Riboseq[[i]] <- list(type = ribo_types[i], plus = ribo_files[[i]]$plus, minus = ribo_files[[i]]$minus)
      }
    }
  }
  RNAseqBamPairorSingle= rna_paired
  assign("RNAseqBamPairorSingle", RNAseqBamPairorSingle, envir = .GlobalEnv)
  assign("Samples", sample_names, envir = .GlobalEnv)

  return(list(RNAseq = if (include_rna) RNAseq else NULL, Riboseq = Riboseq))
}

#' Plot Gene Transcript Model
#'
#' This function creates a gene model plot showing exons, UTRs, CDS, and optional eORFs for a given gene and its isoforms.
#'
#' @param GeneTxInfo A `Gene_info` object containing gene and transcript information.
#' @param eORFTxInfo An optional `eORF_info` object containing eORF information.
#' @param XLIM Not used in the current function (legacy parameter).
#' @param plot_ORF_ranges Logical, whether to plot ORF ranges.
#' @param plot_range An optional numeric vector specifying the genomic range to plot (start and end positions).
#' @param transcript_label_font_size Optional numeric value to control the font size of the transcript ID labels.
#'
#' @return A `ggplot` object representing the gene model.
#' @export
plotGeneTxModel <- function(GeneTxInfo = GeneTxInfo, eORFTxInfo = NULL, XLIM = NULL, plot_ORF_ranges = TRUE, plot_range = NULL, transcript_label_font_size = 10) {
  # Load necessary libraries
  # Extract information from the GeneTxInfo object
  isoforms <- GeneTxInfo$num_isoforms
  genelim <- c(GeneTxInfo$range_left, GeneTxInfo$range_right)
  tx_names <- GeneTxInfo$tx_names
  tx_id <- GeneTxInfo$tx_id
  strand <- GeneTxInfo$strand

  # Prepare lists to store plotting data
  plot_data_list <- list()
  line_data_list <- list()
  idx <- 1

  # Sort transcripts so the main transcript (tx_id) is plotted first at the top,
  # and other transcripts are sorted alphabetically below
  other_tx_names <- setdiff(tx_names, tx_id)
  sorted_tx_names <- c(tx_id, sort(other_tx_names))

  # Assign y-axis positions for each isoform, main transcript at the top
  y_step <- 0.3
  y_positions <- seq(1, by = y_step, length.out = length(sorted_tx_names))
  isoform_positions <- data.frame(
    isoform = sorted_tx_names,
    y = rev(y_positions),
    stringsAsFactors = FALSE
  )

  # Create a mapping from isoform to its y-axis position
  isoform_y_map <- setNames(isoform_positions$y, isoform_positions$isoform)

  # Loop through each isoform to generate plotting data
  for (isoform in isoform_positions$isoform) {
    y_value <- isoform_y_map[isoform]
    isoform_data_list <- list()
    isoform_idx <- 1

    # Get exon ranges for this isoform
    exons_gr <- GeneTxInfo$exonByYFGtx[[isoform]]
    if (length(exons_gr) == 0) {
      # If no exons found, print a warning and continue to next isoform
      warning(paste("Exons for isoform", isoform, "not found in exonByYFGtx"))
      next
    }

    # Keep original exon ranges
    exons_gr_original <- exons_gr

    # Extract original CDS, fiveUTR, and threeUTR before truncation
    original_cds_ranges <- GeneTxInfo$xlimCds[[isoform]]
    original_fiveUTR_gr <- NULL
    if (isoform %in% names(GeneTxInfo$fiveUTRByYFGtx)) {
      original_fiveUTR_gr <- unlist(GeneTxInfo$fiveUTRByYFGtx[isoform])
    }
    original_threeUTR_gr <- NULL
    if (isoform %in% names(GeneTxInfo$threeUTRByYFGtx)) {
      original_threeUTR_gr <- unlist(GeneTxInfo$threeUTRByYFGtx[isoform])
    }

    # If a custom plot_range is specified, intersect exons with this range to truncate
    segment_gr <- NULL
    if (!is.null(plot_range)) {
      segment_gr <- GRanges(seqnames = GeneTxInfo$chr,
                            ranges = IRanges(plot_range[1], plot_range[2]),
                            strand = GeneTxInfo$strand)
      exons_gr_truncated <- pintersect(exons_gr, segment_gr)
      exons_gr_truncated <- exons_gr_truncated[width(exons_gr_truncated) > 0]
      if (length(exons_gr_truncated) == 0) {
        # If no exons remain after truncation, skip this isoform
        next
      }
      exons_gr <- exons_gr_truncated
    }

    # Determine the transcript start and end after possible truncation
    transcript_start <- min(start(exons_gr))
    transcript_end <- max(end(exons_gr))

    # Helper function to truncate a given GRanges feature to the plotting segment
    truncate_feature <- function(feature_gr, orig_feature_gr) {
      if (length(feature_gr) == 0) return(NULL)

      # If no truncation range is given, just return original starts/ends
      if (is.null(segment_gr)) {
        df <- data.frame(
          start = start(feature_gr),
          end = end(feature_gr),
          orig_start = start(feature_gr),
          orig_end = end(feature_gr),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        return(df)
      } else {
        # If truncation range is provided, intersect and find truncated segments
        truncated_gr <- pintersect(feature_gr, segment_gr)
        truncated_gr <- truncated_gr[width(truncated_gr) > 0]
        if (length(truncated_gr) == 0) return(NULL)

        # For each truncated range, find original boundaries
        out_list <- list()
        for (i in seq_along(truncated_gr)) {
          tgr <- truncated_gr[i]
          hit <- findOverlaps(tgr, feature_gr)
          if (length(hit) > 0) {
            f_idx <- subjectHits(hit)[1]
            orig_start_val <- start(feature_gr[f_idx])
            orig_end_val <- end(feature_gr[f_idx])
          } else {
            orig_start_val <- start(tgr)
            orig_end_val <- end(tgr)
          }
          out_list[[i]] <- data.frame(
            start = start(tgr),
            end = end(tgr),
            orig_start = orig_start_val,
            orig_end = orig_end_val,
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }
        final_df <- do.call(rbind, out_list)
        return(final_df)
      }
    }

    # Truncate CDS, fiveUTR, and threeUTR features if necessary
    cds_df_raw <- NULL
    if (!is.null(original_cds_ranges) && length(original_cds_ranges) > 0) {
      cds_df_raw <- truncate_feature(original_cds_ranges, original_cds_ranges)
      if (!is.null(cds_df_raw)) {
        cds_df <- data.frame(
          start = cds_df_raw$start,
          end = cds_df_raw$end,
          y = y_value,
          feature = "CDS",
          isoform = isoform,
          height_factor = 1,
          orf_id = NA,
          orig_start = cds_df_raw$orig_start,
          orig_end = cds_df_raw$orig_end,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        isoform_data_list[[isoform_idx]] <- cds_df
        isoform_idx <- isoform_idx + 1
      }
    }

    fiveUTR_df_raw <- NULL
    if (!is.null(original_fiveUTR_gr) && length(original_fiveUTR_gr) > 0) {
      fiveUTR_df_raw <- truncate_feature(original_fiveUTR_gr, original_fiveUTR_gr)
      if (!is.null(fiveUTR_df_raw)) {
        fiveUTR_df <- data.frame(
          start = fiveUTR_df_raw$start,
          end = fiveUTR_df_raw$end,
          y = y_value,
          feature = "5' UTR",
          isoform = isoform,
          height_factor = 1,
          orf_id = NA,
          orig_start = fiveUTR_df_raw$orig_start,
          orig_end = fiveUTR_df_raw$orig_end,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        isoform_data_list[[isoform_idx]] <- fiveUTR_df
        isoform_idx <- isoform_idx + 1
      }
    }

    threeUTR_df_raw <- NULL
    if (!is.null(original_threeUTR_gr) && length(original_threeUTR_gr) > 0) {
      threeUTR_df_raw <- truncate_feature(original_threeUTR_gr, original_threeUTR_gr)
      if (!is.null(threeUTR_df_raw)) {
        threeUTR_df <- data.frame(
          start = threeUTR_df_raw$start,
          end = threeUTR_df_raw$end,
          y = y_value,
          feature = "3' UTR",
          isoform = isoform,
          height_factor = 1,
          orf_id = NA,
          orig_start = threeUTR_df_raw$orig_start,
          orig_end = threeUTR_df_raw$orig_end,
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        isoform_data_list[[isoform_idx]] <- threeUTR_df
        isoform_idx <- isoform_idx + 1
      }
    }

    # If we have no features, treat as ncRNA (moved to before eORF features)
    if (length(isoform_data_list) == 0) {
      exons_raw <- truncate_feature(exons_gr_original, exons_gr_original)
      if (!is.null(exons_raw)) {
        ncRNA_df <- data.frame(
          start=exons_raw$start,
          end=exons_raw$end,
          y=y_value,
          feature="ncRNA",
          isoform=isoform,
          height_factor=1,
          orf_id=NA,
          orig_start = exons_raw$orig_start,
          orig_end = exons_raw$orig_end,
          stringsAsFactors=FALSE,
          row.names = NULL
        )
        isoform_data_list[[isoform_idx]] <- ncRNA_df
        isoform_idx <- isoform_idx+1
      }
    }
    
    # If eORF info is provided and plot_ORF_ranges is TRUE, include eORF features
    if (!is.null(eORFTxInfo)) {
      for (eORF_idx in seq_along(eORFTxInfo$eORF.tx_id)) {
        eORF_ranges <- eORFTxInfo$xlim.eORF[[eORF_idx]]

        ### Skip if the eORF doesn't overlap this isoform's exons
        overlap_exons <- findOverlaps(eORF_ranges, exons_gr_original, type="within")
        if (length(unique(queryHits(overlap_exons))) < length(eORF_ranges)) {
          next
        }

        if (length(eORF_ranges) > 0) {
          original_eORF_ranges <- eORF_ranges

          # Check overlap with CDS, 5'UTR, and 3'UTR to classify eORF type
          overlaps_CDS <- FALSE
          overlaps_fiveUTR <- FALSE
          overlaps_threeUTR <- FALSE

          if (!is.null(original_cds_ranges) && length(original_cds_ranges) > 0) {
            overlap_cds <- findOverlaps(original_eORF_ranges, original_cds_ranges)
            if (length(overlap_cds) > 0) {
              overlaps_CDS <- TRUE
            }
          }

          if (!is.null(original_fiveUTR_gr) && length(original_fiveUTR_gr) > 0) {
            overlap_five <- findOverlaps(original_eORF_ranges, original_fiveUTR_gr)
            if (length(overlap_five) > 0) {
              overlaps_fiveUTR <- TRUE
            }
          }

          if (!is.null(original_threeUTR_gr) && length(original_threeUTR_gr) > 0) {
            overlap_three <- findOverlaps(original_eORF_ranges, original_threeUTR_gr)
            if (length(overlap_three) > 0) {
              overlaps_threeUTR <- TRUE
            }
          }

          # Assign feature label based on overlaps
          if (overlaps_fiveUTR) {
            if (overlaps_CDS) {
              feature_label <- "ouORF"
            } else {
              feature_label <- "uORF"
            }
          } else if (overlaps_threeUTR) {
            if (overlaps_CDS) {
              feature_label <- "odORF"
            } else {
              feature_label <- "dORF"
            }
          } else if (overlaps_CDS) {
            feature_label <- "nORF"
          } else {
            feature_label <- "ORF"
          }

          # Truncate the eORF if necessary
          eORF_df_raw <- truncate_feature(original_eORF_ranges, original_eORF_ranges)
          if (!is.null(eORF_df_raw)) {
            eORF_df <- data.frame(
              start = eORF_df_raw$start,
              end = eORF_df_raw$end,
              y = y_value,
              feature = feature_label,
              isoform = isoform,
              orf_id = feature_label,
              height_factor = ifelse(overlaps_CDS, 2/3, 1),
              orig_start = eORF_df_raw$orig_start,
              orig_end = eORF_df_raw$orig_end,
              stringsAsFactors = FALSE,
              row.names = NULL
            )
            isoform_data_list[[isoform_idx]] <- eORF_df
            isoform_idx <- isoform_idx + 1
          }
        }
      }
    }

    # Combine isoform data into a single data frame
    if (length(isoform_data_list) > 0) {
      isoform_df <- do.call(rbind, isoform_data_list)
      plot_data_list[[idx]] <- isoform_df
      idx <- idx + 1
    }

    # Add intron lines if there is more than one exon
    intron_df_list <- list()
    if (length(exons_gr) > 1) {
      exons_sorted <- exons_gr[order(start(exons_gr))]
      for (i in seq_len(length(exons_sorted) - 1)) {
        intron_start <- end(exons_sorted[i])
        intron_end <- start(exons_sorted[i + 1])
        intron_df <- data.frame(
          xstart = intron_start,
          xend = intron_end,
          y = y_value,
          isoform = isoform,
          row.names = NULL
        )
        intron_df_list[[length(intron_df_list) + 1]] <- intron_df
      }
    }

    # If truncated by plot_range, check for flanking introns beyond segment
    if (!is.null(segment_gr)) {
      segment_left <- start(segment_gr)
      segment_right <- end(segment_gr)
      exons_sorted <- exons_gr[order(start(exons_gr))]

      if (start(exons_sorted[1]) > segment_left) {
        intron_df <- data.frame(
          xstart = segment_left,
          xend = start(exons_sorted[1]),
          y = y_value,
          isoform = isoform,
          row.names = NULL
        )
        intron_df_list[[length(intron_df_list) + 1]] <- intron_df
      }

      if (end(exons_sorted[length(exons_sorted)]) < segment_right) {
        intron_df <- data.frame(
          xstart = end(exons_sorted[length(exons_sorted)]),
          xend = segment_right,
          y = y_value,
          isoform = isoform,
          row.names = NULL
        )
        intron_df_list[[length(intron_df_list) + 1]] <- intron_df
      }
    }

    # Combine intron data into a single data frame if any
    if (length(intron_df_list) > 0) {
      intron_data <- do.call(rbind, intron_df_list)
      line_data_list[[length(line_data_list) + 1]] <- intron_data
    }
  }

  # Combine all isoforms' data
  if (length(plot_data_list) > 0) {
    plot_data <- do.call(rbind, plot_data_list)
  } else {
    stop("No valid exons or features found to plot.")
  }

  # Combine intron data if available
  if (length(line_data_list) > 0) {
    line_data <- do.call(rbind, line_data_list)
  } else {
    line_data <- data.frame()
  }

  # Ensure orf_id column exists
  if (!"orf_id" %in% names(plot_data)) {
    plot_data$orf_id <- NA
  }
  plot_data$orf_id <- as.character(plot_data$orf_id)

  # Define the order of features for plotting
  feature_order <- c("uORF", "ouORF", "nORF", "ORF", "odORF", "dORF", "5' UTR", "CDS", "3' UTR", "ncRNA")
  plot_data$feature <- factor(plot_data$feature, levels = feature_order)

  # Determine rectangle height
  plot_data$height <- 0.08 * plot_data$height_factor
  plot_data$ymin <- plot_data$y - plot_data$height
  plot_data$ymax <- plot_data$y + plot_data$height

  # Identify features that are present
  unique_features <- levels(plot_data$feature)[levels(plot_data$feature) %in% plot_data$feature]

  # Assign colors to different feature types
  feature_colors <- c(
    "uORF" = "yellow",
    "ouORF" = "#FFD700",
    "nORF" = "orange",
    "ORF" = "lightblue",
    "odORF" = "#FFD700",
    "dORF" = "yellow",
    "5' UTR" = "lightgrey",
    "CDS" = "black",
    "3' UTR" = "white",
    "ncRNA" = "#FFB6C1"
  )

  feature_colors <- feature_colors[unique_features]

  # Adjust legend settings based on number of feature types
  num_legend_items <- length(unique_features)
  legend_text_size <- 8
  base_key_size <- 1
  if (num_legend_items > 3) {
    key_size <- base_key_size * 3 / num_legend_items
  } else {
    key_size <- base_key_size
  }
  key_size <- max(0.8, key_size)

  p_gene <- ggplot()

  # Plot intron segments as horizontal lines between exons
  if (nrow(line_data) > 0) {
    p_gene <- p_gene +
      geom_segment(data = line_data, aes(x = xstart, xend = xend, y = y, yend = y),
                   color = "black", inherit.aes = FALSE)
  }

  # Draw filled rectangles for features (exons, UTRs, ORFs) without borders
  p_gene <- p_gene +
    geom_rect(data = plot_data,
              aes(xmin = start - 0.5, xmax = end + 0.5, ymin = ymin, ymax = ymax, fill = feature),
              color = NA, inherit.aes = FALSE)

  # Draw top and bottom borders for features
  p_gene <- p_gene +
    geom_segment(data = plot_data,
                 aes(x = start - 0.5, xend = end + 0.5, y = ymax, yend = ymax),
                 color = "black") +
    geom_segment(data = plot_data,
                 aes(x = start - 0.5, xend = end + 0.5, y = ymin, yend = ymin),
                 color = "black")

  # Draw vertical borders at start and end of features only if not truncated
  left_borders <- plot_data[plot_data$start == plot_data$orig_start, ]
  if (nrow(left_borders) > 0) {
    p_gene <- p_gene +
      geom_segment(data = left_borders,
                   aes(x = start - 0.5, xend = start - 0.5, y = ymin, yend = ymax),
                   color = "black")
  }

  right_borders <- plot_data[plot_data$end == plot_data$orig_end, ]
  if (nrow(right_borders) > 0) {
    p_gene <- p_gene +
      geom_segment(data = right_borders,
                   aes(x = end + 0.5, xend = end + 0.5, y = ymin, yend = ymax),
                   color = "black")
  }

  # Define fill scale for feature colors and set up the legend
  p_gene <- p_gene +
    scale_fill_manual(
      name = "Feature",
      values = feature_colors,
      breaks = unique_features,
      labels = unique_features,
      guide = guide_legend(
        override.aes = list(
          fill = feature_colors,
          color = "black"  # Border color for legend keys
        ),
        ncol = 1,
        keyheight = unit(key_size, "lines"),
        keywidth = unit(1, "lines")
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),
      legend.key.size = unit(key_size, "lines"),
      axis.text.y = element_text(size = transcript_label_font_size),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(size = 12),
      plot.margin = unit(c(-0.2, 0.2, 0, 1.5), "lines"),
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) +
    xlab("") + #Genomic Position
    ylab("") +
    coord_cartesian(clip = "off")

  # Reverse x-axis if gene is on the negative strand
  if (strand == "-") {
    p_gene <- p_gene + scale_x_reverse(limits = c(max(genelim) + 0.5, min(genelim) - 0.5))
  } else {
    p_gene <- p_gene + scale_x_continuous(limits = c(min(genelim) - 0.5, max(genelim) + 0.5))
  }

  # Format isoform labels. The main transcript is bold, others are plain
  labels <- sapply(isoform_positions$isoform, function(x) {
    if (x == tx_id) {
      paste0("bold('", x, "')")
    } else {
      paste0("'", x, "'")
    }
  })
  labels <- parse(text = labels)

  padding <- 0.1
  p_gene <- p_gene + scale_y_continuous(
    breaks = isoform_positions$y,
    labels = labels,
    limits = c(min(isoform_positions$y) - padding,
               max(isoform_positions$y) + padding)
  )

  return(p_gene)
}

#' Plot DNA and Amino Acid Sequences
#'
#' Generates a plot showing DNA nucleotides (if the region is \(\le\)201 bases)
#' and amino acids (always as colored tiles). The AA letters themselves are only shown if
#' \(\le\)201 bases. If over 201 bases, we still color the AA tiles (e.g., start/stop codons),
#' but omit the text labels to avoid overprinting.
#'
#' @param GeneTxInfo A \code{Gene_info} object containing gene-specific information, including exons.
#' @param plot_range Optional numeric vector of length 2 specifying the transcript range (in transcript coordinates).
#'   If provided, we show only that region (exon-based). Otherwise, we show the full transcript sequence.
#' @param FASTA A \code{BSgenome} object containing the reference genome sequences.
#'
#' @return A \code{ggplot2} object representing the DNA (tiles + letters if \(\le\)201nt) and amino acids (colored tiles always, letters if \(\le\)201nt).
#' @param nucleotide_color_scheme A character string specifying which color scheme
#'   to use for the nucleotides. Defaults to \code{"default"}. If set to
#'   \code{"colorblind"}, a color‐blind friendly palette is used (A=green, T=vermillion,
#'   C=blue, G=yellow). Non‐canonical bases (N) remain grey.
#' @export

plotDNAandAA <- function(GeneTxInfo, plot_range = NULL, FASTA = NULL, nucleotide_color_scheme = "default") {
  # Check if FASTA is provided
  if (is.null(FASTA)) {
    stop("FASTA must be provided as a BSgenome object.")
  }
  
  # Extract the gene limits from GeneTxInfo
  genelim <- c(GeneTxInfo$range_left, GeneTxInfo$range_right)
  strand <- GeneTxInfo$strand
  chr <- GeneTxInfo$chr
  
  # Check for NA in gene limits
  if (any(is.na(genelim))) {
    stop(paste("Gene range is undefined for", GeneTxInfo$tx_id, "- check transcript data in GRangeInfo"))
  }
  
  # Adjust genomic limits based on chromosome length
  chrom_length <- seqlengths(FASTA)[chr]
  if (is.null(chrom_length) || is.na(chrom_length)) {
    stop(paste("Chromosome", chr, "not found in FASTA"))
  }
  genelim_adj <- pmax(pmin(genelim, chrom_length), 1)
  
  # If a custom plot_range is given, use that instead
  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)
    genelim_adj <- pmax(pmin(plot_range, chrom_length), 1)
  }
  
  # Determine range length and decide whether to suppress labels
  range_length <- abs(diff(range(genelim_adj))) + 1
  suppress_labels <- range_length > 201
  long_range <- range_length > 201
  
  # Extract DNA sequence from FASTA within the adjusted range
  seq_region <- GRanges(
    seqnames = chr,
    ranges = IRanges(min(genelim_adj), max(genelim_adj)),
    strand = "+"
  )
  seqs <- getSeq(FASTA, seq_region)
  dna_seq <- as.character(seqs)
  
  # Reverse complement if on the negative strand
  if (strand == "-") {
    dna_seq <- as.character(reverseComplement(DNAString(dna_seq)))
  }
  
  # Split DNA sequence into individual nucleotides
  dna_chars <- unlist(strsplit(dna_seq, split = ""))
  
  # Assign positions along the sequence
  if (strand == "+") {
    positions_seq <- seq(min(genelim_adj), max(genelim_adj))
  } else {
    positions_seq <- seq(max(genelim_adj), min(genelim_adj), by = -1)
  }
  
  # Prepare a data frame for nucleotides
  dna_df <- data.frame(position = positions_seq,
                       nucleotide = dna_chars,
                       stringsAsFactors = FALSE,
                       row.names = NULL)
  
  # Define colors for nucleotides based on the selected color scheme
  if (tolower(nucleotide_color_scheme) == "colorblind") {
    nucleotide_colors <- c(
      "A" = "#009E73",   # Colorblind friendly Green
      "T" = "#D55E00",   # Colorblind friendly vermilion
      "C" = "#0072B2",   # Colorblind friendly Blue
      "G" = "#F0E442",   # Colorblind friendly yellow
      "N" = "grey"
    )

  } else {
    nucleotide_colors <- c(
      "A" = "#00FF00",   # Green
      "T" = "#FF0200",   # Red
      "C" = "#4747FF",   # Blue
      "G" = "#FFA503",   # Orange
      "N" = "grey"
    )
  }
  dna_df$fill_value <- dna_df$nucleotide
  
  # Adjust font size based on range
  num_nucleotides <- length(dna_chars)
  plot_width <- abs(diff(range(genelim_adj)))
  font_size <- (plot_width / num_nucleotides) * 1.1
  font_size <- max(min(font_size, 5), 2)
  font_size <- font_size * 1.3
  
  # Generate three-frame translations
  frames <- c(0, 1, 2)
  frame_colors <- c("Annotated" = "#F1F1F1", "+1" = "#E6E6E6", "+2" = "#C9C9C9")
  
  # Determine annotated frame based on CDS start (if available)
  cds_ranges <- GeneTxInfo$xlimCds[[GeneTxInfo$tx_id]]
  if (length(cds_ranges) > 0) {
    if (strand == "+") {
      cds_start <- min(start(cds_ranges), na.rm = TRUE)
      if (!is.finite(cds_start)) {
        annotated_frame <- 1
      } else {
        annotated_frame <- ((cds_start - min(genelim_adj)) %% 3) + 1
      }
    } else {
      cds_start <- max(end(cds_ranges), na.rm = TRUE)
      if (!is.finite(cds_start)) {
        annotated_frame <- 1
      } else {
        annotated_frame <- ((max(genelim_adj) - cds_start) %% 3) + 1
      }
    }
  } else {
    annotated_frame <- 1
  }
  annotated_frame_zero_based <- as.integer(annotated_frame - 1)
  if (is.na(annotated_frame_zero_based) || annotated_frame_zero_based < 0 || annotated_frame_zero_based > 2) {
    annotated_frame_zero_based <- 0
  }
  
  # Reorder frames: annotated first, then +1, then +2
  other_frames <- frames[frames != annotated_frame_zero_based]
  frame_order <- c(annotated_frame_zero_based, sort(other_frames))
  
  # Adjust y positions
  if (long_range) {
    frame_y_positions <- c(0.8, 0.6, 0.4)
  } else {
    frame_y_positions <- c(0.6, 0.4, 0.2)
  }
  frame_y_map <- setNames(frame_y_positions, frame_order)
  
  frame_labels <- c("Annotated", "+1", "+2")
  frame_label_map <- setNames(frame_labels, frame_order)
  frame_color_map <- setNames(frame_colors, frame_labels)
  
  aa_sequences <- list()
  
  # Translate each frame
  for (frame in frames) {
    codon_starts <- seq(frame + 1, nchar(dna_seq) - 2, by = 3)
    codon_middles <- codon_starts + 1
    if (length(codon_starts) == 0) next
    
    dna_subseq <- substr(dna_seq, codon_starts[1], codon_starts[length(codon_starts)] + 2)
    dna_string <- DNAString(dna_subseq)
    aa_seq <- suppressWarnings(as.character(translate(dna_string, genetic.code = GENETIC_CODE,no.init.codon = TRUE, if.fuzzy.codon = "X")))
    aa_chars <- unlist(strsplit(aa_seq, split = ""))
    
    aa_positions <- positions_seq[codon_middles]
    aa_df <- data.frame(
      position = aa_positions,
      amino_acid = aa_chars,
      y = frame_y_map[as.character(frame)],
      frame = frame,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    aa_sequences[[frame + 1]] <- aa_df
  }
  
  aa_df_combined <- do.call(rbind, aa_sequences)
  aa_font_size <- font_size * 1.3
  aa_df_combined$frame_label <- frame_label_map[as.character(aa_df_combined$frame)]
  aa_df_combined$fill_value <- aa_df_combined$frame_label
  aa_df_combined$fill_value[aa_df_combined$amino_acid == 'M'] <- 'Start'
  aa_df_combined$fill_value[aa_df_combined$amino_acid == '*'] <- 'Stop'
  
  aa_label_df <- aa_df_combined[aa_df_combined$amino_acid %in% c('M', '*'), ]
  
  fill_colors <- c(
    nucleotide_colors,
    "Start" = "green",
    "Stop" = "red",
    frame_color_map
  )
  
  available_breaks <- intersect(unique(aa_df_combined$fill_value), c("Start", "Stop"))
  aa_start_stop_df <- aa_df_combined[aa_df_combined$fill_value %in% c("Start", "Stop"), ]
  aa_regular_df <- aa_df_combined[!aa_df_combined$fill_value %in% c("Start", "Stop"), ]
  
  p_dna_aa <- ggplot()
  
  # Plot DNA sequence tiles if range is not too long
  if (!long_range) {
    p_dna_aa <- p_dna_aa +
      suppressWarnings(
        geom_tile(
          data = dna_df,
          aes(x = position, y = 0.8, fill = fill_value),
          width = 1,
          height = 0.2,
          color = "darkgrey",
          linewidth = 0.2,
          show.legend = FALSE,
          na.rm = TRUE
        )
      )
    if (!suppress_labels) {
      p_dna_aa <- p_dna_aa +
        suppressWarnings(
          geom_text(
            data = dna_df,
            aes(x = position, y = 0.8, label = nucleotide),
            size = font_size,
            fontface = "plain",
            color = "black",
            show.legend = FALSE,
            na.rm = TRUE
          )
        )
    }
  }
  
  # Plot amino acid tiles
  p_dna_aa <- p_dna_aa +
    suppressWarnings(
      geom_tile(
        data = aa_regular_df,
        aes(x = position, y = y, fill = fill_value),
        width = 3,
        height = 0.2,
        color = ifelse(long_range, NA, "darkgrey"),
        linewidth = 0.2,
        show.legend = TRUE,
        na.rm = TRUE
      )
    ) +
    suppressWarnings(
      geom_tile(
        data = aa_start_stop_df,
        aes(x = position, y = y, fill = fill_value),
        width = 3,
        height = 0.2,
        color = ifelse(long_range, NA, "darkgrey"),
        linewidth = 0.5,
        show.legend = TRUE,
        na.rm = TRUE
      )
    )
  
  if (!suppress_labels) {
    p_dna_aa <- p_dna_aa +
      suppressWarnings(
        geom_text(
          data = aa_df_combined,
          aes(x = position, y = y, label = amino_acid),
          size = aa_font_size,
          fontface = "plain",
          colour = "black",
          show.legend = FALSE,
          na.rm = TRUE
        )
      )
  }
  
  p_dna_aa <- p_dna_aa +
    scale_fill_manual(
      name = NULL,
      values = fill_colors,
      breaks = available_breaks,
      na.value = "grey",
      guide = guide_legend(
        override.aes = list(colour = NA)
      )
    )
  
  # Adjust the x-axis based on strand
  if (strand == "-") {
    p_dna_aa <- p_dna_aa + scale_x_reverse(limits = c(max(genelim_adj), min(genelim_adj)))
  } else {
    p_dna_aa <- p_dna_aa + scale_x_continuous(limits = c(min(genelim_adj), max(genelim_adj)))
  }
  
  # Adjust y-axis limits
  if (long_range) {
    y_limits <- c(0.2, 1)
  } else {
    y_limits <- c(0, 1)
  }
  p_dna_aa <- p_dna_aa +
    scale_y_continuous(limits = y_limits) +
    theme_void() +
    theme(
      plot.margin = unit(c(-1, 0.2, -0.5, 0.2), "lines")
    )
  
  return(p_dna_aa)
}

#' Helper function to resolve gene_id and tx_id based on provided inputs
#' @param GRangeInfo A genomic range information object (e.g., `Txome_Range`) containing annotations.
#' @param gene_id Character string specifying the gene ID of interest.
#' @param tx_id Character string specifying the transcript ID to be used as the main isoform.
#'
#' @return gene_id or tx_id

get_gene_tx <- function(gene_id = NULL, tx_id = NULL, GRangeInfo) {
  if (is.null(gene_id) && is.null(tx_id)) {
    stop("Either gene_id or tx_id must be provided.")
  }
  
  if (!is.null(gene_id) && !is.null(tx_id)) {
    txByYFG <- GRangeInfo$txByGene[gene_id]
    if (length(txByYFG) == 0 || !tx_id %in% txByYFG[[1]]$tx_name) {
      stop(paste("Transcript", tx_id, "is not associated with gene", gene_id))
    }
    return(list(gene_id = gene_id, tx_id = tx_id))
  }
  
  if (is.null(gene_id)) {
    tx_to_gene <- GRangeInfo$tx_to_gene
    matches <- tx_to_gene[tx_to_gene$tx_id == tx_id, ]
    if (nrow(matches) == 0) {
      stop(paste("Transcript", tx_id, "not found in any gene."))
    }
    if (nrow(matches) > 1) {
      warning(paste("Transcript", tx_id, "found in multiple genes:", paste(matches$gene_id, collapse = ", "), ". Using the first one."))
    }
    gene_id <- matches$gene_id[1]
    return(list(gene_id = gene_id, tx_id = tx_id))
  }
  
  if (is.null(tx_id)) {
    tx_to_gene <- GRangeInfo$tx_to_gene
    matches <- tx_to_gene[tx_to_gene$gene_id == gene_id, ]
    if (nrow(matches) == 0) {
      stop(paste("No transcripts found for gene", gene_id))
    }
    # Select the smallest tx_id (already sorted in tx_to_gene)
    tx_id <- matches$tx_id[1]
    return(list(gene_id = gene_id, tx_id = tx_id))
  }
}

#' Plot RNA-seq coverage for a gene
#'
#' The `ggRNA` function creates a plot displaying RNA-seq coverage for a specified gene and transcript.
#' It optionally includes genomic sequences (DNA and amino acids) and gene models (exons, UTRs, CDS).
#' This function is simpler than `ggRibo` and `ggRibo_decom` as it focuses solely on RNA-seq coverage.
#'
#' @param gene_id Character string specifying the gene ID of interest.
#' @param tx_id Character string specifying the transcript ID to be used as the main isoform.
#' @param Extend Numeric value specifying how many base pairs to extend the plot beyond the gene range. Default is 100.
#' @param NAME Optional. A character string for an additional name/title to display in the plot.
#' @param RNAcoverline Color for the RNA-seq coverage line. Default is "grey".
#' @param RNAbackground Color for the RNA-seq coverage background bars. Default is "#FEFEAE".
#' @param RNAseq A list where each element is either a character string (BAM file path) or a list with type ("bam", "bigwig", "bedgraph") and data (file for BAM, plus/minus for bigwig/bedgraph).
#' @param SampleNames Vector of sample names corresponding to the RNAseq data.
#' @param GRangeInfo A genomic range information object (e.g., `Txome_Range`) containing annotations.
#' @param RNAseqBamPaired A vector indicating whether each RNA-seq BAM file is paired-end ("paired") or single-end ("single"), used only if RNAseq contains BAM paths.
#' @param Y_scale Character string, either "all" or "each", specifying how to scale the Y-axis for RNA-seq coverage. Default is "all".
#' @param RNA_fix_height Numeric to fix the max RNA-seq coverage height. Default is \code{NULL}.
#' @param plot_ORF_ranges Logical indicating whether to plot ORF ranges in the gene model. Default is FALSE.
#' @param plot_range Optional numeric vector of length two specifying a custom genomic range to plot.
#' @param show_seq Logical indicating whether to display the DNA and amino acid sequences. Default is FALSE.
#' @param FASTA A `BSgenome` object containing genomic sequences (required if show_seq=TRUE).
#' @param plot_genomic_direction Logical indicating whether to plot an arrow showing genomic direction. Default is FALSE.
#' @param dna_aa_height_ratio Numeric value to adjust the height of the DNA and amino acid sequence plot. Default is 0.5.
#' @param gene_model_height_ratio Numeric value to adjust the height of the gene model plot. If NULL, it is auto-calculated.
#' @param transcript_label_font_size Numeric controlling the font size of the transcript ID labels in the gene model plot.
#' @param selected_isoforms Optional vector of transcript IDs to plot. If provided, only these isoforms (and `tx_id`) will be shown.
#' @param nucleotide_color_scheme If "default", uses bright colors for the nucleotides in plotDNAandAA. If "colorblind", uses a color‐blind friendly palette.
#' @param rna_linewidth Numeric value to control the thickness of RNA-seq step lines. Default is \code{0.5}.
#'
#' @return A combined ggplot object displaying RNA-seq coverage, gene models, and optionally genomic sequences.
#' @export
ggRNA <- function(gene_id = NULL, tx_id = NULL, Extend = 100, NAME = "",
                  RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                  RNAseq = inputs_full$RNAseq,
                  SampleNames = Samples,
                  GRangeInfo = Txome_Range,
                  RNAseqBamPaired = RNAseqBamPairorSingle,
                  Y_scale = "all",
                  plot_ORF_ranges = TRUE,
                  plot_range = NULL,
                  show_seq = FALSE,
                  FASTA = NULL,
                  dna_aa_height_ratio = 0.5,
                  gene_model_height_ratio = NULL,
                  RNA_fix_height = NULL,
                  transcript_label_font_size = 10,
                  plot_genomic_direction = FALSE,
                  selected_isoforms = NULL,
                  nucleotide_color_scheme = "default",
                  rna_linewidth = 0.5
) {
  # Validate Y_scale parameter
  if (!(Y_scale %in% c("all", "each"))) {
    stop("Invalid Y_scale value. Please choose either 'all' or 'each'.")
  }

  # Check that RNAbackground is correctly specified
  if (length(RNAbackground) == 1) {
    RNAbackground <- rep(RNAbackground, length(SampleNames))
  } else if (length(RNAbackground) != length(SampleNames)) {
    stop("RNAbackground must be either a single color or have the same length as 'SampleNames'.")
  }

  # Ensure GRangeInfo is provided
  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }
  # Obtain gene_id and/or tx_id
  gene_tx <- get_gene_tx(gene_id, tx_id, GRangeInfo)
  gene_id <- gene_tx$gene_id
  tx_id <- gene_tx$tx_id

  # Extract transcripts for the given gene ID
  txByYFG <- GRangeInfo$txByGene[gene_id]

  # Check if any transcripts are found
  if (length(txByYFG) == 0 || length(txByYFG[[1]]) == 0) {
    stop(paste("No transcripts found for gene ID", gene_id))
  }

  # Count isoforms and ensure 'tx_name' metadata exists
  num_isoforms <- length(txByYFG[[1]])
  if (!"tx_name" %in% names(mcols(txByYFG[[1]]))) {
    stop("Transcript names ('tx_name') not found in GRangeInfo$txByGene.")
  }
  tx_names <- txByYFG[[1]]$tx_name

  # Filter isoforms if selected_isoforms is given
  if (!is.null(selected_isoforms)) {
    tx_names <- intersect(tx_names, selected_isoforms)
  }

  # Ensure main transcript is included
  if (!tx_id %in% tx_names) {
    tx_names <- c(tx_id, tx_names)
  }

  # Check which transcripts have CDS
  tx_names_in_cdsByTx <- intersect(tx_names, names(GRangeInfo$cdsByTx))
  if (length(tx_names_in_cdsByTx) == 0) {
    message("This is a noncoding gene (no annotated CDS).")
  } else {
    # If main transcript doesn't have a CDS, just note it
    if (!(tx_id %in% tx_names_in_cdsByTx)) {
      message(paste("The transcript", tx_id, "has no annotated ORF."))
      if (!(tx_id %in% tx_names)) {
        stop(paste("Transcript ID", tx_id, "not found in gene."))
      } else {
        tx_names <- unique(c(tx_id, tx_names))
      }
    }
  }

  # Extract strand and chromosome from gene annotation
  strand_info <- as.character(strand(unlist(txByYFG)))[1]
  chr <- as.character(seqnames(unlist(txByYFG)))[1]

  # Order transcripts: main first, then others sorted
  other_tx_names <- setdiff(tx_names, tx_id)
  tx_names <- c(tx_id, sort(other_tx_names))

  # Extract CDS and exon information for the chosen transcripts
  cdsByYFGtx_all <- GRangeInfo$cdsByTx
  cdsByYFGtx <- cdsByYFGtx_all[intersect(tx_names, names(cdsByYFGtx_all))]
  for (nct in tx_names) {
    if (!nct %in% names(cdsByYFGtx)) {
      cdsByYFGtx[[nct]] <- GRanges()
    }
  }

  exonByYFGtx_all <- GRangeInfo$exonsByTx
  exonByYFGtx <- exonByYFGtx_all[intersect(tx_names, names(exonByYFGtx_all))]
  for (nct in tx_names) {
    if (!nct %in% names(exonByYFGtx)) {
      exonByYFGtx[[nct]] <- GRanges()
    }
  }

  # Prepare a list of CDS ranges per transcript
  xlimCds <- list()
  for (i in seq_along(tx_names)) {
    cds <- cdsByYFGtx[[tx_names[i]]]
    if (length(cds) > 0) {
      xlimCds[[i]] <- cds
    } else {
      xlimCds[[i]] <- GRanges()
    }
  }
  names(xlimCds) <- tx_names

  # Extract UTRs if available
  isoforms_w_3UTR <- tx_names[tx_names %in% names(GRangeInfo$threeUTR)]
  threeUTRByYFGtx <- GRangeInfo$threeUTR[isoforms_w_3UTR]

  isoforms_w_5UTR <- tx_names[tx_names %in% names(GRangeInfo$fiveUTR)]
  fiveUTRByYFGtx <- GRangeInfo$fiveUTR[isoforms_w_5UTR]
  
  # --- NEW CODE: Restrict range to selected isoforms ---
  # Make a subset of the gene transcripts corresponding to the final tx_names
  txByYFG_subset <- txByYFG[[1]][ txByYFG[[1]]$tx_name %in% tx_names ]
  if (length(txByYFG_subset) == 0) {
    stop("No transcripts left after applying selected_isoforms in ggRibo().")
  }

  # Determine genomic plotting range
  if (!is.null(plot_range)) {
    # Use custom plot_range
    plot_range <- sort(plot_range)
    range_left <- plot_range[1]
    range_right <- plot_range[2]
    gene_ranges <- GRanges(seqnames=chr, ranges=IRanges(range_left, range_right), strand=strand_info)
  } else {
    # Extend beyond gene boundaries by Extend if not provided
    gene_ranges <- reduce(txByYFG_subset)
    if (length(Extend) == 1) {
      Extend_left <- Extend
      Extend_right <- Extend
    } else if (length(Extend) == 2) {
      Extend_left <- Extend[1]
      Extend_right <- Extend[2]
    } else {
      stop("Extend must be numeric length 1 or 2.")
    }

    if (strand_info == "+") {
      range_left <- min(start(gene_ranges)) - Extend_left
      range_right <- max(end(gene_ranges)) + Extend_right
    } else if (strand_info == "-") {
      range_left <- min(start(gene_ranges)) - Extend_right
      range_right <- max(end(gene_ranges)) + Extend_left
    } else {
      stop("Invalid strand information.")
    }

    gene_ranges <- GRanges(seqnames=chr, ranges=IRanges(range_left, range_right), strand=strand_info)
  }

  # Process RNA-seq coverage
  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    for (i in seq_along(RNAseq)) {
      sample_input <- RNAseq[[i]]
      if (is.character(sample_input)) {
        # Assume it's a BAM file (backward compatibility)
        if (is.null(RNAseqBamPaired)) {
          stop("RNAseqBamPaired must be provided when using BAM files.")
        }
        paired <- RNAseqBamPaired[i]
        sample_input_list <- list(type="bam", file=sample_input, paired=paired)
        coverage_vec <- get_RNAseq_coverage(sample_input_list, gene_ranges, strand_info)
        RNAseq_list[[i]] <- coverage_vec
      } else if (is.list(sample_input)) {
        coverage_vec <- get_RNAseq_coverage(sample_input, gene_ranges, strand_info)
        RNAseq_list[[i]] <- coverage_vec
      } else {
        stop("Invalid input for RNAseq: each element must be a character string (BAM file path) or a list with type and data.")
      }
    }
  }
  # Cap RNA-seq counts if RNA_fix_height is provided
  if (!is.null(RNA_fix_height)) {
      RNAseq_list <- lapply(RNAseq_list, function(vec) {
          pmin(vec, RNA_fix_height)
     })
  }

  # Determine global max coverage for scaling
  if (length(RNAseq_list)>0) {
    max_Y_global <- max(unlist(RNAseq_list), na.rm=TRUE)
  } else {
    max_Y_global <- 0
  }

  plot_list <- list()

  # Create plots for each RNA-seq sample
  if (!is.null(RNAseq)) {
    global_start <- min(start(gene_ranges))
    global_end <- max(end(gene_ranges))
    positions <- seq(global_start, global_end)

    for (i in seq_len(length(RNAseq))) {
      RNAseq_counts <- RNAseq_list[[i]]
      RNAseq_df <- data.frame(position=positions, count=RNAseq_counts, row.names=NULL)
      RNAseq_df <- RNAseq_df[!is.na(RNAseq_df$count), ]
      RNAseq_df$isoform <- tx_id

      # Determine scaling based on Y_scale
      current_max_Y <- if (Y_scale=="all") max_Y_global else max(RNAseq_counts, na.rm=TRUE)
      # No Ribo-seq data here, so no second axis scaling needed
      scale_factor_Ribo <- 1
      y_limits <- c(0,current_max_Y*1.1)

      # Start plotting
      p <- ggplot() +
        geom_col(data=RNAseq_df, aes(x=position, y=count), fill=RNAbackground[i], color=RNAbackground[i], na.rm=TRUE) +
        geom_step(data=RNAseq_df, aes(x=position, y=count), linewidth = rna_linewidth, color=RNAcoverline, na.rm=TRUE) +
        theme_bw() +
        theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none",
          plot.margin=unit(c(0,0.2,0,0.2),"lines"),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_line(color="lightgrey",linewidth=0.3),
          axis.title.y=element_text(size=10),
          panel.background=element_rect(fill="white",color=NA)
        )

      # Adjust x-axis based on strand
      if (strand_info=="-") {
        p <- p + scale_x_reverse(limits=c(range_right,range_left))
        x_limits <- c(range_right,range_left)
      } else {
        p <- p + scale_x_continuous(limits=c(range_left,range_right))
        x_limits <- c(range_left,range_right)
      }

      p <- p + xlab("")

      # If main transcript has a CDS, show vertical lines for start/stop
      main_cds <- xlimCds[[tx_id]]
      if (length(main_cds) > 0) {
        cds_left <- min(start(main_cds))
        cds_right <- max(end(main_cds))
        x_min <- min(x_limits)
        x_max <- max(x_limits)

        main_orf_start <- if (strand_info=="+") cds_left else cds_right
        main_orf_stop <- if (strand_info=="+") cds_right else cds_left

        if (!is.na(main_orf_start) && main_orf_start>=x_min && main_orf_start<=x_max) {
          p <- p + geom_vline(xintercept=main_orf_start, linetype="dashed", color="black")
        }
        if (!is.na(main_orf_stop) && main_orf_stop>=x_min && main_orf_stop<=x_max) {
          p <- p + geom_vline(xintercept=main_orf_stop, linetype="dashed", color="darkgrey")
        }
      }

      # Set y-axis limits
      p <- p + scale_y_continuous(limits=y_limits, name="RNA-seq \ncoverage")

      # Annotate sample name
      delta_x <- 0
      x_label <- if (strand_info=="-") {
        range_right - delta_x
      } else {
        range_left + delta_x
      }
      hjust_label <- 0
      y_label <- y_limits[2] * 0.9  # 90% of the upper y-limit
      p <- p + annotate("text",
                  x = x_label, y = y_label,
                  label = SampleNames[i],
                  hjust = hjust_label, vjust = 1,  # Change vjust to 1 for bottom alignment
                  size = 3, fontface = "bold")

      p <- p + theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0,0.2,-0.8,0.2),"lines"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(color="lightgrey",linewidth=0.3),
        axis.title.y=element_text(size=10)
      )

      # If requested, plot genomic direction arrow
      if (plot_genomic_direction == TRUE) {
        if (i == 1) {
          x_min <- min(x_limits)
          x_max <- max(x_limits)
          arrow_y <- y_label * 1.05
          arrow_length <- (x_max - x_min)*0.1
          if (strand_info == "+") {
            p <- p + annotate("segment",
                              x = x_max - arrow_length, xend = x_max,
                              y = arrow_y, yend = arrow_y,
                              arrow = arrow(length=unit(0.1,"inches")),
                              color="black")
          } else {
            p <- p + annotate("segment",
                              x = x_min, xend = x_min + arrow_length,
                              y = arrow_y, yend = arrow_y,
                              arrow = arrow(length=unit(0.1,"inches")),
                              color="black")
          }
        }
      }

      plot_list[[i]] <- p
    }
  }

  # If show_seq is TRUE and FASTA is provided, plot DNA/AA below coverage
  if (show_seq && !is.null(FASTA)) {
    dna_aa_plot <- plotDNAandAA(
      GeneTxInfo=Gene_info$new(
        gene_id=gene_id,
        tx_id=tx_id,
        txByGene=txByYFG,
        cdsByYFGtx=cdsByYFGtx,
        chr=chr,
        generanges=gene_ranges,
        generangesplus=gene_ranges,
        range_left=range_left,
        range_right=range_right,
        num_isoforms=length(tx_names),
        tx_names=tx_names,
        isoforms.w.3UTR=isoforms_w_3UTR,
        isoforms.w.5UTR=isoforms_w_5UTR,
        threeUTRByYFGtx=threeUTRByYFGtx,
        fiveUTRByYFGtx=fiveUTRByYFGtx,
        exonByYFGtx=exonByYFGtx,
        Extend=Extend,
        strand=strand_info,
        xlimCds=xlimCds,
        Riboseq_list=NULL,
        cds_left=ifelse(length(xlimCds[[tx_id]])>0,min(start(xlimCds[[tx_id]])),NA),
        cds_right=ifelse(length(xlimCds[[tx_id]])>0,max(end(xlimCds[[tx_id]])),NA)
      ),
      plot_range=plot_range,
      FASTA=FASTA,
      nucleotide_color_scheme = nucleotide_color_scheme
    )
  } else {
    dna_aa_plot <- NULL
  }

  # Plot the gene model at the bottom
  gene_model_plot <- plotGeneTxModel(
    GeneTxInfo = Gene_info$new(
      gene_id=gene_id,
      tx_id=tx_id,
      txByGene=txByYFG,
      cdsByYFGtx=cdsByYFGtx,
      chr=chr,
      generanges=gene_ranges,
      generangesplus=gene_ranges,
      range_left=range_left,
      range_right=range_right,
      num_isoforms=length(tx_names),
      tx_names=tx_names,
      isoforms.w.3UTR=isoforms_w_3UTR,
      isoforms.w.5UTR=isoforms_w_5UTR,
      threeUTRByYFGtx=threeUTRByYFGtx,
      fiveUTRByYFGtx=fiveUTRByYFGtx,
      exonByYFGtx=exonByYFGtx,
      Extend=Extend,
      strand=strand_info,
      xlimCds=xlimCds,
      Riboseq_list=NULL,
      cds_left=ifelse(length(xlimCds[[tx_id]])>0,min(start(xlimCds[[tx_id]])),NA),
      cds_right=ifelse(length(xlimCds[[tx_id]])>0,max(end(xlimCds[[tx_id]])),NA)
    ),
    eORFTxInfo = NULL,
    plot_ORF_ranges = plot_ORF_ranges,
    plot_range = plot_range,
    transcript_label_font_size = transcript_label_font_size
  )

  num_transcripts <- length(tx_names)
  num_datasets <- ifelse(!is.null(RNAseq),length(RNAseq),0)
  title_height <-0.2
  rna_ribo_height <-0.8

  # If gene_model_height_ratio not set, auto-calculate
  if (is.null(gene_model_height_ratio)) {
    gene_model_height_ratio <-0.2+(num_transcripts)*0.1
  }
  gene_model_height <- gene_model_height_ratio

  # If show_seq=TRUE and FASTA provided, set dna_aa_height
  if (show_seq && !is.null(FASTA)) {
    dna_aa_height <- dna_aa_height_ratio
  } else {
    dna_aa_height <-0
  }

  # Define spacer plot and height
  spacer_plot <- ggplot() + theme_void()
  spacer_height <- 0.03  # Adjust gap size as needed

# Update total height and rel_heights
  total_height_units <- title_height + (num_datasets * rna_ribo_height) + spacer_height + dna_aa_height + gene_model_height
  rel_heights <- c(
    title_height,
    rep(rna_ribo_height, num_datasets),
    spacer_height,  # Add spacer height
    dna_aa_height,
    gene_model_height
  ) / total_height_units

  # Create title plot
  title_plot <- ggplot()+
    theme_void()+
    theme(
      plot.margin=unit(c(0,0,0,0),"lines")
    )+
    annotate("text",
             x=0.5,y=0.5,
             label=paste(gene_id," ",NAME),
             hjust=0.5,vjust=0.5,
             fontface="italic",size=5)

  # Combine all: title, RNAseq coverage plots, dna/aa plot, gene model
  combined_plot <- cowplot::plot_grid(
  title_plot,
  plotlist = c(plot_list, list(spacer_plot), list(dna_aa_plot), list(gene_model_plot)),
  ncol = 1,
  align = "v",
  rel_heights = rel_heights,
  axis = "lr",
  labels = NULL,
  label_size = 10,
  label_fontface = "plain"
  )

  return(combined_plot)
}

#' Plot RNA-seq and Ribo-seq coverage for a gene
#'
#' The `ggRibo` function creates a comprehensive plot displaying RNA-seq coverage and Ribo-seq read counts
#' for a specified gene and transcript. It includes options for displaying extended ORFs (eORFs), genomic sequences,
#' and gene models, and can scale data according to various parameters.
#'
#' @param gene_id Character string specifying the gene ID of interest.
#' @param tx_id Character string specifying the transcript ID to be used as the main isoform.
#' @param eORF.tx_id Optional. Vector of eORF transcript IDs to include in the plot.
#' @param eORFRangeInfo Optional. eORF range information, e.g., an object like `eORF_Range`.
#' @param Extend Numeric value specifying the number of base pairs to extend the plot beyond the gene range. Default is 100.
#' @param NAME Optional. Character string for an additional name or title to display in the plot.
#' @param RNAcoverline Color for the RNA-seq coverage line. Default is "grey".
#' @param RNAbackground Color for the RNA-seq coverage background. Default is "#FEFEAE".
#' @param fExtend Numeric value specifying the number of nucleotides to extend the frame assignment into the 5' UTR. Default is 0.
#' @param tExtend Numeric value specifying the number of nucleotides to extend the frame assignment into the 3' UTR. Default is 0.
#' @param RNAseq List where each element is either a character string (BAM file path) or a list with type ("bam", "bigwig", "bedgraph") and data (file for BAM, plus/minus for bigwig/bedgraph).
#' @param Riboseq List where each element is either a data frame (tabular data) or a list with type ("tabular", "bigwig", "bedgraph") and data.
#' @param SampleNames Vector of sample names corresponding to the RNAseq and Riboseq data.
#' @param GRangeInfo Genomic range information, typically an object like `Txome_Range`.
#' @param RNAseqBamPaired Vector indicating whether each RNA-seq BAM file is paired-end ("paired") or single-end ("single"), used only if RNAseq contains BAM paths.
#' @param Y_scale Character string, either "all" or "each", specifying how to scale the Y-axis for RNA-seq coverage. Default is "all".
#' @param Ribo_fix_height Numeric value to fix the maximum height of Ribo-seq counts in the plot.
#' @param RNA_fix_height Numeric to fix the max RNA-seq coverage height. Default is \code{NULL}.
#' @param plot_ORF_ranges Logical indicating whether to plot ORF ranges in the gene model. Default is FALSE.
#' @param oORF_coloring Character string specifying coloring method for overlapping ORFs ("oORF_colors" or "extend_mORF").
#' @param frame_colors Named vector of colors for the reading frames (0,1,2). Default is c("0"="#FF0000", "1"="#3366FF", "2"="#009900").
#' @param plot_range Optional numeric vector specifying a custom genomic range to plot.
#' @param sample_color Vector specifying colors for each sample or "color" to use default coloring.
#' @param show_seq Logical indicating whether to display the DNA and amino acid sequence. Default is FALSE.
#' @param FASTA A `BSgenome` object with genomic sequences.
#' @param plot_genomic_direction Logical whether to plot the genomic direction arrow. Default is FALSE.
#' @param dna_aa_height_ratio Numeric value adjusting DNA/AA plot height. Default is 0.5.
#' @param gene_model_height_ratio Numeric value adjusting gene model plot height or NULL to auto-adjust.
#' @param transcript_label_font_size Numeric for transcript label font size.
#' @param data_types Vector of sample data type names for each sample (e.g., "Ribo-seq").
#' @param selected_isoforms Optional vector of transcript IDs to plot. If provided, only these isoforms plus `tx_id` are shown.
#' @param nucleotide_color_scheme If "default", uses bright colors for the nucleotides in plotDNAandAA. If "colorblind", uses a color‐blind friendly palette.
#' @param ribo_linewidth Numeric value to control the thickness of Ribo-seq read count lines. Default is \code{0.5}.
#' @param rna_linewidth Numeric value to control the thickness of RNA-seq step lines. Default is \code{0.5}.
#'
#' @return A combined ggplot object displaying RNA-seq coverage, Ribo-seq data, gene models, and optional sequences.
#' @export
ggRibo <- function(gene_id = NULL, tx_id = NULL, eORF.tx_id = NULL,
                   eORFRangeInfo = NULL, Extend = 100, NAME = "",
                   RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                   fExtend = 0,
                   tExtend = 0,
                   RNAseq = inputs_full$RNAseq,
                   Riboseq = inputs_full$Riboseq,
                   SampleNames = Samples,
                   GRangeInfo = Txome_Range,
                   RNAseqBamPaired = RNAseqBamPairorSingle,
                   Y_scale = "all",
                   Ribo_fix_height = NULL,
                   RNA_fix_height = NULL,
                   plot_ORF_ranges = TRUE,
                   oORF_coloring = "extend_mORF",
                   frame_colors = c("0"="#FF0000", "1"="#3366FF", "2"="#009900"),
                   plot_range = NULL,
                   sample_color = rep("color", length(Riboseq)),
                   show_seq = FALSE,
                   FASTA = NULL,
                   dna_aa_height_ratio = 0.5,
                   gene_model_height_ratio = NULL,
                   transcript_label_font_size = 10,
                   plot_genomic_direction = FALSE,
                   data_types = rep("Ribo-seq", length(SampleNames)),
                   selected_isoforms = NULL,
                   nucleotide_color_scheme = "default",
                   ribo_linewidth = 0.5,
                   rna_linewidth = 0.5) {

  # Validate that data_types matches number of samples
  if (length(data_types) != length(SampleNames)) {
    stop("The length of data_types must match the number of samples.")
  }

  # Validate Y_scale parameter
  if (!(Y_scale %in% c("all", "each"))) {
    stop("Invalid Y_scale value. Please choose either 'all' or 'each'.")
  }

  # Ensure RNAbackground is correctly specified (single or matching number of samples)
  if (length(RNAbackground) == 1) {
    RNAbackground <- rep(RNAbackground, length(SampleNames))
  } else if (length(RNAbackground) != length(SampleNames)) {
    stop("RNAbackground must be either a single color or a vector of the same length as 'Samples'.")
  }

  # Check that GRangeInfo is provided
  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }
  # Obtain gene_id and/or tx_id
  gene_tx <- get_gene_tx(gene_id, tx_id, GRangeInfo)
  gene_id <- gene_tx$gene_id
  tx_id <- gene_tx$tx_id

  # Handle eORF annotation if provided
  has_overlapping_ORF <- FALSE
  if (!is.null(eORF.tx_id)) {
    # If eORFRangeInfo is not provided but "eORF_Range" exists globally, use it
    if (is.null(eORFRangeInfo)) {
      if (exists("eORF_Range", envir = .GlobalEnv)) {
        eORFRangeInfo <- get("eORF_Range", envir = .GlobalEnv)
      } else {
        stop("eORFRangeInfo (e.g., eORF_Range) must be provided when eORF.tx_id is specified.")
      }
    }
    # Check that all eORF IDs are found in eORFRangeInfo
    missing_tx_ids <- setdiff(eORF.tx_id, names(eORFRangeInfo$eORFByTx))
    if (length(missing_tx_ids) > 0) {
      stop(paste("eORF Transcript IDs", paste(missing_tx_ids, collapse = ", "), "not found in eORFRangeInfo$eORFByTx."))
    }
  }

  # Retrieve transcripts for the specified gene
  txByYFG <- GRangeInfo$txByGene[gene_id]
  if (length(txByYFG) == 0 || length(txByYFG[[1]]) == 0) {
    stop(paste("No transcripts found for gene ID", gene_id))
  }

  # Extract transcript names and verify presence of tx_name metadata
  num_isoforms <- length(txByYFG[[1]])
  if (!"tx_name" %in% names(mcols(txByYFG[[1]]))) {
    stop("Transcript names ('tx_name') not found in GRangeInfo$txByGene. Please ensure 'tx_name' is a metadata column.")
  }
  tx_names <- txByYFG[[1]]$tx_name

  # Filter isoforms if user specified selected_isoforms
  if(!is.null(selected_isoforms)) {
    tx_names <- intersect(tx_names, selected_isoforms)
  }

  # Ensure main transcript is included
  if(!tx_id %in% tx_names) {
    tx_names <- c(tx_id, tx_names)
  }

  # Check which transcripts have CDS
  tx_names_in_cdsByTx <- intersect(tx_names, names(GRangeInfo$cdsByTx))
  if (length(tx_names_in_cdsByTx) == 0) {
    message("This is a noncoding gene (no annotated CDS for any isoforms). Frame will be from start of transcript.")
  } else {
    if (!(tx_id %in% tx_names_in_cdsByTx)) {
      message(paste("The transcript", tx_id, "has no annotated ORF. Frame is from the transcript start."))
      if (!(tx_id %in% tx_names)) {
        stop(paste("Transcript ID", tx_id, "not found in gene."))
      } else {
        tx_names <- unique(c(tx_id, tx_names))
      }
    }
  }

  # Extract basic gene information (strand, chr)
  strand_info <- as.character(strand(unlist(txByYFG)))[1]
  chr <- as.character(seqnames(unlist(txByYFG)))[1]

  # Order transcripts: main at top, others sorted
  other_tx_names <- setdiff(tx_names, tx_id)
  tx_names <- c(tx_id, sort(other_tx_names))

  # Extract CDS and exon info for chosen transcripts
  cdsByYFGtx_all <- GRangeInfo$cdsByTx
  cdsByYFGtx <- cdsByYFGtx_all[intersect(tx_names, names(cdsByYFGtx_all))]
  for (nct in tx_names) {
    if (!nct %in% names(cdsByYFGtx)) {
      cdsByYFGtx[[nct]] <- GRanges()
    }
  }

  exonByYFGtx_all <- GRangeInfo$exonsByTx
  exonByYFGtx <- exonByYFGtx_all[intersect(tx_names, names(exonByYFGtx_all))]
  for (nct in tx_names) {
    if (!nct %in% names(exonByYFGtx)) {
      exonByYFGtx[[nct]] <- GRanges()
    }
  }

  # Prepare CDS ranges list
  xlimCds <- list()
  for (i in seq_along(tx_names)) {
    cds <- cdsByYFGtx[[tx_names[i]]]
    if (length(cds) > 0) {
      xlimCds[[i]] <- cds
    } else {
      xlimCds[[i]] <- GRanges()
    }
  }
  names(xlimCds) <- tx_names

  # Extract UTR info if available
  isoforms_w_3UTR <- tx_names[tx_names %in% names(GRangeInfo$threeUTR)]
  threeUTRByYFGtx <- GRangeInfo$threeUTR[isoforms_w_3UTR]

  isoforms_w_5UTR <- tx_names[tx_names %in% names(GRangeInfo$fiveUTR)]
  fiveUTRByYFGtx <- GRangeInfo$fiveUTR[isoforms_w_5UTR]

  # --- NEW CODE: Restrict range to selected isoforms ---
  # Make a subset of the gene transcripts corresponding to the final tx_names
  txByYFG_subset <- txByYFG[[1]][ txByYFG[[1]]$tx_name %in% tx_names ]
  if (length(txByYFG_subset) == 0) {
    stop("No transcripts left after applying selected_isoforms in ggRibo().")
  }

  # Define plotting range (with extensions or a custom plot_range if given)
  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)
    range_left <- plot_range[1]
    range_right <- plot_range[2]
    gene_ranges <- GRanges(seqnames=chr,
                           ranges=IRanges(range_left, range_right),
                           strand=strand_info)
  } else {
    gene_ranges <- reduce(txByYFG_subset)
    if (length(Extend) == 1) {
      Extend_left <- Extend
      Extend_right <- Extend
    } else if (length(Extend) == 2) {
      Extend_left <- Extend[1]
      Extend_right <- Extend[2]
    } else {
      stop("Extend must be a numeric value or a vector of two numeric values.")
    }

    if (strand_info == "+") {
      range_left <- min(start(gene_ranges)) - Extend_left
      range_right <- max(end(gene_ranges)) + Extend_right
    } else if (strand_info == "-") {
      range_left <- min(start(gene_ranges)) - Extend_right
      range_right <- max(end(gene_ranges)) + Extend_left
    } else {
      stop("Invalid strand information.")
    }

    gene_ranges <- GRanges(seqnames=chr,
                           ranges=IRanges(range_left, range_right),
                           strand=strand_info)
  }

  # Subset Ribo-seq data to the region of interest
  Riboseq_list <- list()
  if (!is.null(Riboseq)) {
    for (i in seq_along(Riboseq)) {
      sample_input <- Riboseq[[i]]
      if (is.data.frame(sample_input)) {
        # Filter to gene range and strand (backward compatibility)
        df <- sample_input[sample_input$chr == chr &
                           sample_input$position >= range_left &
                           sample_input$position <= range_right &
                           sample_input$strand == strand_info, ]
        Riboseq_list[[i]] <- df
      } else if (is.list(sample_input)) {
        df <- get_Riboseq_data(sample_input, gene_ranges, strand_info)
        Riboseq_list[[i]] <- df
      } else {
        stop("Invalid input for Riboseq: each element must be a data frame or a list with type and data.")
      }
    }
  }

  # Identify CDS boundaries of the main transcript
  main_cds <- xlimCds[[tx_id]]
  if (length(main_cds) > 0) {
    cds_left <- min(start(main_cds))
    cds_right <- max(end(main_cds))
  } else {
    cds_left <- NA
    cds_right <- NA
  }

  # Create a Gene_info object storing details about the gene and transcripts
  GeneTxInfo <- Gene_info$new(
    gene_id=gene_id,
    tx_id=tx_id,
    txByGene=txByYFG,
    cdsByYFGtx=cdsByYFGtx,
    chr=chr,
    generanges=gene_ranges,
    generangesplus=gene_ranges,
    range_left=range_left,
    range_right=range_right,
    num_isoforms=length(tx_names),
    tx_names=tx_names,
    isoforms.w.3UTR=isoforms_w_3UTR,
    isoforms.w.5UTR=isoforms_w_5UTR,
    threeUTRByYFGtx=threeUTRByYFGtx,
    fiveUTRByYFGtx=fiveUTRByYFGtx,
    exonByYFGtx=exonByYFGtx,
    Extend=Extend,
    strand=strand_info,
    xlimCds=xlimCds,
    Riboseq_list=Riboseq_list,
    cds_left=cds_left,
    cds_right=cds_right
  )

  # Handle eORFs if provided
  if (!is.null(eORF.tx_id) && length(tx_names)>0) {
    xlim.eORF <- eORFRangeInfo$eORFByTx[eORF.tx_id]
    eORF_left <- sapply(xlim.eORF,function(gr) min(start(gr)))
    eORF_right <- sapply(xlim.eORF,function(gr) max(end(gr)))

    if (length(Riboseq_list)>0) {
      eORF_Riboseq_list <- lapply(seq_along(Riboseq_list), function(i) {
        lapply(seq_along(xlim.eORF), function(e_idx) {
          eORF_gr <- xlim.eORF[[e_idx]]
          Riboseq_list[[i]][Riboseq_list[[i]]$position >= min(start(eORF_gr)) &
                              Riboseq_list[[i]]$position <= max(end(eORF_gr)),]
        })
      })
    } else {
      eORF_Riboseq_list <- list()
    }

    eORFTxInfo <- eORF_info$new(
      eORF.tx_id = eORF.tx_id,
      eORF_Riboseq_list = eORF_Riboseq_list,
      xlim.eORF = xlim.eORF,
      eORF_left = eORF_left,
      eORF_right = eORF_right
    )

    # Check if any eORF overlaps with the main CDS
    main_cds_ranges <- GeneTxInfo$cdsByYFGtx[[tx_id]]
    if (length(main_cds_ranges)>0) {
      for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
        eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
        overlaps_CDS <- findOverlaps(eORF_ranges, main_cds_ranges)
        if (length(overlaps_CDS)>0) {
          has_overlapping_ORF <- TRUE
          break
        }
      }
    }
  } else {
    eORFTxInfo <- NULL
  }

  # Process RNA-seq coverage if provided
  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    for (i in seq_along(RNAseq)) {
      sample_input <- RNAseq[[i]]
      if (is.character(sample_input)) {
        # Assume it's a BAM file (backward compatibility)
        if (is.null(RNAseqBamPaired)) {
          stop("RNAseqBamPaired must be provided when using BAM files.")
        }
        paired <- RNAseqBamPaired[i]
        sample_input_list <- list(type="bam", file=sample_input, paired=paired)
        coverage_vec <- get_RNAseq_coverage(sample_input_list, gene_ranges, strand_info)
        RNAseq_list[[i]] <- coverage_vec
      } else if (is.list(sample_input)) {
        coverage_vec <- get_RNAseq_coverage(sample_input, gene_ranges, strand_info)
        RNAseq_list[[i]] <- coverage_vec
      } else {
        stop("Invalid input for RNAseq: each element must be a character string (BAM file path) or a list with type and data.")
      }
    }
  }

  # Cap RNA-seq counts if RNA_fix_height is specified
  if (!is.null(RNA_fix_height)) {
    RNAseq_list <- lapply(RNAseq_list, function(vec) {
      pmin(vec, RNA_fix_height)
    })
  }

  # Determine global maxima for RNAseq and Riboseq data for scaling
  if (length(RNAseq_list)>0) {
    max_Y_global <- max(unlist(RNAseq_list), na.rm=TRUE)
  } else {
    max_Y_global <- 0
  }

  # If Ribo_fix_height is given, override Y_scale behavior
  if (!is.null(Ribo_fix_height)) {
    message("Note: Y_scale parameter is disabled when Ribo_fix_height is not NULL.")
    Y_scale <- NULL
    if (length(Riboseq_list)>0) {
      Riboseq_list <- lapply(Riboseq_list, function(df) {
        df$count <- pmin(df$count,Ribo_fix_height)
        df
      })
    }
  }

  if (length(Riboseq_list)>0) {
    all_counts <- unlist(lapply(Riboseq_list, function(df) df$count))
    if (length(all_counts)>0) {
      max_P_global <- max(all_counts, na.rm=TRUE)
      max_P_plot_global <- max_P_global + (1/10)*max_P_global
    } else {
      max_P_global <-0
      max_P_plot_global<-0
    }
  } else {
    max_P_global<-0
    max_P_plot_global<-0
  }

  plot_list <- list()

  # Generate coverage plots for each sample
  if (!is.null(RNAseq)) {
    global_start <- min(start(GeneTxInfo$generangesplus))
    global_end <- max(end(GeneTxInfo$generangesplus))
    positions <- seq(global_start, global_end)

    for (i in seq_len(length(RNAseq))) {
      RNAseq_counts <- RNAseq_list[[i]]
      RNAseq_df <- data.frame(position=positions, count=RNAseq_counts, row.names=NULL)
      RNAseq_df <- RNAseq_df[!is.na(RNAseq_df$count), ]
      RNAseq_df$isoform <- tx_id

      # Extract corresponding Ribo-seq data for this sample
      if (length(Riboseq_list)>0) {
        RiboRslt <- Riboseq_list[[i]]
      } else {
        RiboRslt <- data.frame()
      }

      # Determine scaling based on Y_scale or fixed height
      if (!is.null(Ribo_fix_height)) {
        current_max_Y <- max(RNAseq_counts,na.rm=TRUE)
        scale_factor_Ribo <- if (current_max_Y==0) 1 else current_max_Y / Ribo_fix_height
        y_limits <- c(0, current_max_Y*1.1)
      } else if (Y_scale=="all") {
        current_max_Y <- max_Y_global
        current_max_P <- max_P_global
        scale_factor_Ribo <- if (current_max_P>0) max_Y_global / max_P_global else 1
        y_limits <- c(0,current_max_Y*1.1)
      } else if (Y_scale=="each") {
        current_max_Y <- max(RNAseq_counts, na.rm=TRUE)
        if (nrow(RiboRslt)>0) {
          current_max_P <- max(RiboRslt$count,na.rm=TRUE)
          scale_factor_Ribo <- if (current_max_P>0) current_max_Y / current_max_P else 1
        } else {
          scale_factor_Ribo <-1
        }
        y_limits <- c(0,current_max_Y*1.1)
      }

      # Scale Ribo-seq counts
      if (nrow(RiboRslt)>0 && !is.null(scale_factor_Ribo)) {
        RiboRslt$count_scaled <- RiboRslt$count*scale_factor_Ribo
      }

      sample_color_i <- sample_color[i]

      # Create a ggplot object for this sample’s coverage
      p <- ggplot() +
        geom_col(data=RNAseq_df, aes(x=position, y=count), fill=RNAbackground[i], color=RNAbackground[i], na.rm=TRUE)

      # Add a step line to represent RNA coverage
      RNAseq_df_line <- RNAseq_df
      if (GeneTxInfo$strand == "+") {
        RNAseq_df_line$position <- RNAseq_df_line$position - 0.5
      } else {
        RNAseq_df_line$position <- RNAseq_df_line$position + 0.5
      }
      p <- p + geom_step(data=RNAseq_df_line, aes(x=position, y=count),linewidth =rna_linewidth, color=RNAcoverline, na.rm=TRUE)

      # Basic theming
      p <- p + theme_bw() +
        theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none",
          plot.margin=unit(c(0,0.2,0,0.2),"lines"),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_line(color="lightgrey",linewidth=0.3),
          axis.title.y=element_text(size=10),
          panel.background=element_rect(fill="white",color=NA)
        )

      # Adjust x-axis based on strand
      if (GeneTxInfo$strand=="-") {
        p <- p + scale_x_reverse(limits=c(GeneTxInfo$range_right,GeneTxInfo$range_left))
        x_limits <- c(GeneTxInfo$range_right,GeneTxInfo$range_left)
      } else {
        p <- p + scale_x_continuous(limits=c(GeneTxInfo$range_left,GeneTxInfo$range_right))
        x_limits <- c(GeneTxInfo$range_left,GeneTxInfo$range_right)
      }

      p <- p + xlab("")

      # Add vertical lines for eORF boundaries
      if (!is.null(eORFTxInfo)) {
        x_min <- min(x_limits)
        x_max <- max(x_limits)
        # Add vertical lines for eORF boundaries
        for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
          eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
          eORF_left_pos <- if (length(eORF_ranges)>0) min(start(eORF_ranges)) else NA
          eORF_right_pos <- if (length(eORF_ranges)>0) max(end(eORF_ranges)) else NA

          overlaps_CDS <- FALSE
          if (length(GeneTxInfo$xlimCds[[tx_id]])>0) {
            cds_ranges <- GeneTxInfo$xlimCds[[tx_id]]
            overlap_cds <- findOverlaps(eORF_ranges, cds_ranges)
            if (length(overlap_cds)>0) {
              overlaps_CDS <- TRUE
            }
          }
          line_color <- if (overlaps_CDS) "orange" else "orange" #pink was green
          start_pos <- if (GeneTxInfo$strand=="+") eORF_left_pos else eORF_right_pos
          end_pos <- if (GeneTxInfo$strand=="+") eORF_right_pos else eORF_left_pos
          if (!is.null(start_pos) && !is.na(start_pos) && start_pos>=x_min && start_pos<=x_max) {
            p <- p + geom_vline(xintercept=start_pos, linetype="solid", color=line_color, alpha=0.5)
          }
          if (!is.null(end_pos) && !is.na(end_pos) && end_pos>=x_min && end_pos<=x_max) {
            p <- p + geom_vline(xintercept=end_pos, linetype="dashed", color=line_color, alpha=0.5)
          }
        }
      }
      
      # Check if main ORF is annotated and add vertical lines for ORF start/stop and extensions
      main_has_cds <- length(GeneTxInfo$xlimCds[[tx_id]])>0
      if (main_has_cds) {
        # Add vertical lines for main ORF start/stop
        main_orf_start <- if (GeneTxInfo$strand=="+") GeneTxInfo$cds_left else GeneTxInfo$cds_right
        main_orf_stop <- if (GeneTxInfo$strand=="+") GeneTxInfo$cds_right else GeneTxInfo$cds_left
        x_min <- min(x_limits)
        x_max <- max(x_limits)
        if (!is.na(main_orf_start) && main_orf_start>=x_min && main_orf_start<=x_max) {
          p <- p + geom_vline(xintercept=main_orf_start, linetype="dashed", color="black")
        }
        if (!is.na(main_orf_stop) && main_orf_stop>=x_min && main_orf_stop<=x_max) {
          p <- p + geom_vline(xintercept=main_orf_stop, linetype="dashed", color="darkgrey")
        }

        # Add vertical lines for extended ORF boundaries if fExtend/tExtend > 0
        if (fExtend>0) {
          fExtend_start <- if(GeneTxInfo$strand=="+") main_orf_start - fExtend else main_orf_start + fExtend
          if (!is.na(fExtend_start) && fExtend_start>=x_min && fExtend_start<=x_max) {
            p <- p + geom_vline(xintercept=fExtend_start, linetype="dashed", color="blue")
          }
        }

        if (tExtend>0) {
          tExtend_end <- if(GeneTxInfo$strand=="+") main_orf_stop + tExtend else main_orf_stop - tExtend
          if (!is.na(tExtend_end) && tExtend_end>=x_min && tExtend_end<=x_max) {
            p <- p + geom_vline(xintercept=tExtend_end, linetype="dashed", color="blue")
          }
        }
      }

      # Add Ribo-seq segments according to chosen frame assignment strategy
      # The code below handles various scenarios (no CDS, with CDS, overlapping ORFs)
      # and assigns frames or extends frames into UTRs based on user choices.
      if (nrow(RiboRslt)>0) {
        cds_ranges <- GeneTxInfo$cdsByYFGtx[[tx_id]]
        exons <- GeneTxInfo$exonByYFGtx[[tx_id]]

        # Different coloring methods: "oORF_colors", "extend_mORF", frames with extension
        # or default coding transcript frame assignment
        if (!main_has_cds) {
          # Noncoding transcript: assign frames from transcript start
          if (GeneTxInfo$strand=="+") {
            exons_sorted <- sort(exons, decreasing=FALSE)
          } else {
            exons_sorted <- sort(exons, decreasing=TRUE)
          }
          positions_all <- integer(0)
          tx_positions <- integer(0)
          cum_len <- 0
          for (exn in seq_along(exons_sorted)) {
            exon <- exons_sorted[exn]
            pos <- seq(start(exon), end(exon))
            if (GeneTxInfo$strand=="-") {
              pos <- rev(pos)
            }
            len <- length(pos)
            tx_pos <- seq_len(len) + cum_len
            positions_all <- c(positions_all, pos)
            tx_positions <- c(tx_positions, tx_pos)
            cum_len <- cum_len + len
          }
          position_df <- data.frame(position=positions_all, tx_pos=tx_positions)
          position_df$frame <- factor((position_df$tx_pos - 1) %% 3, levels=c(0,1,2))
          RiboRslt <- merge(RiboRslt, position_df[, c("position","frame")], by="position", all.x=TRUE)

          if (!is.null(Ribo_fix_height)) {
            RiboRslt$count <- pmin(RiboRslt$count, Ribo_fix_height)
            RiboRslt$count_scaled <- RiboRslt$count * scale_factor_Ribo
          }

          # Plot Ribo-seq using frame colors or single color
          if (sample_color_i=="color") {
            p <- p + geom_segment(data=RiboRslt, aes(x=position, xend=position, y=0, yend=count_scaled, color=frame), linewidth=ribo_linewidth)
            p <- p + scale_color_manual(values=frame_colors, na.value="grey")
          } else {
            p <- p + geom_segment(data=RiboRslt, aes(x=position, xend=position, y=0, yend=count_scaled), color=sample_color_i, linewidth=ribo_linewidth)
          }

        } else {
          # -----------------------------
          # MINIMAL FIX: REORDERED so fExtend/tExtend comes before "extend_mORF"
          # -----------------------------
          if (!is.null(oORF_coloring) && oORF_coloring == "oORF_colors") {

            Ribo_main <- RiboRslt
            Ribo_main <- assign_frames(Ribo_main, cds_ranges, GeneTxInfo$strand)

            if (!is.null(eORFTxInfo)) {
              Ribo_gr <- GRanges(
                seqnames = Ribo_main$chr,
                ranges = IRanges(Ribo_main$position, Ribo_main$position),
                strand = Ribo_main$strand
              )

              overlapping_orfs <- GRangesList()
              for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
                overlaps_CDS <- findOverlaps(eORF_ranges, cds_ranges)
                if (length(overlaps_CDS)>0) {
                  overlapping_orfs[[length(overlapping_orfs) + 1]] <- eORF_ranges
                }
              }

              if (length(overlapping_orfs) > 0) {
                overlapping_orfs_gr <- unlist(overlapping_orfs)
                overlaps <- findOverlaps(Ribo_gr, overlapping_orfs_gr)
                Ribo_main$region_type <- 'non_overlapping'
                Ribo_main$region_type[queryHits(overlaps)] <- 'overlapping'
              } else {
                Ribo_main$region_type <- 'non_overlapping'
              }
            } else {
              Ribo_main$region_type <- 'non_overlapping'
            }

            if (!is.null(Ribo_fix_height)) {
              Ribo_main$count <- pmin(Ribo_main$count, Ribo_fix_height)
            }
            Ribo_main$count_scaled <- Ribo_main$count * scale_factor_Ribo

            if (sample_color_i == "color") {
              p <- p + geom_segment(data=Ribo_main[Ribo_main$region_type=='non_overlapping',],
                                    aes(x=position, xend=position, y=0, yend=count_scaled), color='grey', linewidth=ribo_linewidth)
              p <- p + geom_segment(data=Ribo_main[Ribo_main$region_type=='overlapping',],
                                    aes(x=position, xend=position, y=0, yend=count_scaled, color=frame), linewidth=ribo_linewidth)

              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]
                  eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
                  if (nrow(eORF_Riboseq) > 0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count, Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    eORF_Riboseq <- assign_frames(eORF_Riboseq, eORF_ranges, GeneTxInfo$strand)
                    p <- p + geom_segment(data=eORF_Riboseq,
                                          aes(x=position, xend=position, y=0, yend=count_scaled, color=frame), linewidth=ribo_linewidth)
                  }
                }
              }
              p <- p + scale_color_manual(values=frame_colors, na.value='grey')
            } else {
              p <- p + geom_segment(data=Ribo_main,
                                    aes(x=position, xend=position, y=0, yend=count_scaled), color=sample_color_i, linewidth=ribo_linewidth)
              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]
                  if (nrow(eORF_Riboseq)>0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    p <- p + geom_segment(data=eORF_Riboseq,
                                          aes(x=position, xend=position, y=0, yend=count_scaled), color=sample_color_i, linewidth=ribo_linewidth)
                  }
                }
              }
            }

          } else if (fExtend>0 || tExtend>0) {   # <--- Moved up to come BEFORE "extend_mORF"
            Ribo_main <- RiboRslt
            if (!is.null(eORFTxInfo)) {
              Ribo_main <- exclude_eORF_reads(Ribo_main, eORFTxInfo, GeneTxInfo$strand)
            }

            Ribo_main <- assign_frames_with_extension(Ribo_main, cds_ranges, exons, fExtend, tExtend, GeneTxInfo$strand)

            if (!is.null(Ribo_fix_height)) {
              Ribo_main$count <- pmin(Ribo_main$count,Ribo_fix_height)
              Ribo_main$count_scaled <- Ribo_main$count * scale_factor_Ribo
            }

            if (sample_color_i=="color") {
              p <- p +
                geom_segment(data=Ribo_main,
                             aes(x=position,xend=position,y=0,yend=count_scaled,color=frame), linewidth=ribo_linewidth)

              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]
                  eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
                  if (nrow(eORF_Riboseq)>0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    eORF_Riboseq <- assign_frames(eORF_Riboseq, eORF_ranges, GeneTxInfo$strand)

                    p <- p +
                      geom_segment(data=eORF_Riboseq,
                                   aes(x=position,xend=position,y=0,yend=count_scaled,color=frame), linewidth=ribo_linewidth)
                  }
                }
              }

              p <- p + scale_color_manual(values=frame_colors, na.value="grey")

            } else {
              p <- p +
                geom_segment(data=Ribo_main,
                             aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i, linewidth=ribo_linewidth)
              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]
                  if (nrow(eORF_Riboseq)>0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    p <- p +
                      geom_segment(data=eORF_Riboseq,
                                   aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i, linewidth=ribo_linewidth)
                  }
                }
              }
            }

          } else if (!is.null(oORF_coloring) && oORF_coloring=="extend_mORF") {
            extended_cds_ranges <- cds_ranges

            if (!is.null(eORFTxInfo) && has_overlapping_ORF) {
              for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
                overlaps_CDS <- findOverlaps(eORF_ranges, cds_ranges)
                if (length(overlaps_CDS)>0) {
                  extended_cds_ranges <- reduce(c(extended_cds_ranges, eORF_ranges))
                }
              }
            }

            RiboRslt <- assign_frames_extended(RiboRslt, extended_cds_ranges, GeneTxInfo$strand, cds_ranges)

            if (!is.null(Ribo_fix_height)) {
              RiboRslt$count <- pmin(RiboRslt$count,Ribo_fix_height)
              RiboRslt$count_scaled <- RiboRslt$count * scale_factor_Ribo
            }

            if (sample_color_i=="color") {
              p <- p +
                geom_segment(data=RiboRslt,
                             aes(x=position,xend=position,y=0,yend=count_scaled,color=frame), linewidth=ribo_linewidth)

              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
                  eORF_Riboseq <- eORFTxInfo$eORF_Riboseq_list[[i]][[j]]

                  overlaps_CDS <- findOverlaps(eORF_ranges, cds_ranges)
                  overlaps_fiveUTR <- FALSE
                  fiveUTR_ranges <- GeneTxInfo$fiveUTRByYFGtx[[tx_id]]
                  if (!is.null(fiveUTR_ranges) && length(fiveUTR_ranges)>0) {
                    if (length(findOverlaps(eORF_ranges, fiveUTR_ranges))>0) {
                      overlaps_fiveUTR <- TRUE
                    }
                  }

                  overlaps_threeUTR <- FALSE
                  threeUTR_ranges <- GeneTxInfo$threeUTRByYFGtx[[tx_id]]
                  if (!is.null(threeUTR_ranges) && length(threeUTR_ranges)>0) {
                    if (length(findOverlaps(eORF_ranges, threeUTR_ranges))>0) {
                      overlaps_threeUTR <- TRUE
                    }
                  }
                  
                  if ((length(overlaps_CDS)==0 && overlaps_fiveUTR) | (length(overlaps_CDS)==0 && overlaps_threeUTR)) {
                    if (nrow(eORF_Riboseq)>0) {
                      if (!is.null(Ribo_fix_height)) {
                        eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                      }
                      eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                      eORF_Riboseq <- assign_frames(eORF_Riboseq, eORF_ranges, GeneTxInfo$strand)

                      p <- p +
                        geom_segment(data=eORF_Riboseq,
                                     aes(x=position,xend=position,y=0,yend=count_scaled,color=frame), linewidth=ribo_linewidth)
                    }
                  }
                }
              }

              p <- p + scale_color_manual(values=frame_colors, na.value='grey')

            } else {
              p <- p +
                geom_segment(data=RiboRslt,
                             aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i, linewidth=ribo_linewidth)
              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORFTxInfo$eORF_Riboseq_list[[i]][[j]]
                  if (nrow(eORF_Riboseq)>0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    p <- p +
                      geom_segment(data=eORF_Riboseq,
                                   aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i, linewidth=ribo_linewidth)
                  }
                }
              }
            }

          } else {
            # Default coding transcripts frame assignment
            Ribo_main <- assign_frames(Ribo_main, cds_ranges, GeneTxInfo$strand)
            if (frame_logic != "tx_start") {                       # keep tx‑start frames for CDS genes
            Ribo_main <- assign_frames(Ribo_main, cds_ranges, GeneTxInfo$strand)
            }

            Ribo_main <- assign_frames(Ribo_main, cds_ranges, GeneTxInfo$strand)

            if (!is.null(Ribo_fix_height)) {
              Ribo_main$count <- pmin(Ribo_main$count,Ribo_fix_height)
              Ribo_main$count_scaled <- Ribo_main$count * scale_factor_Ribo
            }

            if (sample_color_i=="color") {
              p <- p +
                geom_segment(data=Ribo_main,
                             aes(x=position,xend=position,y=0,yend=count_scaled,color=frame), linewidth=ribo_linewidth)

              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]
                  eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
                  if (nrow(eORF_Riboseq)>0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    eORF_Riboseq <- assign_frames(eORF_Riboseq, eORF_ranges, GeneTxInfo$strand)
                    p <- p +
                      geom_segment(data=eORF_Riboseq,
                                   aes(x=position,xend=position,y=0,yend=count_scaled,color=frame), linewidth=ribo_linewidth)
                  }
                }
              }

              p <- p + scale_color_manual(values=frame_colors, na.value="grey")

            } else {
              p <- p +
                geom_segment(data=Ribo_main,
                             aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i, linewidth=ribo_linewidth)
              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]
                  if (nrow(eORF_Riboseq)>0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    p <- p +
                      geom_segment(data=eORF_Riboseq,
                                   aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i, linewidth=ribo_linewidth)
                  }
                }
              }
            }
          }
        }
      }

      # Set Y-axis scale and labels
      p <- p + scale_y_continuous(
        limits=y_limits,
        name="RNA-seq \ncoverage",
        sec.axis=sec_axis(~ . / scale_factor_Ribo, name = paste0(data_types[i], "\n count"))
      )

      p <- p + xlab("")

      # Add sample name as annotation
      delta_x <- 0
      x_label <- if (GeneTxInfo$strand=="-") {
        GeneTxInfo$range_right - delta_x
      } else {
        GeneTxInfo$range_left + delta_x
      }
      y_label <- y_limits[2] * 0.95
      p <- p + annotate("text",
                  x = x_label, y = y_label,
                  label = SampleNames[i],
                  hjust = 0.15, vjust = 1,
                  size = 3, fontface = "bold")
      p <- p + theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0,0.2,-0.8,0.2),"lines"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(color="lightgrey",linewidth=0.3),
        axis.title.y=element_text(size=10)
      )

      # If requested, plot genomic direction arrow on the first plot
      if (plot_genomic_direction == TRUE){
        if (i == 1) {
          x_min <- min(x_limits)
          x_max <- max(x_limits)
          arrow_y <- y_label * 1.05
          arrow_length <- (x_max - x_min)*0.1
          if (strand_info == "+") {
            p <- p + annotate("segment",
                              x = x_max - arrow_length, xend = x_max,
                              y = arrow_y, yend = arrow_y,
                              arrow = arrow(length=unit(0.1,"inches")),
                              color="black")
          } else {
            p <- p + annotate("segment",
                              x = x_min, xend = x_min + arrow_length,
                              y = arrow_y, yend = arrow_y,
                              arrow = arrow(length=unit(0.1,"inches")),
                              color="black")
          }
        }
      }

      # Store the ggplot object for this sample
      plot_list[[i]] <- ggplotGrob(p)
    }
  }

  # If show_seq = TRUE and FASTA is provided, plot DNA/AA sequence below coverage
  if (show_seq && !is.null(FASTA)) {
    dna_aa_plot <- plotDNAandAA(
      GeneTxInfo=GeneTxInfo,
      plot_range=plot_range,
      FASTA=FASTA,
      nucleotide_color_scheme = nucleotide_color_scheme
    )
  } else {
    dna_aa_plot <- NULL
  }

  # Plot gene model using plotGeneTxModel
  gene_model_plot <- plotGeneTxModel(
    GeneTxInfo = GeneTxInfo,
    eORFTxInfo = eORFTxInfo,
    plot_ORF_ranges = plot_ORF_ranges,
    plot_range = plot_range,
    transcript_label_font_size = transcript_label_font_size
  )

  # Calculate relative heights for combined plot
  num_transcripts <- GeneTxInfo$num_isoforms
  num_datasets <- ifelse(!is.null(RNAseq),length(RNAseq),0)

  title_height <-0.2
  rna_ribo_height <-0.7

  if (is.null(gene_model_height_ratio)) {
    gene_model_height_ratio <-0.1+(num_transcripts)*0.15
  }
  gene_model_height <- gene_model_height_ratio

  if (show_seq && !is.null(FASTA)) {
    dna_aa_height <- dna_aa_height_ratio
  } else {
    dna_aa_height <-0
  }

  # Define spacer plot and height
  spacer_plot <- ggplot() + theme_void()
  spacer_height <- 0.03

  # Update total height and rel_heights
  total_height_units <- title_height + (num_datasets * rna_ribo_height) + spacer_height + dna_aa_height + gene_model_height
  rel_heights <- c(
    title_height,
    rep(rna_ribo_height, num_datasets),
    spacer_height,
    dna_aa_height,
    gene_model_height
  ) / total_height_units
  # Create a title plot
  title_plot <- ggplot()+
    theme_void()+
    theme(
      plot.margin=unit(c(0,0,0,0),"lines")
    )+
    annotate("text",
             x=0.5,y=0.5,
             label=paste(gene_id," ",NAME),
             hjust=0.5,vjust=0.5,
             fontface="italic",size=5)

  # Combine all plots: title, coverage plots, dna/aa plot, gene model
  combined_plot <- cowplot::plot_grid(
  title_plot,
  plotlist = c(plot_list, list(spacer_plot), list(dna_aa_plot), list(gene_model_plot)),
  ncol = 1,
  align = "v",
  rel_heights = rel_heights,
  axis = "lr",
  labels = NULL,
  label_size = 10,
  label_fontface = "plain"
  )
  return(combined_plot)
}


#' ggRibo_decom creates a combined visualization of RNA-Seq coverage and frame-specific Ribo-Seq counts for a specified gene and transcript.
#' It generates three separate plots corresponding to the three reading frames (0, 1, and 2) of Ribo-Seq data, optionally including reads
#' that do not fall within any annotated ORF regions. It can also display extended ORFs (eORFs), genomic sequences, and gene models.
#'
#' @param gene_id Character. The identifier for the gene of interest.
#' @param tx_id Character. The transcript identifier within the gene for which the main ORF is annotated.
#' @param eORF.tx_id Character vector, optional. Transcript identifiers for extended ORFs (eORFs) associated with the gene.
#' @param eORFRangeInfo List, optional. Contains eORF range information. Required if eORF.tx_id is specified.
#' @param Extend Numeric or numeric vector of length 2. Extends the plotting range upstream and downstream of the gene. Defaults to 100.
#' @param NAME Character. An optional label or title.
#' @param RNAcoverline Character. Color for the RNA-Seq coverage line. Defaults to "grey".
#' @param RNAbackground Character or vector of length equal to number of samples. Fill color for RNA-Seq coverage bars. Defaults to "#FEFEAE".
#' @param fExtend Numeric. Extends the ORF upstream by this many bases. Defaults to 0.
#' @param tExtend Numeric. Extends the ORF downstream by this many bases. Defaults to 0.
#' @param RNAseq List. List of RNA-Seq input descriptors (BAM file paths or named lists with plus/minus for bigWig/bedGraph).
#' @param Riboseq List. List of Ribo-Seq input descriptors (tabular data frames or named lists with plus/minus for bigWig/bedGraph).
#' @param SampleNames Character vector. Names corresponding to the samples.
#' @param GRangeInfo List. Genomic range information such as transcripts, exons, CDS, etc.
#' @param RNAseqBamPaired Character vector. Indicates if each RNA-Seq BAM is paired-end ("paired") or single-end ("single").
#' @param Y_scale Character. Either "all" or "each", controlling the Y-axis scaling for RNA-seq coverage. Defaults to "all".
#' @param Ribo_fix_height Numeric, optional. Caps Ribo-Seq counts at a fixed height, ignoring Y_scale.
#' @param RNA_fix_height Numeric to fix the max RNA-seq coverage height. Default is \code{NULL}.
#' @param plot_ORF_ranges Logical. If TRUE, highlights annotated ORF ranges on the gene model. Defaults to FALSE.
#' @param oORF_coloring Character, optional. Method for coloring overlapping ORFs. "oORF_colors" or "extend_mORF".
#' @param frame_colors Named character vector. Colors for frames 0, 1, and 2. Defaults provided.
#' @param plot_range Numeric vector of length 2, optional. Custom genomic range to plot.
#' @param sample_color Character. If "color", uses frame-specific colors. Otherwise, uses a single color for Ribo reads.
#' @param show_seq Logical. If TRUE, displays DNA and AA sequences below the coverage plots. Defaults to FALSE.
#' @param FASTA Optional. Path to a FASTA file or BSgenome object with genomic sequences.
#' @param dna_aa_height_ratio Numeric. Adjusts DNA/AA plot height relative to gene model height. Defaults to 0.5.
#' @param gene_model_height_ratio Numeric, optional. Adjusts gene model plot height. If NULL, auto-scales.
#' @param transcript_label_font_size Numeric. Font size for transcript ID labels in gene model. Defaults to 10.
#' @param plot_genomic_direction Logical. If TRUE, draws an arrow indicating genomic direction on the top plot. Defaults to FALSE.
#' @param data_types Character vector. Describes data type(s) for samples (e.g., "Ribo-seq"). Must match SampleNames length.
#' @param plot_unassigned_reads Logical. If TRUE, plots Ribo-Seq reads not assigned to any ORF as grey segments.
#' @param selected_isoforms Optional. Vector of transcript IDs to plot. If specified, only these isoforms and tx_id are shown.
#' @param frame_logic Character. Determines how reading frames are assigned:
#'   - "tx_start": Frame 0 starts at the beginning of the transcript.
#'   - "CDS_start": Frame 0 starts at the annotated ORF start.
#'   - "CDS_extend": Frame 0 starts at the annotated ORF start and extends to both sides.
#'   Defaults to "tx_start" if no ORF is annotated, otherwise "CDS_start".
#' @param ribo_linewidth Numeric value to control the thickness of Ribo-seq read count lines. Default is \code{0.5}.
#' @param rna_linewidth Numeric value to control the thickness of RNA-seq step lines. Default is \code{0.5}.
#' @return A combined ggplot object with RNA-Seq coverage, three frame-specific Ribo-Seq plots, gene model, and optionally DNA/AA sequences.
#'
#' @export
ggRibo_decom <- function(gene_id = NULL, tx_id = NULL, eORF.tx_id = NULL,
                         eORFRangeInfo = NULL, Extend = 100, NAME = "",
                         RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                         fExtend = 0,
                         tExtend = 0,
                         RNAseq = inputs_full$RNAseq,
                         Riboseq = inputs_full$Riboseq,
                         SampleNames = Samples,
                         GRangeInfo = Txome_Range,
                         RNAseqBamPaired = RNAseqBamPairorSingle,
                         Y_scale = "all",
                         Ribo_fix_height = NULL,
                         RNA_fix_height = NULL,
                         plot_ORF_ranges = FALSE,
                         oORF_coloring = NULL,
                         frame_colors = c("0"="#FF0000", "1"="#3366FF", "2"="#009900"),
                         plot_range = NULL,
                         sample_color = "color",
                         show_seq = FALSE,
                         FASTA = NULL,
                         dna_aa_height_ratio = 0.5,
                         gene_model_height_ratio = NULL,
                         transcript_label_font_size = 10,
                         plot_genomic_direction = FALSE,
                         data_types = "Ribo-seq",
                         plot_unassigned_reads = TRUE,
                         selected_isoforms = NULL,
                         frame_logic = NULL,
                         nth_sample = 1,
                         ribo_linewidth = 0.5,
                         rna_linewidth = 0.5
) {

  # Validate Y_scale
  if (!(Y_scale %in% c("all", "each"))) {
    stop("Invalid Y_scale value. Please choose either 'all' or 'each'.")
  }

  # Ensure RNAbackground is correct length
  if (length(RNAbackground) == 1) {
    RNAbackground <- rep(RNAbackground, length(SampleNames))
  } else if (length(RNAbackground) != length(SampleNames)) {
    stop("RNAbackground must be either a single color or match the length of 'Samples'.")
  }

  # Check GRangeInfo provided
  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }

  # Obtain gene_id and/or tx_id
  gene_tx <- get_gene_tx(gene_id, tx_id, GRangeInfo)
  gene_id <- gene_tx$gene_id
  tx_id <- gene_tx$tx_id

  has_overlapping_ORF <- FALSE
  if (!is.null(eORF.tx_id)) {
    # If eORFRangeInfo not provided, try global eORF_Range, else error
    if (is.null(eORFRangeInfo)) {
      if (exists("eORF_Range", envir = .GlobalEnv)) {
        eORFRangeInfo <- get("eORF_Range", envir = .GlobalEnv)
      } else {
        stop("eORFRangeInfo must be provided when eORF.tx_id is specified.")
      }
    }
    # Check that all provided eORF IDs exist in eORFRangeInfo
    missing_tx_ids <- setdiff(eORF.tx_id, names(eORFRangeInfo$eORFByTx))
    if (length(missing_tx_ids) > 0) {
      stop(paste("eORF Transcript IDs", paste(missing_tx_ids, collapse = ", "),
                 "not found in eORFRangeInfo$eORFByTx."))
    }
  }

  # Get transcripts for the gene
  txByYFG <- GRangeInfo$txByGene[gene_id]
  if (length(txByYFG) == 0 || length(txByYFG[[1]]) == 0) {
    stop(paste("No transcripts found for gene ID", gene_id))
  }

  num_isoforms <- length(txByYFG[[1]])
  if (!"tx_name" %in% names(mcols(txByYFG[[1]]))) {
    stop("Transcript names ('tx_name') not found in GRangeInfo$txByGene.")
  }
  tx_names <- txByYFG[[1]]$tx_name

  # Filter isoforms if selected_isoforms given
  if (!is.null(selected_isoforms)) {
    tx_names <- intersect(tx_names, selected_isoforms)
  }

  # Ensure main transcript included
  if (!tx_id %in% tx_names) {
    tx_names <- c(tx_id, tx_names)
  }

  tx_names_in_cdsByTx <- intersect(tx_names, names(GRangeInfo$cdsByTx))
  if (length(tx_names_in_cdsByTx) == 0) {
    message("Noncoding gene: no annotated CDS for all isoforms. Frame from start of transcript.")
  } else {
    if (!(tx_id %in% tx_names_in_cdsByTx)) {
      message(paste("Transcript", tx_id, "has no annotated ORF. Frame from transcript start."))
      if (!(tx_id %in% tx_names)) {
        stop(paste("Transcript ID", tx_id, "not found in gene."))
      } else {
        tx_names <- unique(c(tx_id, tx_names))
      }
    }
  }

  RNAseq <- RNAseq[nth_sample]
  Riboseq <- Riboseq[nth_sample]
  SampleNames <- SampleNames[nth_sample]
  RNAseqBamPaired <- RNAseqBamPaired[nth_sample]

  strand_info <- as.character(strand(unlist(txByYFG)))[1]
  chr <- as.character(seqnames(unlist(txByYFG)))[1]

  # Order transcripts: main first, others sorted
  other_tx_names <- setdiff(tx_names, tx_id)
  tx_names <- c(tx_id, sort(other_tx_names))

  # Extract CDS and exon info
  cdsByYFGtx_all <- GRangeInfo$cdsByTx
  cdsByYFGtx <- cdsByYFGtx_all[intersect(tx_names, names(cdsByYFGtx_all))]
  for (nct in tx_names) {
    if (!nct %in% names(cdsByYFGtx)) {
      cdsByYFGtx[[nct]] <- GRanges()
    }
  }

  exonByYFGtx_all <- GRangeInfo$exonsByTx
  exonByYFGtx <- exonByYFGtx_all[intersect(tx_names, names(exonByYFGtx_all))]
  for (nct in tx_names) {
    if (!nct %in% names(exonByYFGtx)) {
      exonByYFGtx[[nct]] <- GRanges()
    }
  }

  # Build xlimCds list
  xlimCds <- list()
  for (i in seq_along(tx_names)) {
    cds <- cdsByYFGtx[[tx_names[i]]]
    if (length(cds) > 0) {
      xlimCds[[i]] <- cds
    } else {
      xlimCds[[i]] <- GRanges()
    }
  }
  names(xlimCds) <- tx_names

  isoforms_w_3UTR <- tx_names[tx_names %in% names(GRangeInfo$threeUTR)]
  threeUTRByYFGtx <- GRangeInfo$threeUTR[isoforms_w_3UTR]

  isoforms_w_5UTR <- tx_names[tx_names %in% names(GRangeInfo$fiveUTR)]
  fiveUTRByYFGtx <- GRangeInfo$fiveUTR[isoforms_w_5UTR]

  # Make a subset of the gene transcripts corresponding to the final tx_names
  txByYFG_subset <- txByYFG[[1]][ txByYFG[[1]]$tx_name %in% tx_names ]
  if (length(txByYFG_subset) == 0) {
    stop("No transcripts left after applying selected_isoforms in ggRibo_decom().")
  }

  # Determine plotting region
  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)
    range_left <- plot_range[1]
    range_right <- plot_range[2]
    gene_ranges <- GRanges(seqnames=chr, ranges=IRanges(range_left, range_right), strand=strand_info)
  } else {
    gene_ranges <- reduce(txByYFG_subset)
    if (length(Extend) == 1) {
      Extend_left <- Extend
      Extend_right <- Extend
    } else if (length(Extend) == 2) {
      Extend_left <- Extend[1]
      Extend_right <- Extend[2]
    } else {
      stop("Extend must be numeric length 1 or 2.")
    }

    if (strand_info == "+") {
      range_left <- min(start(gene_ranges)) - Extend_left
      range_right <- max(end(gene_ranges)) + Extend_right
    } else if (strand_info == "-") {
      range_left <- min(start(gene_ranges)) - Extend_right
      range_right <- max(end(gene_ranges)) + Extend_left
    } else {
      stop("Invalid strand info.")
    }

    gene_ranges <- GRanges(seqnames=chr, ranges=IRanges(range_left, range_right), strand=strand_info)
  }

  # Process Riboseq data
  Riboseq_list <- list()
  if (!is.null(Riboseq)) {
    Riboseq_list[[1]] <- get_Riboseq_data(Riboseq[[1]], gene_ranges, strand_info)
  }

  main_cds <- xlimCds[[tx_id]]
  if (length(main_cds) > 0) {
    cds_left <- min(start(main_cds))
    cds_right <- max(end(main_cds))
  } else {
    cds_left <- NA
    cds_right <- NA
  }

  GeneTxInfo <- Gene_info$new(
    gene_id=gene_id,
    tx_id=tx_id,
    txByGene=txByYFG,
    cdsByYFGtx=cdsByYFGtx,
    chr=chr,
    generanges=gene_ranges,
    generangesplus=gene_ranges,
    range_left=range_left,
    range_right=range_right,
    num_isoforms=length(tx_names),
    tx_names=tx_names,
    isoforms.w.3UTR=isoforms_w_3UTR,
    isoforms.w.5UTR=isoforms_w_5UTR,
    threeUTRByYFGtx=threeUTRByYFGtx,
    fiveUTRByYFGtx=fiveUTRByYFGtx,
    exonByYFGtx=exonByYFGtx,
    Extend=Extend,
    strand=strand_info,
    xlimCds=xlimCds,
    Riboseq_list=Riboseq_list,
    cds_left=cds_left,
    cds_right=cds_right
  )

  if (!is.null(eORF.tx_id) && length(tx_names)>0) {
    xlim.eORF <- eORFRangeInfo$eORFByTx[eORF.tx_id]
    eORF_left <- sapply(xlim.eORF,function(gr) min(start(gr)))
    eORF_right <- sapply(xlim.eORF,function(gr) max(end(gr)))

    if (length(Riboseq_list)>0) {
      eORF_Riboseq_list <- lapply(seq_along(Riboseq_list), function(i) {
        lapply(seq_along(xlim.eORF), function(e_idx) {
          eORF_gr <- xlim.eORF[[e_idx]]
          Riboseq_list[[i]][Riboseq_list[[i]]$position >= min(start(eORF_gr)) &
                              Riboseq_list[[i]]$position <= max(end(eORF_gr)),]
        })
      })
    } else {
      eORF_Riboseq_list <- list()
    }

    eORFTxInfo <- eORF_info$new(
      eORF.tx_id = eORF.tx_id,
      eORF_Riboseq_list = eORF_Riboseq_list,
      xlim.eORF = xlim.eORF,
      eORF_left = eORF_left,
      eORF_right = eORF_right
    )

    main_cds_ranges <- GeneTxInfo$cdsByYFGtx[[tx_id]]
    if (length(main_cds_ranges)>0) {
      for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
        eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
        overlaps_CDS <- findOverlaps(eORF_ranges, main_cds_ranges)
        if (length(overlaps_CDS)>0) {
          has_overlapping_ORF <- TRUE
          break
        }
      }
    }
  } else {
    eORFTxInfo <- NULL
  }

  # Determine if transcript has annotated ORF
  has_annotated_ORF <- length(main_cds) > 0

  # Set default frame_logic if not provided
  if (is.null(frame_logic)) {
    frame_logic <- if (has_annotated_ORF) "CDS_start" else "tx_start"
  } else {
    valid_logics <- c("tx_start", "CDS_start", "CDS_extend")
    if (!frame_logic %in% valid_logics) {
      stop("Invalid frame_logic. Choose from 'tx_start', 'CDS_start', 'CDS_extend'.")
    }
    if (!has_annotated_ORF && frame_logic %in% c("CDS_start", "CDS_extend")) {
      warning("Transcript has no annotated ORF. Falling back to 'tx_start'.")
      frame_logic <- "tx_start"
    }
  }

  # Process RNAseq coverage
  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    RNAseq_list[[1]] <- get_RNAseq_coverage(RNAseq[[1]], GeneTxInfo$generangesplus, strand_info)
  }

  # Cap RNA-seq counts if RNA_fix_height is provided
  if (!is.null(RNA_fix_height)) {
    RNAseq_list <- lapply(RNAseq_list, function(vec) {
        pmin(vec, RNA_fix_height)
    })
  }

  if (length(RNAseq_list)>0) {
    max_Y_global <- max(unlist(RNAseq_list), na.rm=TRUE)
  } else {
    max_Y_global <- 0
  }

  if (!is.null(Ribo_fix_height)) {
    message("Y_scale ignored because Ribo_fix_height is set.")
    Y_scale <- NULL
    if (length(Riboseq_list)>0) {
      Riboseq_list <- lapply(Riboseq_list, function(df) {
        df$count <- pmin(df$count,Ribo_fix_height)
        df
      })
    }
  }

  if (length(Riboseq_list)>0) {
    all_counts <- unlist(lapply(Riboseq_list, function(df) df$count))
    if (length(all_counts)>0) {
      max_P_global <- max(all_counts, na.rm=TRUE)
      max_P_plot_global <- max_P_global + (1/10)*max_P_global
    } else {
      max_P_global<-0
      max_P_plot_global<-0
    }
  } else {
    max_P_global<-0
    max_P_plot_global<-0
  }

  if (!is.null(RNAseq)) {
    global_start <- min(start(GeneTxInfo$generangesplus))
    global_end <- max(end(GeneTxInfo$generangesplus))
    positions <- seq(global_start, global_end)
    RNAseq_counts <- RNAseq_list[[1]]
    RNAseq_df <- data.frame(position=positions, count=RNAseq_counts, row.names=NULL)
    RNAseq_df <- RNAseq_df[!is.na(RNAseq_df$count), ]
    RNAseq_df$isoform <- tx_id

    if (length(Riboseq_list)>0) {
      RiboRslt <- Riboseq_list[[1]]
    } else {
      RiboRslt <- data.frame()
    }

    # Determine scaling
    if (!is.null(Ribo_fix_height)) {
      current_max_Y <- max(RNAseq_counts,na.rm=TRUE)
      scale_factor_Ribo <- if (current_max_Y==0) 1 else current_max_Y / Ribo_fix_height
      y_limits <- c(0, current_max_Y*1.1)
    } else if (Y_scale=="all") {
      current_max_Y <- max_Y_global
      current_max_P <- max_P_global
      scale_factor_Ribo <- if (current_max_P>0) max_Y_global / max_P_global else 1
      y_limits <- c(0,current_max_Y*1.1)
    } else if (Y_scale=="each") {
      current_max_Y <- max(RNAseq_counts, na.rm=TRUE)
      if (length(RiboRslt)>0 && nrow(RiboRslt)>0) {
        current_max_P <- max(RiboRslt$count,na.rm=TRUE)
        scale_factor_Ribo <- if (current_max_P>0) current_max_Y / current_max_P else 1
      } else {
        scale_factor_Ribo <-1
      }
      y_limits <- c(0,current_max_Y*1.1)
    }

    if (length(RiboRslt)>0 && !is.null(scale_factor_Ribo)) {
      RiboRslt$count_scaled <- RiboRslt$count*scale_factor_Ribo
    }

    sample_color_i <- sample_color
    x_limits <- if (GeneTxInfo$strand=="-") {
      c(GeneTxInfo$range_right,GeneTxInfo$range_left)
    } else {
      c(GeneTxInfo$range_left,GeneTxInfo$range_right)
    }

    main_has_cds <- length(GeneTxInfo$xlimCds[[tx_id]])>0
    cds_ranges <- GeneTxInfo$cdsByYFGtx[[tx_id]]
    exons <- GeneTxInfo$exonByYFGtx[[tx_id]]

    # Define a helper function to create frame-specific plots
    make_frame_plot <- function(RNAseq_df, Ribo_df, frame_color, frame_label, y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned) {
      RNAseq_df_line <- RNAseq_df
      if (GeneTxInfo$strand == "+") {
        RNAseq_df_line$position <- RNAseq_df_line$position - 0.5
      } else {
        RNAseq_df_line$position <- RNAseq_df_line$position + 0.5
      }

      p <- ggplot() +
        geom_col(data=RNAseq_df, aes(x=position, y=count), fill=RNAbackground[1], color=RNAbackground[1], na.rm=TRUE) +
        geom_step(data=RNAseq_df_line, aes(x=position, y=count), linewidth = rna_linewidth, color=RNAcoverline, na.rm=TRUE) +
        theme_bw() +
        theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none",
          plot.margin=unit(c(0,0.2,-0.8,0.2),"lines"),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_line(color="lightgrey",linewidth=0.3),
          axis.title.y=element_text(size=10),
          panel.background=element_rect(fill="white",color=NA)
        )

      if (GeneTxInfo$strand=="-") {
        p <- p + scale_x_reverse(limits=c(GeneTxInfo$range_right,GeneTxInfo$range_left))
        x_limits <- c(GeneTxInfo$range_right,GeneTxInfo$range_left)
      } else {
        p <- p + scale_x_continuous(limits=c(GeneTxInfo$range_left,GeneTxInfo$range_right))
        x_limits <- c(GeneTxInfo$range_left,GeneTxInfo$range_right)
      }

      p <- p + xlab("")

      # Add main ORF vertical lines if present
      if (main_has_cds) {
        main_orf_start <- if (GeneTxInfo$strand=="+") GeneTxInfo$cds_left else GeneTxInfo$cds_right
        main_orf_stop <- if (GeneTxInfo$strand=="+") GeneTxInfo$cds_right else GeneTxInfo$cds_left
        x_min <- min(x_limits)
        x_max <- max(x_limits)

        if (!is.na(main_orf_start) && main_orf_start>=x_min && main_orf_start<=x_max) {
          p <- p + geom_vline(xintercept=main_orf_start, linetype="dashed", color="black")
        }
        if (!is.na(main_orf_stop) && main_orf_stop>=x_min && main_orf_stop<=x_max) {
          p <- p + geom_vline(xintercept=main_orf_stop, linetype="dashed", color="darkgrey")
        }

        if (fExtend>0) {
          fExtend_start <- if(GeneTxInfo$strand=="+") main_orf_start - fExtend else main_orf_start + fExtend
          if (!is.na(fExtend_start) && fExtend_start>=x_min && fExtend_start<=x_max) {
            p <- p + geom_vline(xintercept=fExtend_start, linetype="dashed", color="blue")
          }
        }

        if (tExtend>0) {
          tExtend_end <- if(GeneTxInfo$strand=="+") main_orf_stop + tExtend else main_orf_stop - tExtend
          if (!is.na(tExtend_end) && tExtend_end>=x_min && tExtend_end<=x_max) {
            p <- p + geom_vline(xintercept=tExtend_end, linetype="dashed", color="blue")
          }
        }
      }

      # Plot unassigned reads if requested
      if (plot_unassigned && nrow(na_data)>0) {
        p <- p + geom_segment(data=na_data, aes(x=position, xend=position, y=0, yend=count_scaled), color="grey", linewidth=ribo_linewidth, na.rm=TRUE)
      }

      if (nrow(Ribo_df)>0) {
        p <- p + geom_segment(data=Ribo_df, aes(x=position, xend=position, y=0, yend=count_scaled), color=frame_color, linewidth=ribo_linewidth, na.rm=TRUE)
      }

      # Add eORF vertical lines if eORFTxInfo available
      if (!is.null(eORFTxInfo)) {
        x_min <- min(x_limits)
        x_max <- max(x_limits)
        cds_ranges <- GeneTxInfo$xlimCds[[GeneTxInfo$tx_id]]
        for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
          eORF_ranges <- eORFTxInfo$xlim.eORF[[j]]
          eORF_left_pos <- min(start(eORF_ranges))
          eORF_right_pos <- max(end(eORF_ranges))
          overlaps_CDS <- FALSE
          if (length(cds_ranges) > 0 && length(findOverlaps(eORF_ranges, cds_ranges)) > 0) {
            overlaps_CDS <- TRUE
          }
          line_color <- if (overlaps_CDS) "orange" else "orange"

          if (GeneTxInfo$strand=="+") {
            start_pos <- eORF_left_pos
            end_pos <- eORF_right_pos
          } else {
            start_pos <- eORF_right_pos
            end_pos <- eORF_left_pos
          }

          if (!is.na(start_pos) && start_pos>=x_min && start_pos<=x_max) {
            p <- p + geom_vline(xintercept=start_pos, linetype="solid", color=line_color, alpha=0.5)
          }
          if (!is.na(end_pos) && end_pos>=x_min && end_pos<=x_max) {
            p <- p + geom_vline(xintercept=end_pos, linetype="dashed", color=line_color, alpha=0.5)
          }
        }
      }

      p <- p + scale_y_continuous(
        limits=y_limits,
        name="RNA-seq \ncoverage",
        sec.axis=sec_axis(~ . / scale_factor_Ribo, name = paste0(data_types[1], "\n count"))
      )

      return(p)
    }

    # Assign frames based on frame_logic
    if (frame_logic == "tx_start" || (!has_annotated_ORF && frame_logic %in% c("CDS_start", "CDS_extend"))) {
      # Frames from start of transcript
      if (GeneTxInfo$strand=="+") {
        exons_sorted <- sort(exons, decreasing=FALSE)
      } else {
        exons_sorted <- sort(exons, decreasing=TRUE)
      }
      positions_all <- integer(0)
      tx_positions <- integer(0)
      cum_len <- 0
      for (exn in seq_along(exons_sorted)) {
        exon <- exons_sorted[exn]
        pos <- seq(start(exon), end(exon))
        if (GeneTxInfo$strand=="-") {
          pos <- rev(pos)
        }
        len <- length(pos)
        tx_pos <- seq_len(len) + cum_len
        positions_all <- c(positions_all, pos)
        tx_positions <- c(tx_positions, tx_pos)
        cum_len <- cum_len + len
      }
      position_df <- data.frame(position=positions_all, tx_pos=tx_positions)
      position_df$frame <- factor((position_df$tx_pos - 1) %% 3, levels=c(0,1,2))
      RiboRslt <- merge(RiboRslt, position_df[, c("position","frame")], by="position", all.x=TRUE)
    } else if (frame_logic == "CDS_start") {
      # Frames from CDS start, only within CDS
      RiboRslt <- assign_frames(RiboRslt, cds_ranges, GeneTxInfo$strand)
    } else if (frame_logic == "CDS_extend") {
      # Frames from CDS start, extended to entire transcript
      if (GeneTxInfo$strand == "+") {
        exons_sorted <- sort(exons, decreasing=FALSE)
      } else {
        exons_sorted <- sort(exons, decreasing=TRUE)
      }
      positions_all <- integer(0)
      tx_positions <- integer(0)
      cum_len <- 0
      for (exn in seq_along(exons_sorted)) {
        exon <- exons_sorted[exn]
        pos <- seq(start(exon), end(exon))
        if (GeneTxInfo$strand=="-") {
          pos <- rev(pos)
        }
        len <- length(pos)
        tx_pos <- seq_len(len) + cum_len
        positions_all <- c(positions_all, pos)
        tx_positions <- c(tx_positions, tx_pos)
        cum_len <- cum_len + len
      }
      position_df <- data.frame(position=positions_all, tx_pos=tx_positions)

      # Find the transcript position of the CDS start
      if (GeneTxInfo$strand == "+") {
        cds_start_genomic <- min(start(cds_ranges))
      } else {
        cds_start_genomic <- max(end(cds_ranges))
      }
      cds_start_tx <- position_df$tx_pos[position_df$position == cds_start_genomic][1]

      # Compute frames relative to CDS start
      position_df$frame <- factor((position_df$tx_pos - cds_start_tx) %% 3, levels=c(0,1,2))
      RiboRslt <- merge(RiboRslt, position_df[, c("position","frame")], by="position", all.x=TRUE)
    }

    if (!is.null(Ribo_fix_height)) {
      RiboRslt$count <- pmin(RiboRslt$count,Ribo_fix_height)
      RiboRslt$count_scaled <- RiboRslt$count * scale_factor_Ribo
    }

    # Split reads by frame
    frame0_data <- RiboRslt[RiboRslt$frame=="0", ]
    frame1_data <- RiboRslt[RiboRslt$frame=="1", ]
    frame2_data <- RiboRslt[RiboRslt$frame=="2", ]
    na_data <- RiboRslt[is.na(RiboRslt$frame),]

    # Create three separate frame plots
    p0 <- make_frame_plot(RNAseq_df, frame0_data, frame_colors["0"], "Frame0", y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned_reads)
    p1 <- make_frame_plot(RNAseq_df, frame1_data, frame_colors["1"], "Frame1", y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned_reads)
    p2 <- make_frame_plot(RNAseq_df, frame2_data, frame_colors["2"], "Frame2", y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned_reads)

    # Annotate frame labels
    delta_x <- 0
    x_label <- if (GeneTxInfo$strand=="-") {
      GeneTxInfo$range_right - delta_x
    } else {
      GeneTxInfo$range_left + delta_x
    }
    y_label <- y_limits[2] * 0.9

    p0 <- p0 + annotate("text",
                        x=x_label,y=y_label,
                        label="Frame 0",
                        hjust=-0.3,vjust=0.4,
                        size=3,fontface="bold")
    p1 <- p1 + annotate("text",
                        x=x_label,y=y_label,
                        label="Frame 1",
                        hjust=-0.3,vjust=0.4,
                        size=3,fontface="bold")
    p2 <- p2 + annotate("text",
                        x=x_label,y=y_label,
                        label="Frame 2",
                        hjust=-0.3,vjust=0.4,
                        size=3,fontface="bold")

    p2 <- p2 + theme(plot.margin=unit(c(0,0.2,-0.8,0.2),"lines"))

    # Add genomic direction arrow on first plot if requested
    if (plot_genomic_direction == TRUE) {
      x_min <- min(x_limits)
      x_max <- max(x_limits)
      arrow_y <- y_label * 1.05
      arrow_length <- (x_max - x_min)*0.1
      if (strand_info == "+") {
        p0 <- p0 + annotate("segment",
                            x = x_max - arrow_length, xend = x_max,
                            y = arrow_y, yend = arrow_y,
                            arrow = arrow(length=unit(0.1,"inches")),
                            color="black")
      } else {
        p0 <- p0 + annotate("segment",
                            x = x_min, xend = x_min + arrow_length,
                            y = arrow_y, yend = arrow_y,
                            arrow = arrow(length=unit(0.1,"inches")),
                            color="black")
      }
    }

    # If show_seq = TRUE and FASTA provided, plot DNA/AA below
    if (show_seq && !is.null(FASTA)) {
      dna_aa_plot <- plotDNAandAA(
        GeneTxInfo=GeneTxInfo,
        plot_range=plot_range,
        FASTA=FASTA
      )
    } else {
      dna_aa_plot <- NULL
    }

    # Plot gene model at bottom
    gene_model_plot <- plotGeneTxModel(
      GeneTxInfo = GeneTxInfo,
      eORFTxInfo = eORFTxInfo,
      plot_ORF_ranges = plot_ORF_ranges,
      plot_range = plot_range,
      transcript_label_font_size = transcript_label_font_size
    )

    title_height <-0.2
    frame_plot_height <-0.8
    if (is.null(gene_model_height_ratio)) {
      gene_model_height_ratio <-0.2+(GeneTxInfo$num_isoforms)*0.1
    }
    gene_model_height <- gene_model_height_ratio
    if (show_seq && !is.null(FASTA)) {
      dna_aa_height <- dna_aa_height_ratio
    } else {
      dna_aa_height <-0
    }

    # Define spacer plot and height
    spacer_plot <- ggplot() + theme_void()
    spacer_height <- 0.03

    # Update total height and rel_heights
    total_height_units <- title_height + (3 * frame_plot_height) + spacer_height + dna_aa_height + gene_model_height
    rel_heights <- c(
      title_height,
      rep(frame_plot_height, 3),
      spacer_height,  # Add spacer height
      dna_aa_height,
      gene_model_height
    ) / total_height_units

    title_plot <- ggplot()+
      theme_void()+
      theme(
        plot.margin=unit(c(0,0,0,0),"lines")
      )+
      annotate("text",
               x=0.5,y=0.5,
               label=paste(gene_id," ",NAME),
               hjust=0.5,vjust=0.5,
               fontface="italic",size=5)

    # Combine all: title, 3 frame plots, optional DNA/AA, gene model
    combined_plot <- cowplot::plot_grid(
      title_plot,
      p0,
      p1,
      p2,
      spacer_plot,
      dna_aa_plot,
      gene_model_plot,
      ncol=1,
      align="v",
      rel_heights=rel_heights,
      axis="lr",
      labels=NULL,
      label_size=10,
      label_fontface="plain"
    )

    return(combined_plot)
  } else {
    stop("No RNA-seq data provided.")
  }
}

#' Plot RNA-seq and Ribo-seq Coverage in Transcript Coordinates
#'
#' This function plots RNA-seq coverage and frame-specific Ribo-seq read counts
#' for a specified transcript in *transcript* coordinates (exon-spliced).
#' It optionally plots eORFs, genomic sequences, and a transcript model.
#'
#' @param gene_id Character string specifying the gene ID of interest.
#' @param tx_id Character string specifying the transcript ID to plot.
#' @param eORF.tx_id Optional. Vector of eORF transcript IDs to include in the plot.
#' @param eORFRangeInfo Optional. An object (e.g., \code{eORF_Range}) providing eORF ranges. If \code{NULL}, eORFs won't be plotted.
#' @param Extend Number of bases to extend the plotting region beyond the transcript. Default is \code{100}.
#' @param NAME Optional. Additional name or title for the plot.
#' @param RNAcoverline Color for the RNA-seq coverage step line. Default is \code{"grey"}.
#' @param RNAbackground Fill color for RNA-seq coverage bars. Default is \code{"#FEFEAE"}.
#' @param fExtend Number of nucleotides to extend into 5' UTR for frame assignment. Default is \code{0}.
#' @param tExtend Number of nucleotides to extend into 3' UTR for frame assignment. Default is \code{0}.
#' @param RNAseq A list describing RNA-seq input(s). Each entry can be a file path or coverage descriptor.
#' @param Riboseq A list describing Ribo-seq input(s). Each entry can be a file path or coverage descriptor.
#' @param SampleNames Vector of sample names. Must match \code{RNAseq} and/or \code{Riboseq} length.
#' @param GRangeInfo A \code{Txome_Range}-like object containing genome/transcript annotations.
#' @param RNAseqBamPaired Vector indicating pairing for BAM (e.g. \code{"paired"} or \code{"single"}). Not used for bigWig.
#' @param Y_scale Either \code{"all"} or \code{"each"}, controlling y-scaling across samples. Default is \code{"all"}.
#' @param Ribo_fix_height Numeric to fix the max Ribo-seq coverage height. Default is \code{NULL}.
#' @param RNA_fix_height Numeric to fix the max RNA-seq coverage height. Default is \code{NULL}.
#' @param plot_ORF_ranges Logical; if \code{TRUE}, attempt to plot eORFs in the gene model. Default is \code{TRUE}.
#' @param frame_colors Named vector of colors for reading frames 0,1,2.
#' @param sample_color Either \code{"color"} or a vector of colors for each sample, controlling how Ribo-seq reads are drawn.
#' @param show_seq Logical; if \code{TRUE}, also plot the transcript DNA/AA via \code{plotDNAandAA_tx}. Default is \code{FALSE}.
#' @param FASTA A \code{BSgenome} object with the reference sequences. Needed if \code{show_seq=TRUE}.
#' @param dna_aa_height_ratio Numeric adjusting the vertical space for DNA/AA. Default is \code{0.5}.
#' @param gene_model_height_ratio Numeric adjusting the height of the transcript model. Default is \code{1.3}.
#' @param transcript_label_font_size Numeric controlling transcript label font size in the gene model. Default is \code{10}.
#' @param plot_genomic_direction If \code{TRUE}, attempts to draw an arrow for the genomic direction in coverage plots.
#' @param data_types Vector describing the type of each sample (e.g. \code{"Ribo-seq"} or \code{"RNA-seq"}). Must match \code{SampleNames} length.
#' @param plot_range Optional numeric \code{c(start,end)} specifying the transcript coordinate range to plot.
#' @param nucleotide_color_scheme If \code{"colorblind"}, uses a color-blind-friendly palette for nucleotides in DNA/AA. Otherwise uses \code{"default"}.
#' @param oORF_coloring Coloring scheme for overlapping ORFs: \code{"extend_mORF"} (use main CDS frame for overlapping eORFs) or \code{"oORF_colors"} (use eORF-specific frames). Default is \code{"extend_mORF"}.
#' @param ribo_linewidth Numeric value to control the thickness of Ribo-seq read count lines. Default is \code{0.5}.
#' @param rna_linewidth Numeric value to control the thickness of RNA-seq step lines. Default is \code{0.5}.
#'
#' @return A combined \code{ggplot} object displaying RNA-seq coverage, Ribo-seq coverage, optional eORFs, a transcript model,
#'   and (if requested) the spliced DNA/AA sequences, all in transcript coordinates.
#'
#' @export
ggRibo_tx <- function(gene_id = NULL, tx_id = NULL, eORF.tx_id = NULL,
                      eORFRangeInfo=eORF_Range, Extend = 100, NAME = "",
                      RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                      fExtend = 0, tExtend = 0,
                      RNAseq = inputs_full$RNAseq,
                      Riboseq = inputs_full$Riboseq,
                      SampleNames = Samples,
                      GRangeInfo = Txome_Range,
                      RNAseqBamPaired = RNAseqBamPairorSingle,
                      Y_scale = "all",
                      Ribo_fix_height = NULL,
                      RNA_fix_height = NULL,
                      plot_ORF_ranges = TRUE,
                      frame_colors = c("0"="#FF0000","1"="#3366FF","2"="#009900"),
                      sample_color = rep("color", length(Riboseq)),
                      show_seq = FALSE,
                      FASTA = NULL,
                      dna_aa_height_ratio = 0.5,
                      gene_model_height_ratio = 1.3,
                      transcript_label_font_size = 10,
                      plot_genomic_direction = FALSE,
                      data_types = rep("Ribo-seq", length(SampleNames)),
                      plot_range = NULL,
                      nucleotide_color_scheme = "default",
                      oORF_coloring = "extend_mORF",
                      ribo_linewidth = 0.5,
                      rna_linewidth = 0.5)
{
  # Basic checks
  if (!is.null(eORF.tx_id) && is.null(eORFRangeInfo) && exists("eORF_Range", envir = .GlobalEnv)) {
    eORFRangeInfo <- get("eORF_Range", envir = .GlobalEnv)
  }
  if (length(data_types) != length(SampleNames)) {
    stop("Length of data_types must match length of SampleNames.")
  }
  if (!(Y_scale %in% c("all","each"))) {
    stop("Y_scale must be 'all' or 'each'.")
  }
  if (length(RNAbackground) == 1) {
    RNAbackground <- rep(RNAbackground, length(SampleNames))
  } else if (length(RNAbackground) != length(SampleNames)) {
    stop("RNAbackground must be a single color or match SampleNames length.")
  }
  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }

  # Obtain gene_id and/or tx_id
  gene_tx <- get_gene_tx(gene_id, tx_id, GRangeInfo)
  gene_id <- gene_tx$gene_id
  tx_id <- gene_tx$tx_id

  # Retrieve transcripts
  txByYFG <- GRangeInfo$txByGene[gene_id]
  if (length(txByYFG)==0 || length(txByYFG[[1]])==0) {
    stop(paste("No transcripts found for gene",gene_id))
  }
  tx_names <- txByYFG[[1]]$tx_name
  if (!tx_id %in% tx_names) {
    stop(paste("Transcript",tx_id,"not found in gene",gene_id))
  }

  # Get strand and chr from exonsByTx if gene_id is NA or txByYFG is empty
  if (is.na(gene_id) || length(txByYFG) == 0) {
    exon_gr <- GRangeInfo$exonsByTx[tx_id][[1]]
    if (length(exon_gr) == 0) {
      stop(paste("No exons found for transcript", tx_id))
    }
    strand_info <- as.character(strand(exon_gr)[1])
    chr <- as.character(seqnames(exon_gr)[1])
    txByYFG <- GRangesList(GRanges(chr, IRanges(min(start(exon_gr)), max(end(exon_gr))), strand = strand_info, tx_name = tx_id))
    names(txByYFG) <- "NA_gene"  # Placeholder
  } else {
    strand_info <- as.character(strand(unlist(txByYFG)))[1]
    chr <- as.character(seqnames(unlist(txByYFG)))[1]
  }

  # Exons
  exonByYFGtx <- GRangeInfo$exonsByTx[tx_id]
  exons       <- exonByYFGtx[[1]]
  if (length(exons)==0) {
    stop(paste("No exons for transcript",tx_id))
  }
  if (strand_info=="+") {
    exons <- sort(exons)
  } else {
    exons <- sort(exons, decreasing=TRUE)
  }

  # Build (genomic_pos -> tx_pos) map
  positions    <- integer(0)
  tx_positions <- integer(0)
  cum_len <- 0
  for (e_idx in seq_along(exons)) {
    e <- exons[e_idx]
    pos_vec <- seq(start(e), end(e))
    if (strand_info=="-") pos_vec <- rev(pos_vec)
    len_e <- length(pos_vec)
    tx_pos<- seq(cum_len+1, cum_len+len_e)
    positions    <- c(positions, pos_vec)
    tx_positions <- c(tx_positions, tx_pos)
    cum_len      <- cum_len + len_e
  }
  position_map <- data.frame(genomic_pos=positions, tx_pos=tx_positions)

  # Determine transcript plot range
  if (!is.null(plot_range)) {
    x_min <- plot_range[1]
    x_max <- plot_range[2]
  } else {
    x_min <- min(tx_positions)
    x_max <- max(tx_positions)
  }

  # define gene_ranges in genomic coords for coverage retrieval
  range_left  <- min(start(exons)) - Extend
  range_right <- max(end(exons))   + Extend
  gene_ranges <- GRanges(seqnames=chr, ranges=IRanges(range_left, range_right), strand=strand_info)

  # Convert Ribo-seq data to transcript coords
  Riboseq_list <- list()
  if (!is.null(Riboseq)) {
    for (i in seq_along(Riboseq)) {
      df <- get_Riboseq_data(Riboseq[[i]], gene_ranges, strand_info)
      merged <- merge(df, position_map, by.x="position", by.y="genomic_pos", all.x=FALSE)
      out_df <- data.frame(
        position=merged$tx_pos,
        count=merged$count,
        strand=merged$strand,
        chr=merged$chr
      )
      out_df <- out_df[out_df$position >= x_min & out_df$position <= x_max, ]
      Riboseq_list[[i]] <- out_df
    }
  }

  # Convert RNA-seq data
  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    for (i in seq_along(RNAseq)) {
      coverage_vec <- get_RNAseq_coverage(RNAseq[[i]], gene_ranges, strand_info)
      global_start <- start(gene_ranges)
      pos_vec <- seq(global_start, global_start + length(coverage_vec) - 1)
      df <- data.frame(genomic_pos=pos_vec, count=coverage_vec)
      merged <- merge(df, position_map, by="genomic_pos", all.x=FALSE)
      sums <- aggregate(count ~ tx_pos, data=merged, sum)
      sums <- sums[order(sums$tx_pos), ]
      RNAseq_list[[i]] <- sums
    }
  }

  # Cap RNA-seq counts if RNA_fix_height is provided
  if (!is.null(RNA_fix_height)) {
    RNAseq_list <- lapply(RNAseq_list, function(df) {
      df$count <- pmin(df$count, RNA_fix_height)
      return(df)
    })
  }

  # Gene_info object
  cdsByYFGtx <- GRangeInfo$cdsByTx[tx_id]
  xlimCds <- list(cdsByYFGtx[[1]])
  names(xlimCds) <- tx_id
  fiveUTRByYFGtx  <- GRangeInfo$fiveUTR[tx_id]
  threeUTRByYFGtx <- GRangeInfo$threeUTR[tx_id]

  # Fixed checks for isoforms.w.5UTR and isoforms.w.3UTR
  isoforms.w.5UTR <- if(length(unlist(fiveUTRByYFGtx)) > 0) tx_id else character(0)
  isoforms.w.3UTR <- if(length(unlist(threeUTRByYFGtx)) > 0) tx_id else character(0)

  GeneTxInfo <- Gene_info$new(
    gene_id=gene_id,
    tx_id=tx_id,
    txByGene=txByYFG,
    cdsByYFGtx=cdsByYFGtx,
    chr=chr,
    generanges=gene_ranges,
    generangesplus=gene_ranges,
    range_left=min(tx_positions),
    range_right=max(tx_positions),
    num_isoforms=1,
    tx_names=tx_id,
    isoforms.w.3UTR= isoforms.w.3UTR,
    isoforms.w.5UTR= isoforms.w.5UTR,
    threeUTRByYFGtx=threeUTRByYFGtx,
    fiveUTRByYFGtx= fiveUTRByYFGtx,
    exonByYFGtx=exonByYFGtx,
    Extend=Extend,
    strand=strand_info,
    xlimCds=xlimCds,
    Riboseq_list=Riboseq_list,
    cds_left= if(length(cdsByYFGtx[[1]])>0) min(start(cdsByYFGtx[[1]])) else NA,
    cds_right=if(length(cdsByYFGtx[[1]])>0) max(end(cdsByYFGtx[[1]])) else NA
  )

  # eORF handling (no message about global environment)
  if (!is.null(eORF.tx_id)) {
    # if user doesn't pass eORFRangeInfo, we do not forcibly load from .GlobalEnv
    if (is.null(eORFRangeInfo)) {
      warning("eORF.tx_id provided but eORFRangeInfo is NULL; skipping eORFs.")
      eORFTxInfo <- NULL
    } else {
      xlim.eORF <- eORFRangeInfo$eORFByTx[eORF.tx_id]
      eORF_Riboseq_list <- lapply(seq_along(Riboseq_list), function(i){
        lapply(xlim.eORF, function(eorf_gr){
          epos <- seq(min(start(eorf_gr)), max(end(eorf_gr)))
          txpos_eorf <- position_map$tx_pos[match(epos, position_map$genomic_pos)]
          Riboseq_list[[i]][Riboseq_list[[i]]$position %in% txpos_eorf, ]
        })
      })
      eORFTxInfo <- eORF_info$new(
        eORF.tx_id = eORF.tx_id,
        eORF_Riboseq_list = eORF_Riboseq_list,
        xlim.eORF = xlim.eORF,
        eORF_left  = sapply(xlim.eORF, function(gr) min(start(gr))),
        eORF_right = sapply(xlim.eORF, function(gr) max(end(gr)))
      )
    }
  } else {
    eORFTxInfo <- NULL
  }

  # Build coverage plots
  plot_list <- list()
  if (!is.null(RNAseq)) {
    for (i in seq_along(RNAseq)) {
      df_cov <- RNAseq_list[[i]]
      df_cov <- df_cov[df_cov$tx_pos>=x_min & df_cov$tx_pos<=x_max, ]
      RNAseq_df <- data.frame(position=df_cov$tx_pos, count=df_cov$count)
      RiboRslt  <- Riboseq_list[[i]]

      # Y-scaling
      if (Y_scale=="each") {
        current_max_Y    <- if(nrow(RNAseq_df)>0) max(RNAseq_df$count,na.rm=TRUE) else 0
        current_max_Ribo <- if(nrow(RiboRslt)>0) max(RiboRslt$count,na.rm=TRUE) else 0
        scale_factor <- if (current_max_Ribo>0) current_max_Y / current_max_Ribo else 1
        y_limits     <- c(0, current_max_Y*1.1)
      } else {
        all_rna_max  <- max(unlist(lapply(RNAseq_list,  function(xx) max(xx$count, na.rm=TRUE))),na.rm=TRUE)
        all_ribo_max <- max(unlist(lapply(Riboseq_list, function(xx) if(nrow(xx)>0) max(xx$count, na.rm=TRUE) else 0)),na.rm=TRUE)
        scale_factor <- if(all_ribo_max>0) all_rna_max / all_ribo_max else 1
        y_limits     <- c(0, all_rna_max*1.1)
      }
      if (!is.null(Ribo_fix_height)) {
        scale_factor <- if(nrow(RiboRslt)>0){
          cmY <- if(nrow(RNAseq_df)>0) max(RNAseq_df$count,na.rm=TRUE) else 1
          cmY / Ribo_fix_height
        } else 1
        y_limits <- c(0, if(nrow(RNAseq_df)>0) max(RNAseq_df$count,na.rm=TRUE)*1.1 else 1)
      }

      y_label <- y_limits[2] * 0.90  # Lower to 90% of the upper limit
      p <- ggplot() +
        geom_col(data=RNAseq_df, aes(x=position, y=count),
                 fill=RNAbackground[i], color=RNAbackground[i], na.rm=TRUE) +
        geom_step(data=RNAseq_df, aes(x=position-0.5, y=count),linewidth = rna_linewidth,
                  color=RNAcoverline, na.rm=TRUE) +
        theme_bw() +
        theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none",
          plot.margin=unit(c(0,0.2,-0.8,0.2),"lines"),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_line(color="lightgrey",linewidth=0.3),
          axis.title.y=element_text(size=10)
        ) +
        scale_x_continuous(limits=c(x_min,x_max)) +
        scale_y_continuous(limits=y_limits,
                           name="RNA-seq\ncoverage",
                           sec.axis=sec_axis(~. / scale_factor, name=data_types[i])) +
        xlab("") +
        annotate("text",
                  x = x_min, y = y_label,
                  label = SampleNames[i],
                  hjust = 0, vjust = 1,
                  size = 3, fontface = "bold")

      # Ribo-seq + frame assignment
      if (nrow(RiboRslt)>0) {
        cds_ranges <- cdsByYFGtx[[1]]
        cds_tx_positions <- integer(0)
        if (length(cds_ranges)>0) {
          if (strand_info=="+") {
            cds_ranges <- sort(cds_ranges)
          } else {
            cds_ranges <- sort(cds_ranges, decreasing=TRUE)
          }
          for (cgr in seq_along(cds_ranges)) {
            cpos <- seq(start(cds_ranges[cgr]), end(cds_ranges[cgr]))
            if (strand_info=="-") cpos <- rev(cpos)
            tpos <- position_map$tx_pos[match(cpos, position_map$genomic_pos)]
            tpos <- tpos[!is.na(tpos)]
            cds_tx_positions <- c(cds_tx_positions, tpos)
          }
        }

        # Assign frames for main CDS
        RiboRslt$main_frame <- assign_frames_tx(RiboRslt, cds_tx_positions)$frame

        # Identify main ORF positions
        main_orf_positions <- if (length(cds_tx_positions) > 0) cds_tx_positions else integer(0)

        # Identify eORF positions and determine overlaps
        if (!is.null(eORFTxInfo)) {
          eORF_list <- eORFTxInfo$eORF_Riboseq_list[[i]]
          # Determine which eORFs overlap with main ORF
          eORF_overlaps <- logical(length(eORF_list))
          if (length(main_orf_positions) > 0) {
            main_orf_start <- min(main_orf_positions)
            main_orf_end <- max(main_orf_positions)
            for (j in seq_along(eORF_list)) {
              eORF_data <- eORF_list[[j]]
              if (nrow(eORF_data) > 0) {
                eORF_start <- min(eORF_data$position)
                eORF_end <- max(eORF_data$position)
                if (eORF_start <= main_orf_end && eORF_end >= main_orf_start) {
                  eORF_overlaps[j] <- TRUE
                }
              }
            }
          }
          overlapping_eORF_positions <- unique(unlist(lapply(eORF_list[eORF_overlaps], function(df) df$position)))
          non_overlapping_eORF_positions <- unique(unlist(lapply(eORF_list[!eORF_overlaps], function(df) df$position)))
        } else {
          overlapping_eORF_positions <- integer(0)
          non_overlapping_eORF_positions <- integer(0)
        }

        # Frame assignment based on oORF_coloring
        RiboRslt$plot_frame <- factor(NA, levels = c(0,1,2))  # Default to NA (grey)

        if (oORF_coloring == "extend_mORF") {
          if (length(cds_tx_positions) > 0) {
            cds_start_tx <- min(cds_tx_positions)
            all_positions <- seq(min(RiboRslt$position), max(RiboRslt$position))
            extended_frame <- (all_positions - cds_start_tx) %% 3
            frame_df <- data.frame(position = all_positions, extended_frame = factor(extended_frame, levels = c(0,1,2)))
            RiboRslt <- merge(RiboRslt, frame_df, by = "position", all.x = TRUE)
          } else {
            RiboRslt$extended_frame <- factor(NA, levels = c(0,1,2))
          }

          # Assign plot_frame for main ORF positions
          idx_main <- RiboRslt$position %in% main_orf_positions
          RiboRslt$plot_frame[idx_main] <- RiboRslt$main_frame[idx_main]

          # Assign plot_frame for overlapping eORF positions (excluding main ORF)
          idx_overlap <- RiboRslt$position %in% overlapping_eORF_positions & !idx_main
          RiboRslt$plot_frame[idx_overlap] <- RiboRslt$extended_frame[idx_overlap]

          # Assign plot_frame for non-overlapping eORF positions (e.g., uORFs)
          if (!is.null(eORFTxInfo) && length(non_overlapping_eORF_positions) > 0) {
            for (j in which(!eORF_overlaps)) {
              eORF_data <- eORF_list[[j]]
              if (nrow(eORF_data) > 0) {
                eORF_start <- min(eORF_data$position)
                eORF_positions <- eORF_data$position
                eORF_frame <- (eORF_positions - eORF_start) %% 3
                eORF_frame_df <- data.frame(position = eORF_positions, eORF_frame = factor(eORF_frame, levels = c(0,1,2)))
                # Assign to RiboRslt where plot_frame is still NA to avoid overwriting
                idx_eORF <- RiboRslt$position %in% eORF_positions & is.na(RiboRslt$plot_frame)
                if (any(idx_eORF)) {
                  temp_df <- merge(RiboRslt[idx_eORF, ], eORF_frame_df, by = "position", all.x = TRUE)
                  RiboRslt$plot_frame[idx_eORF] <- temp_df$eORF_frame
                }
              }
            }
          }
        } else if (oORF_coloring == "oORF_colors") {
          if (!is.null(eORFTxInfo)) {
            for (j in seq_along(eORF_list)) {
              eORF_data <- eORF_list[[j]]
              if (nrow(eORF_data) > 0) {
                ref <- min(eORF_data$position)
                eORF_frame <- (eORF_data$position - ref) %% 3
                RiboRslt$plot_frame[RiboRslt$position %in% eORF_data$position] <- factor(eORF_frame, levels = c(0,1,2))
              }
            }
          }
        } else {
          stop("Invalid oORF_coloring option.")
        }

        if (!is.null(Ribo_fix_height)) {
          RiboRslt$count <- pmin(RiboRslt$count, Ribo_fix_height)
        }
        RiboRslt$count_scaled <- RiboRslt$count * scale_factor

        # Plot all reads with appropriate coloring
        if (sample_color[i] == "color") {
          p <- p + geom_segment(
            data=RiboRslt,
            aes(x=position, xend=position, y=0, yend=count_scaled, color=plot_frame),
            linewidth=ribo_linewidth, na.rm=TRUE
          ) +
            scale_color_manual(values=frame_colors, na.value="grey", drop=FALSE)
        } else {
          p <- p + geom_segment(
            data=RiboRslt,
            aes(x=position, xend=position, y=0, yend=count_scaled),
            color=sample_color[i], linewidth=ribo_linewidth, na.rm=TRUE
          )
        }

        # Possibly add dashed lines for the main CDS if wide enough
        if ((x_max - x_min) >= 50 && length(cds_tx_positions) > 0) {
          cds_start_tx <- min(cds_tx_positions)
          cds_stop_tx  <- max(cds_tx_positions)
          if (cds_start_tx >= x_min && cds_start_tx <= x_max) {
            p <- p + geom_vline(xintercept=cds_start_tx, linetype="dashed", color="black", alpha=0.5)
          }
          if (cds_stop_tx >= x_min && cds_stop_tx <= x_max) {
            p <- p + geom_vline(xintercept=cds_stop_tx, linetype="dashed", color="darkgrey", alpha=0.5)
          }
        }

        # Add vertical lines for eORF boundaries
        if (!is.null(eORFTxInfo)) {
          for (j in seq_along(eORF.tx_id)) {
            if (length(eORFRangeInfo$eORFByTx[[ eORF.tx_id[j] ]]) > 0) {
              eorf_gr <- eORFRangeInfo$eORFByTx[[ eORF.tx_id[j] ]]
              eorf_gen_start <- min(start(eorf_gr))
              eorf_gen_stop  <- max(end(eorf_gr))
              eorf_tx_start  <- position_map$tx_pos[match(eorf_gen_start, position_map$genomic_pos)]
              eorf_tx_stop   <- position_map$tx_pos[match(eorf_gen_stop,  position_map$genomic_pos)]
              if (!is.na(eorf_tx_start) && eorf_tx_start >= x_min && eorf_tx_start <= x_max) {
                p <- p + geom_vline(xintercept=eorf_tx_start, linetype="solid", color="orange", alpha=0.5)
              }
              if (!is.na(eorf_tx_stop) && eorf_tx_stop >= x_min && eorf_tx_stop <= x_max) {
                p <- p + geom_vline(xintercept=eorf_tx_stop,  linetype="dashed", color="orange", alpha=0.5)
              }
            }
          }
        }
      }
      plot_list[[i]] <- ggplotGrob(p)
    }
  }

  # Gene model
  gene_model_plot <- plotGeneTxModel_tx(
    GeneTxInfo                = GeneTxInfo,
    eORFTxInfo                = if(exists("eORFTxInfo")) eORFTxInfo else NULL,
    plot_ORF_ranges           = plot_ORF_ranges,
    transcript_label_font_size= transcript_label_font_size,
    plot_range                = plot_range
  )

  # Optional DNA/AA
  dna_aa_plot <- NULL
  if (show_seq && !is.null(FASTA)) {
    dna_aa_plot <- plotDNAandAA_tx(
      GeneTxInfo              = GeneTxInfo,
      plot_range              = if(!is.null(plot_range)) plot_range else c(x_min,x_max),
      FASTA                   = FASTA,
      nucleotide_color_scheme = nucleotide_color_scheme
    )
  }

  # Combine
  title_height      <- 0.2
  rna_ribo_height   <- 0.8
  num_samples <- length(RNAseq)
  gene_model_height <- gene_model_height_ratio * (0.3 + 0.1 * num_samples)
  if (show_seq && !is.null(FASTA)) {
    gene_model_height <- gene_model_height * 1.1
  }
  if (!is.null(eORF.tx_id)) {
    gene_model_height <- gene_model_height * 1.1
  }
  dna_aa_height     <- if(!is.null(dna_aa_plot)) dna_aa_height_ratio else 0
  spacer_height     <- if(!is.null(dna_aa_plot)) 0.03 else 0

  total_height <- title_height + (num_samples * rna_ribo_height) + spacer_height + dna_aa_height + gene_model_height
  rel_heights  <- c(title_height,
                    rep(rna_ribo_height, num_samples),
                    spacer_height,
                    dna_aa_height,
                    gene_model_height) / total_height

  title_plot <- ggplot() + theme_void() +
    annotate("text", x=0.5, y=0.5, label=paste(gene_id, NAME),
             hjust=0.5, vjust=0.5, fontface="italic", size=5)

  spacer_plot <- if(!is.null(dna_aa_plot)) ggplot() + theme_void() else NULL

  combined_plot <- cowplot::plot_grid(
  title_plot,
  plotlist = c(plot_list, list(spacer_plot), list(dna_aa_plot), list(gene_model_plot)),
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = rel_heights)
  return(combined_plot)
}
                                
#' Assign Frames in Transcript Coordinates for Main CDS
#'
#' Assigns frames to Ribo-seq reads (with transcript coordinate "position")
#' using the earliest CDS position as reference.
#'
#' @param Ribo_data A data frame with column "position" (transcript coordinate) and "count".
#' @param cds_tx_positions An integer vector of transcript positions that belong to the CDS.
#'
#' @return The same data frame with an added "frame" column (factor with levels 0,1,2).
assign_frames_tx <- function(Ribo_data, cds_tx_positions) {
  if (length(cds_tx_positions) == 0) {
    Ribo_data$frame <- factor(NA, levels = c(0, 1, 2))
    return(Ribo_data)
  }
  cds_start_tx <- min(cds_tx_positions)
  Ribo_data$frame <- factor(NA, levels = c(0, 1, 2))
  idx <- which(Ribo_data$position >= cds_start_tx)
  if (length(idx) > 0) {
    frames_mod <- (Ribo_data$position[idx] - cds_start_tx) %% 3
    Ribo_data$frame[idx] <- factor(frames_mod, levels = c(0, 1, 2))
  }
  return(Ribo_data)
}


#' Assign Frames for eORF Reads in Transcript Coordinates
#'
#' For eORF reads, assign frames using the minimum transcript coordinate in the eORF data
#' as the reference. This provides a frame assignment relative to each eORF's own start.
#'
#' @param Ribo_data A data frame with column "position" (transcript coordinate) and "count".
#'
#' @return The same data frame with an added "frame" column (factor with levels 0,1,2).
assign_frames_tx_eORF <- function(Ribo_data) {
  if (nrow(Ribo_data) == 0) {
    Ribo_data$frame <- factor(NA, levels = c(0, 1, 2))
    return(Ribo_data)
  }
  ref <- min(Ribo_data$position)
  Ribo_data$frame <- factor((Ribo_data$position - ref) %% 3, levels = c(0, 1, 2))
  return(Ribo_data)
}


#' Plot Transcript Model in Exon Coordinates
#'
#' Creates a transcript model plot showing exons, UTRs, CDS, and optional eORFs
#' in transcript coordinates (1..N). Features are drawn with appropriate borders.
#'
#' @param GeneTxInfo A `Gene_info` object containing gene/transcript information.
#' @param eORFTxInfo An optional `eORF_info` object.
#' @param plot_ORF_ranges Logical, whether to plot ORF ranges.
#' @param transcript_label_font_size Numeric to control transcript label font size.
#' @param plot_range Optional numeric vector (start, end) in transcript coordinates.
#'
#' @return A `ggplot` object representing the transcript model.
#'
#' @export
plotGeneTxModel_tx <- function(GeneTxInfo,
                               eORFTxInfo = NULL,
                               plot_ORF_ranges = TRUE,
                               transcript_label_font_size = 10,
                               plot_range = NULL) {

  tx_id   <- GeneTxInfo$tx_id
  strand  <- GeneTxInfo$strand
  exons_gr<- GeneTxInfo$exonByYFGtx[[tx_id]]
  if (length(exons_gr) == 0) {
    stop(paste("No exons found for transcript", tx_id))
  }

  if (strand == "+") {
    exons_gr <- sort(exons_gr, decreasing = FALSE)
  } else {
    exons_gr <- sort(exons_gr, decreasing = TRUE)
  }

  positions   <- integer(0)
  tx_positions<- integer(0)
  cum_len     <- 0
  for (i in seq_along(exons_gr)) {
    exon <- exons_gr[i]
    pos  <- seq(start(exon), end(exon))
    if (strand=="-") pos <- rev(pos)
    len <- length(pos)
    tx_pos <- seq(cum_len+1, cum_len+len)
    positions   <- c(positions, pos)
    tx_positions<- c(tx_positions, tx_pos)
    cum_len     <- cum_len + len
  }
  position_map <- data.frame(genomic_pos=positions, tx_pos=tx_positions)

  if (!is.null(plot_range)) {
    x_min <- plot_range[1]
    x_max <- plot_range[2]
  } else {
    x_min <- min(tx_positions)
    x_max <- max(tx_positions)
  }

  plot_data_list <- list()
  y_transcript   <- 1
  y_eorf         <- 1.2

  map_and_truncate <- function(feature_gr, feature_type, y_pos) {
    if (length(feature_gr)==0) return(NULL)
    out_list <- list()
    for (j in seq_along(feature_gr)) {
      f   <- feature_gr[j]
      fpos<- seq(start(f), end(f))
      txp <- position_map$tx_pos[match(fpos, position_map$genomic_pos)]
      txp <- txp[!is.na(txp)]
      if (length(txp)==0) next

      stx <- min(txp)
      etx <- max(txp)
      if (etx < x_min || stx > x_max) next
      stx_clamped <- max(stx, x_min)
      etx_clamped <- min(etx, x_max)

      out_list[[length(out_list)+1]] <- data.frame(
        start = stx_clamped,
        end   = etx_clamped,
        feature = feature_type,
        y = y_pos,
        start_truncated = (stx_clamped > stx),
        end_truncated   = (etx_clamped < etx)
      )
    }
    if (length(out_list)==0) return(NULL)
    do.call(rbind, out_list)
  }

  cds_gr <- GeneTxInfo$cdsByYFGtx[[tx_id]]
  fiveUTR_gr <- if (tx_id %in% names(GeneTxInfo$fiveUTRByYFGtx)) {
    unlist(GeneTxInfo$fiveUTRByYFGtx[tx_id])
  } else GRanges()
  threeUTR_gr<- if (tx_id %in% names(GeneTxInfo$threeUTRByYFGtx)) {
    unlist(GeneTxInfo$threeUTRByYFGtx[tx_id])
  } else GRanges()

  exon_df <- map_and_truncate(exons_gr, "exon", y_transcript)
  if (!is.null(exon_df)) plot_data_list[[length(plot_data_list)+1]] <- exon_df

  five_df <- map_and_truncate(fiveUTR_gr, "5' UTR", y_transcript)
  if (!is.null(five_df)) plot_data_list[[length(plot_data_list)+1]] <- five_df

  three_df<- map_and_truncate(threeUTR_gr, "3' UTR", y_transcript)
  if (!is.null(three_df)) plot_data_list[[length(plot_data_list)+1]] <- three_df

  if (length(cds_gr)>0) {
    cds_df <- map_and_truncate(cds_gr, "CDS", y_transcript)
    if (!is.null(cds_df)) plot_data_list[[length(plot_data_list)+1]] <- cds_df
  }

  if (plot_ORF_ranges && !is.null(eORFTxInfo)) {
    for (e_idx in seq_along(eORFTxInfo$eORF.tx_id)) {
      eORF_ranges <- eORFTxInfo$xlim.eORF[[e_idx]]
      #Check if eORF completely included in the transcript range
      overlap_exons <- findOverlaps(eORF_ranges, exons_gr, type="within")
      if (length(unique(queryHits(overlap_exons))) < length(eORF_ranges)) {
        next
      }
      if (length(eORF_ranges)==0) next
      overlaps_5prime <- length(findOverlaps(eORF_ranges, fiveUTR_gr))   > 0
      overlaps_3prime <- length(findOverlaps(eORF_ranges, threeUTR_gr))  > 0
      overlaps_CDS    <- length(findOverlaps(eORF_ranges, cds_gr))       > 0
      feature_label <-
        if (overlaps_5prime && overlaps_CDS)    "ouORF"
        else if (overlaps_5prime)               "uORF"
        else if (overlaps_3prime && overlaps_CDS) "odORF"
        else if (overlaps_3prime)               "dORF"
        else if (overlaps_CDS)                  "nORF"
        else                                    "ORF"

      eORF_df <- map_and_truncate(eORF_ranges, feature_label, y_eorf)
      if (!is.null(eORF_df)) {
        plot_data_list[[length(plot_data_list)+1]] <- eORF_df
      }
    }
  }

  plot_data <- do.call(rbind, plot_data_list)
  if (is.null(plot_data) || nrow(plot_data)==0) {
    stop("No features to plot within the specified range.")
  }

  # Define the desired order for the legend
  desired_order <- c("5' UTR", "CDS", "3' UTR", "uORF", "ouORF", "nORF", "ORF", "odORF", "dORF")

  # Get present features excluding "exon"
  present_features <- setdiff(unique(plot_data$feature), "exon")

  # Order legend_features: first those in desired_order, then others
  legend_features <- c(desired_order[desired_order %in% present_features],
                       setdiff(present_features, desired_order))

  # Check if any eORFs are plotted
  eorf_features <- c("uORF", "ouORF", "nORF", "ORF", "odORF", "dORF")
  has_eorfs <- any(eorf_features %in% plot_data$feature)

  plot_data$height <- 0.08
  plot_data$ymin   <- plot_data$y - plot_data$height
  plot_data$ymax   <- plot_data$y + plot_data$height

  feature_colors <- c(
    "uORF"   = "yellow",
    "ouORF"  = "#FFD700",
    "nORF"   = "orange",
    "ORF"    = "lightblue",
    "odORF"  = "#FFD700",
    "dORF"   = "yellow",
    "5' UTR" = "lightgrey",
    "CDS"    = "black",
    "3' UTR" = "white",
    "exon"   = "lightgrey"
  )

  # Assign "unknown" to any feature not in feature_colors
  known_levels <- names(feature_colors)
  plot_data$feature <- ifelse(plot_data$feature %in% known_levels, plot_data$feature, "unknown")

  # Add "unknown" to feature_colors if necessary
  if ("unknown" %in% plot_data$feature && !"unknown" %in% names(feature_colors)) {
    feature_colors <- c(feature_colors, "unknown" = "grey")
  }

  # Adjust y-axis labels and limits based on whether eORFs are plotted
  if (has_eorfs) {
    y_limits <- c(0.8, 1.3)
    y_breaks <- c(y_transcript, y_eorf)
    y_labels <- c(tx_id, "eORFs")
  } else {
    y_limits <- c(0.8, 1.1)
    y_breaks <- y_transcript
    y_labels <- tx_id
  }

  p_gene <- ggplot() +
    geom_rect(data = plot_data,
              aes(xmin = start - 0.5, xmax = end + 0.5, ymin = ymin, ymax = ymax, fill = feature),
              color = "black", linewidth = 0.5) +
    scale_fill_manual(values = feature_colors, breaks = legend_features) +
    scale_x_continuous(limits = c(x_min - 0.5, x_max + 0.5), name = "Transcript Position") +
    scale_y_continuous(limits = y_limits, breaks = y_breaks, labels = y_labels) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = transcript_label_font_size),
      axis.text.x = element_text(size = 8),
      axis.ticks.x = element_line(linewidth = 0.5),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(1, "lines"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )

  # Add white lines for truncated starts and ends
  truncated_data <- plot_data[plot_data$start_truncated | plot_data$end_truncated, ]
  if (nrow(truncated_data) > 0) {
    start_truncated <- truncated_data[truncated_data$start_truncated, ]
    if (nrow(start_truncated) > 0) {
      start_segments <- data.frame(
        x = start_truncated$start - 0.5,
        xend = start_truncated$start - 0.5,
        y = start_truncated$ymin - 0.01,
        yend = start_truncated$ymax + 0.01
      )
      p_gene <- p_gene + geom_segment(data = start_segments, aes(x = x, xend = xend, y = y, yend = yend), color = "white", linewidth = 1.5)
    }

    end_truncated <- truncated_data[truncated_data$end_truncated, ]
    if (nrow(end_truncated) > 0) {
      end_segments <- data.frame(
        x = end_truncated$end + 0.5,
        xend = end_truncated$end + 0.5,
        y = end_truncated$ymin - 0.01,
        yend = end_truncated$ymax + 0.01
      )
      p_gene <- p_gene + geom_segment(data = end_segments, aes(x = x, xend = xend, y = y, yend = yend), color = "white", linewidth = 1.5)
    }
  }

  return(p_gene)
}

#' @title Plot DNA and Amino Acid Sequences in Transcript Coordinates
#'
#' @description
#' Plots nucleotides (and their translated amino acids) for a spliced transcript
#' within a specified transcript coordinate range. Nucleotides/AA are drawn in
#' 5'->3' order, ignoring the genomic strand orientation for the final x-axis.
#'
#' @param GeneTxInfo A \code{Gene_info} object containing exons and other transcript data.
#' @param plot_range Optional numeric \code{c(start, end)} specifying the transcript coordinate range.
#'   If \code{NULL}, plots the full spliced transcript.
#' @param FASTA A \code{BSgenome} or \code{FaFile} object with the reference genome sequences.
#' @param nucleotide_color_scheme Either \code{"default"} or \code{"colorblind"} for the nucleotide palette.
#'
#' @return A \code{ggplot2} object with tiles for the DNA bases and codons,
#'   optionally labeled with base/AA letters if the region is \(\le\)201 nt.
#'
#' @export
plotDNAandAA_tx <- function(GeneTxInfo, plot_range = NULL, FASTA = NULL, nucleotide_color_scheme = "default") {
  if (is.null(FASTA)) {
    stop("FASTA must be provided as a BSgenome object.")
  }
  tx_id <- GeneTxInfo$tx_id
  exons_gr <- GeneTxInfo$exonByYFGtx[[tx_id]]
  strand <- GeneTxInfo$strand
  if (length(exons_gr) == 0) {
    stop(paste("No exons found for transcript", tx_id))
  }

  # Sort exons
  if (strand=="+") {
    exons_gr <- sort(exons_gr, decreasing=FALSE)
  } else {
    exons_gr <- sort(exons_gr, decreasing=TRUE)
  }

  # Build spliced transcript
  positions   <- integer(0)
  tx_positions<- integer(0)
  cum_len     <- 0
  exon_seqs   <- character(length(exons_gr))

  for (i in seq_along(exons_gr)) {
    e <- exons_gr[i]
    e_rng <- GRanges(seqnames=seqnames(e),
                     ranges=IRanges(start(e),end(e)),
                     strand=strand(e))
    exon_dna <- getSeq(FASTA, e_rng)
    exon_str <- paste0(exon_dna, collapse="")
    exon_len <- nchar(exon_str)
    pos_vec  <- seq(start(e), end(e))
    if (strand=="-") {
      pos_vec <- rev(pos_vec)
    }
    tx_vec <- seq(cum_len+1, cum_len+exon_len)
    positions   <- c(positions, pos_vec)
    tx_positions<- c(tx_positions, tx_vec)
    cum_len     <- cum_len+exon_len
    exon_seqs[i]<- exon_str
  }

  full_tx_dna <- paste0(exon_seqs, collapse="")
  tx_length   <- nchar(full_tx_dna)

  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)
    plot_range[1] <- max(plot_range[1],1)
    plot_range[2] <- min(plot_range[2],tx_length)
  } else {
    plot_range <- c(1, tx_length)
  }

  region_length <- plot_range[2] - plot_range[1] + 1
  suppress_labels <- (region_length>201)
  long_range_flag <- (region_length>201)

  dna_subseq <- substring(full_tx_dna, plot_range[1], plot_range[2])
  sub_positions <- seq(plot_range[1], plot_range[2])
  dna_chars <- unlist(strsplit(dna_subseq, split=""))

  dna_df <- data.frame(
    position   = sub_positions,
    nucleotide = dna_chars,
    stringsAsFactors=FALSE,
    row.names=NULL
  )

  # Color scheme
  if (tolower(nucleotide_color_scheme)=="colorblind") {
    nucleotide_colors <- c(
      "A"="#009E73",
      "T"="#D55E00",
      "C"="#0072B2",
      "G"="#F0E442",
      "N"="grey"
    )
  } else {
    nucleotide_colors <- c(
      "A"="#00FF00",
      "T"="#FF0200",
      "C"="#4747FF",
      "G"="#FFA503",
      "N"="grey"
    )
  }
  dna_df$fill_value <- dna_df$nucleotide

  font_size  <- max(min((region_length/length(dna_chars))*1.2,5),2)*1.3
  aa_font_sz <- font_size *1.3

  # Figure out main CDS
  cds_gr <- GeneTxInfo$cdsByYFGtx[[tx_id]]
  cds_positions <- integer(0)
  if (length(cds_gr)>0) {
    if (strand=="+") {
      cds_gr <- sort(cds_gr, decreasing=FALSE)
    } else {
      cds_gr <- sort(cds_gr, decreasing=TRUE)
    }
    for (cgr in seq_along(cds_gr)) {
      cpos <- seq(start(cds_gr[cgr]), end(cds_gr[cgr]))
      if (strand=="-") cpos <- rev(cpos)
      txp <- tx_positions[positions %in% cpos]
      cds_positions <- c(cds_positions, txp)
    }
  }
  if (length(cds_positions)>0) {
    cds_start_tx   <- min(cds_positions)
    annotated_frame<- ((cds_start_tx -1) %% 3)+1
  } else {
    annotated_frame<- 1
  }
  annotated_frame_zb <- as.integer(annotated_frame-1)
  if (is.na(annotated_frame_zb) || annotated_frame_zb<0 || annotated_frame_zb>2) {
    annotated_frame_zb <- 0
  }

  frames <- c(0,1,2)
  other_frames <- frames[frames != annotated_frame_zb]
  frame_order  <- c(annotated_frame_zb, sort(other_frames))
  if (long_range_flag) {
    frame_y_positions <- c(0.8,0.6,0.4)
  } else {
    frame_y_positions <- c(0.6,0.4,0.2)
  }
  names(frame_y_positions) <- frame_order
  frame_labels <- c("Annotated","+1","+2")
  names(frame_labels) <- frame_order

  # Translate in each frame
  dna_len <- nchar(dna_subseq)
  aa_list <- list()
  for (frame in frames) {
    codon_starts <- seq(frame+1, dna_len-2, by=3)
    if (length(codon_starts)==0) next
    codon_middles <- codon_starts+1
    dna_coding    <- substring(dna_subseq, codon_starts[1], codon_starts[length(codon_starts)]+2)
    aa_str  <- as.character(translate(DNAString(dna_coding),genetic.code = GENETIC_CODE, no.init.codon = TRUE, if.fuzzy.codon="X"))
    aa_chars<- unlist(strsplit(aa_str, split=""))
    aa_positions <- sub_positions[codon_middles]
    aa_df <- data.frame(
      position   = aa_positions,
      amino_acid = aa_chars,
      y          = frame_y_positions[as.character(frame)],
      frame      = frame,
      stringsAsFactors=FALSE,
      row.names=NULL
    )
    aa_list[[frame+1]] <- aa_df
  }
  aa_df_combined <- do.call(rbind, aa_list)
  if (is.null(aa_df_combined)) {
    aa_df_combined <- data.frame(position=integer(0), amino_acid=character(0),
                                 y=numeric(0), frame=integer(0),
                                 stringsAsFactors=FALSE, row.names=NULL)
  }
  aa_df_combined$frame_label <- frame_labels[as.character(aa_df_combined$frame)]
  aa_df_combined$fill_value  <- aa_df_combined$frame_label
  aa_df_combined$fill_value[aa_df_combined$amino_acid=="M"] <- "Start"
  aa_df_combined$fill_value[aa_df_combined$amino_acid=="*"] <- "Stop"

  p <- ggplot()

  if (!long_range_flag) {
    p <- p + geom_tile(
      data=dna_df,
      aes(x=position, y=0.8, fill=fill_value),
      width=1, height=0.2, color="darkgrey", linewidth=0.2,
      show.legend=FALSE, na.rm=TRUE
    )
    if (!suppress_labels) {
      p <- p + geom_text(
        data=dna_df,
        aes(x=position, y=0.8, label=nucleotide),
        size=font_size, fontface="plain", color="black",
        show.legend=FALSE, na.rm=TRUE
      )
    }
  }

  aa_ss <- subset(aa_df_combined, fill_value %in% c("Start","Stop"))
  aa_rg <- subset(aa_df_combined, !fill_value %in% c("Start","Stop"))

  p <- p +
    geom_tile(
      data=aa_rg,
      aes(x=position, y=y, fill=fill_value),
      width=3, height=0.2,
      color=ifelse(long_range_flag,NA,"darkgrey"),
      linewidth=0.2, show.legend=TRUE, na.rm=TRUE
    ) +
    geom_tile(
      data=aa_ss,
      aes(x=position, y=y, fill=fill_value),
      width=3, height=0.2,
      color=ifelse(long_range_flag,NA,"darkgrey"),
      linewidth=0.5, show.legend=TRUE, na.rm=TRUE
    )

  if (!suppress_labels) {
    p <- p + geom_text(
      data=aa_df_combined,
      aes(x=position, y=y, label=amino_acid),
      size=aa_font_sz, fontface="plain", color="black",
      show.legend=FALSE, na.rm=TRUE
    )
  }

  frame_colors <- c("Annotated"="#F1F1F1", "+1"="#E6E6E6", "+2"="#C9C9C9")
  fill_colors  <- c(nucleotide_colors, "Start"="green","Stop"="red", frame_colors)

  # Condition for breaks
  possible_breaks <- intersect(unique(aa_df_combined$fill_value), c("Start","Stop"))
  if (length(possible_breaks)==0) {
    p <- p + scale_fill_manual(
      name=NULL, values=fill_colors,
      na.value="grey",
      guide=guide_legend(override.aes=list(colour=NA)),
      breaks=NULL
    )
  } else {
    p <- p + scale_fill_manual(
      name=NULL, values=fill_colors,
      na.value="grey",
      guide=guide_legend(override.aes=list(colour=NA)),
      breaks=possible_breaks
    )
  }

  p <- p + scale_x_continuous(limits=c(plot_range[1], plot_range[2])) +
    theme_void() +
    theme(plot.margin=unit(c(-1,0.2,-0.5,0.2),"lines"))

  if (long_range_flag) {
    p <- p + scale_y_continuous(limits=c(0.2,1))
  } else {
    p <- p + scale_y_continuous(limits=c(0,1))
  }
  return(p)
}
