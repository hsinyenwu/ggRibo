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
  # Create a TxDb object from the given annotation file using txdbmaker
  txdb <- txdbmaker::makeTxDbFromGFF(annotation, format = format, dataSource = dataSource, organism = organism)

  # Extract exonic regions grouped by transcript
  exonsByTx <- exonsBy(txdb, by = 'tx', use.names = TRUE)
  # Extract transcript ranges grouped by gene
  txByGene <- transcriptsBy(txdb, by = 'gene')
  # Extract coding DNA sequence (CDS) ranges grouped by transcript
  cdsByTx <- cdsBy(txdb, by = "tx", use.names = TRUE)
  # Extract 5' UTR ranges grouped by transcript
  fiveUTR <- fiveUTRsByTranscript(txdb, use.names = TRUE)
  # Extract 3' UTR ranges grouped by transcript
  threeUTR <- threeUTRsByTranscript(txdb, use.names = TRUE)

  # Create a Range_info object to store all extracted genomic range information
  Txome_Range <- Range_info$new(
    exonsByTx = exonsByTx,
    txByGene = txByGene,
    cdsByTx = cdsByTx,
    fiveUTR = fiveUTR,
    threeUTR = threeUTR
  )

  # Assign the newly created Range_info object to a global variable Txome_Range for easy access
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
  txdb <- txdbmaker::makeTxDbFromGFF(annotation, format = format, dataSource = dataSource, organism = organism)

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

# Function to plot DNA and Amino Acid sequences
#' @title Plot DNA and Amino Acid Sequences
#' @description
#' Generates a plot showing DNA nucleotides and their corresponding amino acids for the genomic region of interest.
#'
#' @param GeneTxInfo A \code{Gene_info} object containing gene-specific information.
#' @param plot_range Optional vector specifying the genomic range to plot.
#' @param FASTA A \code{BSgenome} object containing the reference genome sequences.
#' @return A \code{ggplot2} object representing the DNA and amino acid sequences.
plotDNAandAA <- function(GeneTxInfo, plot_range = NULL, FASTA = NULL) {
  # This function extracts the DNA sequence for the specified region, translates it into amino acids
  # for all three reading frames, and then uses ggplot2 to create a visual representation of both the
  # nucleotide and amino acid sequence. If the region is large, it suppresses labels for clarity.

  # Extract the gene limits from GeneTxInfo
  genelim <- c(GeneTxInfo$range_left, GeneTxInfo$range_right)
  strand <- GeneTxInfo$strand

  # Adjust genomic limits based on chromosome length
  chrom_length <- seqlengths(FASTA)[GeneTxInfo$chr]
  genelim_adj <- pmax(pmin(genelim, chrom_length), 1)

  # If a custom plot_range is given, use that instead
  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)
    genelim_adj <- pmax(pmin(plot_range, chrom_length), 1)
  }

  # Determine range length and decide whether to suppress labels
  range_length <- abs(diff(range(genelim_adj))) + 1
  suppress_labels <- FALSE
  if (range_length > 201) {
    suppress_labels <- TRUE
  }
  long_range <- range_length > 201

  # Extract DNA sequence from FASTA within the adjusted range
  seq_region <- GRanges(
    seqnames = GeneTxInfo$chr,
    ranges = IRanges(min(genelim_adj), max(genelim_adj)),
    strand = "+"
  )
  seqs <- getSeq(FASTA, seq_region)
  dna_seq <- as.character(seqs)

  # Reverse complement if on the negative strand
  if (GeneTxInfo$strand == "-") {
    dna_seq <- as.character(reverseComplement(DNAString(dna_seq)))
  }

  # Split DNA sequence into individual nucleotides
  dna_chars <- unlist(strsplit(dna_seq, split = ""))

  # Assign positions along the sequence
  if (GeneTxInfo$strand == "+") {
    positions_seq <- seq(min(genelim_adj), max(genelim_adj))
  } else {
    positions_seq <- seq(max(genelim_adj), min(genelim_adj), by = -1)
  }

  # Prepare a data frame for nucleotides
  dna_df <- data.frame(position = positions_seq,
                       nucleotide = dna_chars,
                       stringsAsFactors = FALSE,
                       row.names = NULL)

  # Define colors for nucleotides
  nucleotide_colors <- c(
    "A" = "#00FF00",   # Green
    "T" = "#FF0200",   # Red
    "C" = "#4747FF",   # Blue
    "G" = "#FFA503",   # Orange
    "N" = "grey"
  )
  dna_df$fill_value <- dna_df$nucleotide

  # Adjust font size based on range
  num_nucleotides <- length(dna_chars)
  plot_width <- abs(diff(range(genelim_adj)))
  font_size <- (plot_width / num_nucleotides) * 1.2
  font_size <- max(min(font_size, 5), 2)
  font_size <- font_size * 1.5

  # Generate three-frame translations
  frames <- c(0, 1, 2)
  frame_colors <- c("Annotated" = "#F1F1F1", "+1" = "#E6E6E6", "+2" = "#C9C9C9")

  # Determine annotated frame based on CDS start (if available)
  cds_ranges <- GeneTxInfo$xlimCds[[GeneTxInfo$tx_id]]
  if (length(cds_ranges) > 0) {
    if (GeneTxInfo$strand == "+") {
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
    aa_seq <- suppressWarnings(as.character(translate(dna_string, if.fuzzy.codon = "X")))
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
  aa_font_size <- font_size*1.25
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
  aa_start_stop_df <- aa_df_combined[aa_df_combined$fill_value %in% c("Start","Stop"),]
  aa_regular_df <- aa_df_combined[!aa_df_combined$fill_value %in% c("Start","Stop"),]

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
  if (GeneTxInfo$strand == "-") {
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
plotGeneTxModel <- function(GeneTxInfo = GeneTxInfo, eORFTxInfo = NULL, XLIM = NULL, plot_ORF_ranges = FALSE, plot_range = NULL, transcript_label_font_size = 10) {
  # Load necessary libraries
  library(ggplot2)
  library(GenomicRanges)

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

        ### NEW CODE START: Skip if the eORF doesn't overlap this isoform's exons
        overlap_exons <- findOverlaps(eORF_ranges, exons_gr_original)
        if (length(overlap_exons) == 0) {
          # If no overlap with any exon of this isoform, skip plotting eORF for this isoform
          next
        }
        ### NEW CODE END

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
              height_factor = ifelse(overlaps_CDS, 8/15, 1),
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

  # If plot_ORF_ranges and eORFTxInfo are set, add vertical lines indicating ORF start/end
  if (plot_ORF_ranges && !is.null(eORFTxInfo)) {
    for (eORF_idx in seq_along(eORFTxInfo$eORF.tx_id)) {
      eORF_ranges <- eORFTxInfo$xlim.eORF[[eORF_idx]]
      eORF_left <- eORFTxInfo$eORF_left[eORF_idx]
      eORF_right <- eORFTxInfo$eORF_right[eORF_idx]

      # Adjust start and end based on strand
      if (strand == "+") {
        start_pos <- eORF_left
        end_pos <- eORF_right
      } else {
        start_pos <- eORF_right
        end_pos <- eORF_left
      }

      # Check if eORF overlaps with CDS
      overlaps_CDS <- FALSE
      if (!is.null(GeneTxInfo$xlimCds[[tx_id]]) && length(GeneTxInfo$xlimCds[[tx_id]]) > 0) {
        cds_ranges_check <- GeneTxInfo$xlimCds[[tx_id]]
        overlap_cds_check <- findOverlaps(eORF_ranges, cds_ranges_check)
        if (length(overlap_cds_check) > 0) {
          overlaps_CDS <- TRUE
        }
      }
      line_color <- if (overlaps_CDS) "hotpink" else "orange"

      # Add vertical lines for ORF boundaries
      p_gene <- p_gene +
        geom_vline(xintercept = start_pos, linetype = "solid", color = line_color, alpha = 0.5) +
        geom_vline(xintercept = end_pos, linetype = "dashed", color = line_color, alpha = 0.5)
    }
  }

  return(p_gene)
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
#' @param left_margin Numeric value for the left margin extension in the plot. Default is 1.5.
#' @param right_margin Numeric value for the right margin extension in the plot. Default is 1.5.
#' @param fExtend Numeric value specifying the number of nucleotides to extend the frame assignment into the 5' UTR. Default is 0.
#' @param tExtend Numeric value specifying the number of nucleotides to extend the frame assignment into the 3' UTR. Default is 0.
#' @param RNAseq List of file paths to RNA-seq BAM files.
#' @param Riboseq List of data frames containing Ribo-seq data.
#' @param SampleNames Vector of sample names corresponding to the RNAseq and Riboseq data.
#' @param data_types Vector of sample data type names for each sample (e.g., "Ribo-seq").
#' @param GRangeInfo Genomic range information, typically an object like `Txome_Range`.
#' @param RNAseqBamPaired Vector indicating whether each RNA-seq BAM file is paired-end ("paired") or single-end ("single").
#' @param Y_scale Character string, either "all" or "each", specifying how to scale the Y-axis for RNA-seq coverage. Default is "all".
#' @param Ribo_fix_height Numeric value to fix the maximum height of Ribo-seq counts in the plot.
#' @param plot_ORF_ranges Logical indicating whether to plot ORF ranges in the gene model. Default is FALSE.
#' @param oORF_coloring Character string specifying coloring method for overlapping ORFs ("oORF_colors" or "extend_mORF").
#' @param frame_colors Named vector of colors for the reading frames (0,1,2). Default is c("0"="#FF0000", "1"="#3366FF", "2"="#009900").
#' @param plot_range Optional numeric vector specifying a custom genomic range to plot.
#' @param sample_color Vector specifying colors for each sample or "color" to use default coloring.
#' @param show_seq Logical indicating whether to display the DNA and amino acid sequence. Default is FALSE.
#' @param FASTA A `BSgenome` object or path to FASTA with genomic sequences.
#' @param plot_genomic_direction Logical whether to plot the genomic direction arrow. Default is FALSE. Genomic direction indicates the direction you should see on genome browsers such as JBrowse or IGV.
#' @param dna_aa_height_ratio Numeric value adjusting DNA/AA plot height. Default is 0.5.
#' @param gene_model_height_ratio Numeric value adjusting gene model plot height or NULL to auto-adjust.
#' @param transcript_label_font_size Numeric for transcript label font size.
#' @param selected_isoforms Optional vector of transcript IDs to plot. If provided, only these isoforms plus `tx_id` are shown.
#'
#' @return A combined ggplot object displaying RNA-seq coverage, Ribo-seq data, gene models, and optional sequences.
#' @export
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
#' @param left_margin Numeric value for the left margin extension in the plot. Default is 1.5.
#' @param right_margin Numeric value for the right margin extension in the plot. Default is 1.5.
#' @param fExtend Numeric value specifying the number of nucleotides to extend the frame assignment into the 5' UTR. Default is 0.
#' @param tExtend Numeric value specifying the number of nucleotides to extend the frame assignment into the 3' UTR. Default is 0.
#' @param RNAseq List of file paths to RNA-seq BAM files.
#' @param Riboseq List of data frames containing Ribo-seq data.
#' @param SampleNames Vector of sample names corresponding to the RNAseq and Riboseq data.
#' @param data_types Vector of sample data type names for each sample (e.g., "Ribo-seq").
#' @param GRangeInfo Genomic range information, typically an object like `Txome_Range`.
#' @param RNAseqBamPaired Vector indicating whether each RNA-seq BAM file is paired-end ("paired") or single-end ("single").
#' @param Y_scale Character string, either "all" or "each", specifying how to scale the Y-axis for RNA-seq coverage. Default is "all".
#' @param Ribo_fix_height Numeric value to fix the maximum height of Ribo-seq counts in the plot.
#' @param plot_ORF_ranges Logical indicating whether to plot ORF ranges in the gene model. Default is FALSE.
#' @param oORF_coloring Character string specifying coloring method for overlapping ORFs ("oORF_colors" or "extend_mORF").
#' @param frame_colors Named vector of colors for the reading frames (0,1,2). Default is c("0"="#FF0000", "1"="#3366FF", "2"="#009900").
#' @param plot_range Optional numeric vector specifying a custom genomic range to plot.
#' @param sample_color Vector specifying colors for each sample or "color" to use default coloring.
#' @param show_seq Logical indicating whether to display the DNA and amino acid sequence. Default is FALSE.
#' @param FASTA A `BSgenome` object or path to FASTA with genomic sequences.
#' @param plot_genomic_direction Logical whether to plot the genomic direction arrow. Default is FALSE. Genomic direction indicates the direction you should see on genome browsers such as JBrowse or IGV.
#' @param dna_aa_height_ratio Numeric value adjusting DNA/AA plot height. Default is 0.5.
#' @param gene_model_height_ratio Numeric value adjusting gene model plot height or NULL to auto-adjust.
#' @param transcript_label_font_size Numeric for transcript label font size.
#' @param selected_isoforms Optional vector of transcript IDs to plot. If provided, only these isoforms plus `tx_id` are shown.
#'
#' @return A combined ggplot object displaying RNA-seq coverage, Ribo-seq data, gene models, and optional sequences.
#' @export
ggRibo <- function(gene_id, tx_id, eORF.tx_id = NULL,
                   eORFRangeInfo = NULL, Extend = 100, NAME = "",
                   RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                   fExtend = 0,
                   tExtend = 0,
                   RNAseq = RNAseqData,
                   Riboseq = RiboseqData,
                   SampleNames = Samples,
                   GRangeInfo = Txome_Range,
                   RNAseqBamPaired = RNAseqBamPairorSingle,
                   Y_scale = "all",
                   Ribo_fix_height = NULL,
                   plot_ORF_ranges = FALSE,
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
                   selected_isoforms = NULL) {

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
  if (!is.null(Riboseq)) {
    Riboseq_list <- lapply(seq_along(Riboseq), function(x) {
      Riboseq[[x]][Riboseq[[x]]$chr == chr &
                     Riboseq[[x]]$position >= range_left &
                     Riboseq[[x]]$position <= range_right &
                     Riboseq[[x]]$strand == strand_info, ]
    })
  } else {
    Riboseq_list <- list()
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
  if (!is.null(RNAseq)) {
    what1 <- c("rname","strand","pos","qwidth","seq")
    param <- ScanBamParam(which=GeneTxInfo$generangesplus, what=what1)
  }

  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    if (is.null(RNAseqBamPaired)) {
      stop("RNAseqBamPaired must be provided when RNAseq data is supplied.")
    }

    if (length(RNAseqBamPaired) != length(RNAseq)) {
      stop("Length of RNAseqBamPaired must match the number of RNAseq samples.")
    }

    # Compute coverage over gene range
    global_start <- min(start(GeneTxInfo$generangesplus))
    global_end <- max(end(GeneTxInfo$generangesplus))

    RNAseq_list <- lapply(seq_len(length(RNAseq)), function(i) {
      # Read and coverage calculation depending on paired or single-end
      if (RNAseqBamPaired[i] == "paired") {
        readPairs1 <- suppressWarnings(
          readGAlignmentPairs(RNAseq[i], param=param, strandMode=2)
        )
        if (length(readPairs1)==0) {
          warning(paste("No paired-end reads found in", RNAseq[i]))
          Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
        } else {
          readPairs1 <- readPairs1[strand(readPairs1)==GeneTxInfo$strand]
          if (length(readPairs1)==0) {
            warning(paste("No reads matching strand",GeneTxInfo$strand,"found in",RNAseq[i]))
            Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
          } else {
            cvg1 <- coverage(readPairs1)
            if (is.null(cvg1[[GeneTxInfo$chr]])) {
              Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
            } else {
              Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][global_start:global_end])
            }
          }
        }
        Gtx1
      } else if (RNAseqBamPaired[i] == "single") {
        alignments <- suppressWarnings(readGAlignments(RNAseq[i], param=param))
        if (length(alignments)==0) {
          warning(paste("No single-end reads found in", RNAseq[i]))
          Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
        } else {
          alignments <- alignments[strand(alignments)==GeneTxInfo$strand]
          if (length(alignments)==0) {
            warning(paste("No reads matching strand",GeneTxInfo$strand,"found in",RNAseq[i]))
            Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
          } else {
            cvg1 <- coverage(alignments)
            if (is.null(cvg1[[GeneTxInfo$chr]])) {
              Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
            } else {
              Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][global_start:global_end])
            }
          }
        }
        Gtx1
      } else {
        stop(paste("Invalid value in RNAseqBamPaired for sample", SampleNames[i],
                   "- expected 'paired' or 'single', got", RNAseqBamPaired[i]))
      }
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
      p <- p + geom_step(data=RNAseq_df_line, aes(x=position, y=count), color=RNAcoverline, na.rm=TRUE)

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
            p <- p + geom_segment(data=RiboRslt, aes(x=position, xend=position, y=0, yend=count_scaled, color=frame))
            p <- p + scale_color_manual(values=frame_colors, na.value="grey")
          } else {
            p <- p + geom_segment(data=RiboRslt, aes(x=position, xend=position, y=0, yend=count_scaled), color=sample_color_i)
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
                                    aes(x=position, xend=position, y=0, yend=count_scaled), color='grey')
              p <- p + geom_segment(data=Ribo_main[Ribo_main$region_type=='overlapping',],
                                    aes(x=position, xend=position, y=0, yend=count_scaled, color=frame))

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
                                          aes(x=position,xend=position,y=0,yend=count_scaled,color=frame))
                  }
                }
              }
              p <- p + scale_color_manual(values=frame_colors, na.value='grey')

            } else {
              p <- p + geom_segment(data=Ribo_main,
                                    aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
              if (!is.null(eORFTxInfo)) {
                for (j in seq_along(eORFTxInfo$eORF.tx_id)) {
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]
                  if (nrow(eORF_Riboseq)>0) {
                    if (!is.null(Ribo_fix_height)) {
                      eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                    }
                    eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                    p <- p + geom_segment(data=eORF_Riboseq,
                                          aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
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
                             aes(x=position,xend=position,y=0,yend=count_scaled,color=frame))

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
                                   aes(x=position,xend=position,y=0,yend=count_scaled,color=frame))
                  }
                }
              }

              p <- p + scale_color_manual(values=frame_colors, na.value="grey")

            } else {
              p <- p +
                geom_segment(data=Ribo_main,
                             aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
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
                                   aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
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
                             aes(x=position,xend=position,y=0,yend=count_scaled,color=frame))

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

                  if (length(overlaps_CDS)==0 && overlaps_fiveUTR) {
                    if (nrow(eORF_Riboseq)>0) {
                      if (!is.null(Ribo_fix_height)) {
                        eORF_Riboseq$count <- pmin(eORF_Riboseq$count,Ribo_fix_height)
                      }
                      eORF_Riboseq$count_scaled <- eORF_Riboseq$count * scale_factor_Ribo
                      eORF_Riboseq <- assign_frames(eORF_Riboseq, eORF_ranges, GeneTxInfo$strand)

                      p <- p +
                        geom_segment(data=eORF_Riboseq,
                                     aes(x=position,xend=position,y=0,yend=count_scaled,color=frame))
                    }
                  }
                }
              }

              p <- p + scale_color_manual(values=frame_colors, na.value='grey')

            } else {
              p <- p +
                geom_segment(data=RiboRslt,
                             aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
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
                                   aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
                  }
                }
              }
            }

          } else {
            # Default coding transcripts frame assignment
            Ribo_main <- RiboRslt
            if (!is.null(eORFTxInfo)) {
              eORF_positions <- unlist(lapply(seq_along(eORFTxInfo$xlim.eORF), function(j) {
                ranges <- eORFTxInfo$xlim.eORF[[j]]
                seq(min(start(ranges)), max(end(ranges)))
              }))
              Ribo_main <- Ribo_main[!(Ribo_main$position %in% eORF_positions), ]
            }

            Ribo_main <- assign_frames(Ribo_main, cds_ranges, GeneTxInfo$strand)

            if (!is.null(Ribo_fix_height)) {
              Ribo_main$count <- pmin(Ribo_main$count,Ribo_fix_height)
              Ribo_main$count_scaled <- Ribo_main$count * scale_factor_Ribo
            }

            if (sample_color_i=="color") {
              p <- p +
                geom_segment(data=Ribo_main,
                             aes(x=position,xend=position,y=0,yend=count_scaled,color=frame))

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
                                   aes(x=position,xend=position,y=0,yend=count_scaled,color=frame))
                  }
                }
              }

              p <- p + scale_color_manual(values=frame_colors, na.value="grey")

            } else {
              p <- p +
                geom_segment(data=Ribo_main,
                             aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
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
                                   aes(x=position,xend=position,y=0,yend=count_scaled),color=sample_color_i)
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
      y_label <- if (!is.null(Ribo_fix_height) || Y_scale=="all") {
        current_max_Y
      } else {
        current_max_Y + 0.01*current_max_Y
      }

      p <- p + annotate("text",
                        x=x_label,y=y_label,
                        label=SampleNames[i],
                        hjust=0,vjust=0,
                        size=3,fontface="bold")

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
      FASTA=FASTA
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
  rna_ribo_height <-0.8

  if (is.null(gene_model_height_ratio)) {
    gene_model_height_ratio <-0.2+(num_transcripts)*0.1
  }
  gene_model_height <- gene_model_height_ratio

  if (show_seq && !is.null(FASTA)) {
    dna_aa_height <- dna_aa_height_ratio
  } else {
    dna_aa_height <-0
  }

  total_height_units <- title_height+(num_datasets*rna_ribo_height)+dna_aa_height+gene_model_height
  rel_heights <- c(
    title_height,
    rep(rna_ribo_height,num_datasets),
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
    plotlist=c(plot_list,list(dna_aa_plot),list(gene_model_plot)),
    ncol=1,
    align="v",
    rel_heights=rel_heights,
    axis="lr",
    labels=NULL,
    label_size=10,
    label_fontface="plain"
  )
  return(combined_plot)
}

# -----------------------------------
# ggRibo_decom function
# -----------------------------------
#'
#' `ggRibo_decom` creates a combined visualization of RNA-Seq coverage and frame-specific Ribo-Seq counts for a specified gene and transcript.
#' It generates three separate plots corresponding to the three reading frames (0, 1, and 2) of Ribo-Seq data, optionally including reads
#' that do not fall within any annotated ORF regions. It can also display extended ORFs (eORFs), genomic sequences, and gene models.
#'
#' @param gene_id Character. The identifier for the gene of interest.
#' @param tx_id Character. The transcript identifier within the gene for which the main ORF is annotated.
#' @param eORF.tx_id Character vector, optional. Transcript identifiers for extended ORFs (eORFs) associated with the gene.
#' @param eORFRangeInfo List, optional. Contains eORF range information. Required if `eORF.tx_id` is specified.
#' @param Extend Numeric or numeric vector of length 2. Extends the plotting range upstream and downstream of the gene. Defaults to 100.
#' @param NAME Character. An optional label or title.
#' @param RNAcoverline Character. Color for the RNA-Seq coverage line. Defaults to "grey".
#' @param RNAbackground Character or vector of length equal to number of samples. Fill color for RNA-Seq coverage bars. Defaults to "#FEFEAE".
#' @param fExtend Numeric. Extends the ORF upstream by this many bases. Defaults to 0.
#' @param tExtend Numeric. Extends the ORF downstream by this many bases. Defaults to 0.
#' @param RNAseq Character vector. Paths to RNA-Seq BAM files.
#' @param Riboseq Character vector. Paths to Ribo-Seq data files (as data frames).
#' @param SampleNames Character vector. Names corresponding to the samples.
#' @param GRangeInfo List. Genomic range information such as transcripts, exons, CDS, etc.
#' @param RNAseqBamPaired Character vector. Indicates if each RNA-Seq BAM is paired-end ("paired") or single-end ("single").
#' @param Y_scale Character. Either "all" or "each", controlling the Y-axis scaling for RNA-seq coverage. Defaults to "all".
#' @param Ribo_fix_height Numeric, optional. Caps Ribo-Seq counts at a fixed height, ignoring Y_scale.
#' @param plot_ORF_ranges Logical. If TRUE, highlights annotated ORF ranges on the gene model. Defaults to FALSE.
#' @param oORF_coloring Character, optional. Method for coloring overlapping ORFs. "oORF_colors" or "extend_mORF".
#' @param frame_colors Named character vector. Colors for frames 0, 1, and 2. Defaults provided.
#' @param plot_range Numeric vector of length 2, optional. Custom genomic range to plot.
#' @param sample_color Character. If "color", uses frame-specific colors. Otherwise, uses a single color for Ribo reads.
#' @param show_seq Logical. If TRUE, displays DNA and AA sequences below the coverage plots. Defaults to FALSE.
#' @param FASTA Optional. Path to a FASTA file or `BSgenome` object with genomic sequences.
#' @param dna_aa_height_ratio Numeric. Adjusts DNA/AA plot height relative to gene model height. Defaults to 0.5.
#' @param gene_model_height_ratio Numeric, optional. Adjusts gene model plot height. If NULL, auto-scales.
#' @param transcript_label_font_size Numeric. Font size for transcript ID labels in gene model. Defaults to 10.
#' @param plot_genomic_direction Logical. If TRUE, draws an arrow indicating genomic direction on the top plot. Defaults to FALSE.
#' @param data_types Character vector. Describes data type(s) for samples (e.g., "Ribo-seq"). Must match SampleNames length.
#' @param plot_unassigned_reads Logical. If TRUE, plots Ribo-Seq reads not assigned to any ORF as grey segments.
#' @param selected_isoforms Optional. Vector of transcript IDs to plot. If specified, only these isoforms and tx_id are shown.
#'
#' @return A combined `ggplot` object with RNA-Seq coverage, three frame-specific Ribo-Seq plots, gene model, and optionally DNA/AA sequences.
#'
#' @export
ggRibo_decom <- function(gene_id, tx_id, eORF.tx_id = NULL,
                         eORFRangeInfo = NULL, Extend = 100, NAME = "",
                         RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                         fExtend = 0,
                         tExtend = 0,
                         RNAseq = RNAseqData,
                         Riboseq = RiboseqData,
                         SampleNames = Samples,
                         GRangeInfo = Txome_Range,
                         RNAseqBamPaired = RNAseqBamPairorSingle,
                         Y_scale = "all",
                         Ribo_fix_height = NULL,
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
                         selected_isoforms = NULL
) {
  # Validate data_types length matches SampleNames
  if (length(data_types) != length(SampleNames)) {
    stop("The length of data_types must match the number of samples.")
  }

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

  # ggRibo_decom only supports one sample at a time
  if (length(SampleNames) > 1) {
    stop("ggRibo_decom only supports one sample at a time.")
  }

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
    stop("No transcripts left after applying selected_isoforms in ggRibo().")
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

  # Filter Riboseq data to the region
  if (!is.null(Riboseq)) {
    Riboseq_list <- list(Riboseq[[1]][Riboseq[[1]]$chr == chr &
                                        Riboseq[[1]]$position >= range_left &
                                        Riboseq[[1]]$position <= range_right &
                                        Riboseq[[1]]$strand == strand_info, ])
  } else {
    Riboseq_list <- list()
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

  # Process RNAseq coverage
  if (!is.null(RNAseq)) {
    what1 <- c("rname","strand","pos","qwidth","seq")
    param <- ScanBamParam(which=GeneTxInfo$generangesplus, what=what1)
  }

  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    if (is.null(RNAseqBamPaired)) {
      stop("RNAseqBamPaired must be provided.")
    }
    if (length(RNAseqBamPaired) != length(RNAseq)) {
      stop("Length of RNAseqBamPaired must match number of RNAseq samples.")
    }

    global_start <- min(start(GeneTxInfo$generangesplus))
    global_end <- max(end(GeneTxInfo$generangesplus))

    i <- 1
    # Only one sample supported here
    if (RNAseqBamPaired[i] == "paired") {
      readPairs1 <- suppressWarnings(
        readGAlignmentPairs(RNAseq[i], param=param, strandMode=2)
      )
      if (length(readPairs1)==0) {
        Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
      } else {
        readPairs1 <- readPairs1[strand(readPairs1)==GeneTxInfo$strand]
        if (length(readPairs1)==0) {
          Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
        } else {
          cvg1 <- coverage(readPairs1)
          if (is.null(cvg1[[GeneTxInfo$chr]])) {
            Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
          } else {
            Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][global_start:global_end])
          }
        }
      }
      Gtx1
    } else if (RNAseqBamPaired[i] == "single") {
      alignments <- suppressWarnings(readGAlignments(RNAseq[i], param=param))
      if (length(alignments)==0) {
        Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
      } else {
        alignments <- alignments[strand(alignments)==GeneTxInfo$strand]
        if (length(alignments)==0) {
          Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
        } else {
          cvg1 <- coverage(alignments)
          if (is.null(cvg1[[GeneTxInfo$chr]])) {
            Gtx1 <- numeric(length=width(GeneTxInfo$generangesplus))
          } else {
            Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][global_start:global_end])
          }
        }
      }
      Gtx1
    } else {
      stop(paste("Invalid value in RNAseqBamPaired:", RNAseqBamPaired[i]))
    }

    RNAseq_list[[1]] <- Gtx1
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
        geom_step(data=RNAseq_df_line, aes(x=position, y=count), color=RNAcoverline, na.rm=TRUE) +
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
        p <- p + geom_segment(data=na_data, aes(x=position, xend=position, y=0, yend=count_scaled), color="grey", na.rm=TRUE)
      }

      if (nrow(Ribo_df)>0) {
        p <- p + geom_segment(data=Ribo_df, aes(x=position, xend=position, y=0, yend=count_scaled), color=frame_color, na.rm=TRUE)
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

    # Assign frames depending on whether transcript is coding or not
    if (length(cds_ranges)==0) {
      # Noncoding: frames from start of transcript
      if (GeneTxInfo$strand=="+") {
        exons_sorted <- sort(exons,decreasing=FALSE)
      } else {
        exons_sorted <- sort(exons,decreasing=TRUE)
      }
      positions_all <- integer(0)
      tx_positions <- integer(0)
      cum_len <-0
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
    } else {
      # Coding transcript logic (assign frames)
      # If extend_mORF or oORF_colors are set, or we want to exclude eORFs, this is handled similarly to ggRibo but simplified here
      # For ggRibo_decom, we just do a standard assignment and then separate by frame
      # Below we do a similar approach: assign frames, then split by frame
      if (!is.null(oORF_coloring) && oORF_coloring=="extend_mORF") {
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
      } else if (fExtend>0 || tExtend>0) {
        Ribo_main <- RiboRslt
        if (!is.null(eORFTxInfo)) {
          Ribo_main <- exclude_eORF_reads(Ribo_main, eORFTxInfo, GeneTxInfo$strand)
        }
        RiboRslt <- assign_frames_with_extension(Ribo_main, cds_ranges, exons, fExtend, tExtend, GeneTxInfo$strand)
      } else {
        Ribo_main <- RiboRslt
        if (!is.null(eORFTxInfo)) {
          eORF_positions <- unlist(lapply(seq_along(eORFTxInfo$xlim.eORF), function(j) {
            ranges <- eORFTxInfo$xlim.eORF[[j]]
            seq(min(start(ranges)), max(end(ranges)))
          }))
          Ribo_main <- Ribo_main[!(Ribo_main$position %in% eORF_positions),]
        }
        Ribo_main <- assign_frames(Ribo_main, cds_ranges, GeneTxInfo$strand)
        RiboRslt <- Ribo_main
      }
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
    y_label <- current_max_Y + 0.01*current_max_Y

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

    total_height_units <- title_height+(3*frame_plot_height)+dna_aa_height+gene_model_height
    rel_heights <- c(
      title_height,
      rep(frame_plot_height,3),
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

# -----------------------------------
# ggRNA function
# -----------------------------------
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
#' @param RNAseq A list of file paths to RNA-seq BAM files.
#' @param SampleNames Vector of sample names corresponding to the RNAseq data.
#' @param GRangeInfo A genomic range information object (e.g., `Txome_Range`) containing annotations.
#' @param RNAseqBamPaired A vector indicating whether each RNA-seq BAM file is paired-end ("paired") or single-end ("single").
#' @param Y_scale Character string, either "all" or "each", specifying how to scale the Y-axis for RNA-seq coverage. Default is "all".
#' @param plot_ORF_ranges Logical indicating whether to plot ORF ranges in the gene model. Default is FALSE.
#' @param plot_range Optional numeric vector of length two specifying a custom genomic range to plot.
#' @param show_seq Logical indicating whether to display the DNA and amino acid sequences. Default is FALSE.
#' @param FASTA A `BSgenome` object or path to a FASTA file containing genomic sequences (required if show_seq=TRUE).
#' @param plot_genomic_direction Logical indicating whether to plot an arrow showing genomic direction. Default is FALSE.
#' @param dna_aa_height_ratio Numeric value to adjust the height of the DNA and amino acid sequence plot. Default is 0.5.
#' @param gene_model_height_ratio Numeric value to adjust the height of the gene model plot. If NULL, it is auto-calculated.
#' @param transcript_label_font_size Numeric controlling the font size of the transcript ID labels in the gene model plot.
#' @param selected_isoforms Optional vector of transcript IDs to plot. If provided, only these isoforms (and `tx_id`) will be shown.
#'
#' @return A combined ggplot object displaying RNA-seq coverage, gene models, and optionally genomic sequences.
#' @export
ggRNA <- function(gene_id, tx_id, Extend = 100, NAME = "",
                  RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                  RNAseq = RNAseqData,
                  SampleNames = Samples,
                  GRangeInfo = Txome_Range,
                  RNAseqBamPaired = RNAseqBamPairorSingle,
                  Y_scale = "all",
                  plot_ORF_ranges = FALSE,
                  plot_range = NULL,
                  show_seq = FALSE,
                  FASTA = NULL,
                  dna_aa_height_ratio = 0.5,
                  gene_model_height_ratio = NULL,
                  transcript_label_font_size = 10,
                  plot_genomic_direction = FALSE,
                  selected_isoforms = NULL
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
  if (!is.null(RNAseq)) {
    what1 <- c("rname","strand","pos","qwidth","seq")
    param <- ScanBamParam(which=gene_ranges, what=what1)
  }

  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    # Check RNAseqBamPaired parameter
    if (is.null(RNAseqBamPaired)) {
      stop("RNAseqBamPaired must be provided when RNAseq data is supplied.")
    }

    if (length(RNAseqBamPaired) != length(RNAseq)) {
      stop("Length of RNAseqBamPaired must match the number of RNAseq samples.")
    }

    global_start <- min(start(gene_ranges))
    global_end <- max(end(gene_ranges))

    # Compute coverage for each RNA-seq sample
    RNAseq_list <- lapply(seq_len(length(RNAseq)), function(i) {
      if (RNAseqBamPaired[i] == "paired") {
        # Paired-end RNA-seq
        readPairs1 <- suppressWarnings(readGAlignmentPairs(RNAseq[i], param=param, strandMode=2))
        if (length(readPairs1) == 0) {
          warning(paste("No paired-end reads found in", RNAseq[i]))
          Gtx1 <- numeric(length=width(gene_ranges))
        } else {
          readPairs1 <- readPairs1[strand(readPairs1)==strand_info]
          if (length(readPairs1)==0) {
            warning(paste("No reads matching strand",strand_info,"found in",RNAseq[i]))
            Gtx1 <- numeric(length=width(gene_ranges))
          } else {
            cvg1 <- coverage(readPairs1)
            if (is.null(cvg1[[chr]])) {
              Gtx1 <- numeric(length=width(gene_ranges))
            } else {
              Gtx1 <- as.numeric(cvg1[[chr]][global_start:global_end])
            }
          }
        }
        Gtx1
      } else if (RNAseqBamPaired[i] == "single") {
        # Single-end RNA-seq
        alignments <- suppressWarnings(readGAlignments(RNAseq[i], param=param))
        if (length(alignments)==0) {
          warning(paste("No single-end reads found in", RNAseq[i]))
          Gtx1 <- numeric(length=width(gene_ranges))
        } else {
          alignments <- alignments[strand(alignments)==strand_info]
          if (length(alignments)==0) {
            warning(paste("No reads matching strand",strand_info,"found in",RNAseq[i]))
            Gtx1 <- numeric(length=width(gene_ranges))
          } else {
            cvg1 <- coverage(alignments)
            if (is.null(cvg1[[chr]])) {
              Gtx1 <- numeric(length=width(gene_ranges))
            } else {
              Gtx1 <- as.numeric(cvg1[[chr]][global_start:global_end])
            }
          }
        }
        Gtx1
      } else {
        stop(paste("Invalid value in RNAseqBamPaired for sample", SampleNames[i]))
      }
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
        geom_step(data=RNAseq_df, aes(x=position, y=count), color=RNAcoverline, na.rm=TRUE) +
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
      y_label <- current_max_Y + 0.01*current_max_Y
      p <- p + annotate("text",
                        x=x_label,y=y_label,
                        label=SampleNames[i],
                        hjust=hjust_label,vjust=0,
                        size=3,fontface="bold")

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
      FASTA=FASTA
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

  total_height_units <- title_height+(num_datasets*rna_ribo_height)+dna_aa_height+gene_model_height
  rel_heights <- c(
    title_height,
    rep(rna_ribo_height,num_datasets),
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
    plotlist=c(plot_list,list(dna_aa_plot),list(gene_model_plot)),
    ncol=1,
    align="v",
    rel_heights=rel_heights,
    axis="lr",
    labels=NULL,
    label_size=10,
    label_fontface="plain"
  )

  return(combined_plot)
}
