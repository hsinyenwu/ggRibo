# ---------------------------
# Function Definitions
# ---------------------------

# Function to import GTF/GFF annotation and create Range_info object
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
  # Create a TxDb object from the annotation file
  txdb <- txdbmaker::makeTxDbFromGFF(annotation, format = format, dataSource = dataSource, organism = organism)

  # Extract genomic ranges by transcript
  exonsByTx <- exonsBy(txdb, by = 'tx', use.names = TRUE)
  txByGene <- transcriptsBy(txdb, by = 'gene')
  cdsByTx <- cdsBy(txdb, by = "tx", use.names = TRUE)
  fiveUTR <- fiveUTRsByTranscript(txdb, use.names = TRUE)
  threeUTR <- threeUTRsByTranscript(txdb, use.names = TRUE)

  # Create a Range_info object with the extracted ranges
  Txome_Range <- Range_info$new(
    exonsByTx = exonsByTx,
    txByGene = txByGene,
    cdsByTx = cdsByTx,
    fiveUTR = fiveUTR,
    threeUTR = threeUTR
  )

  # Assign the Range_info object to the global environment for accessibility
  assign("Txome_Range", Txome_Range, envir = .GlobalEnv)
}

# Function to import eORF annotation and create eORF_Range_info object
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

  # Extract CDS ranges by transcript (representing eORFs)
  cdsByTx <- cdsBy(txdb, by = "tx", use.names = TRUE)

  # Create an eORF_Range_info object with the extracted eORF ranges
  eORF_Range <- eORF_Range_info$new(
    eORFByTx = cdsByTx
  )

  # Assign the eORF_Range_info object to the global environment for accessibility
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
  # Read each Ribo-seq data file and store in a list
  Ribo_data_list <- lapply(RiboseqData, function(file) {
    # Read the file as a tab-delimited file without headers
    Ribo1 <- read.delim(file = file, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    # Assign column names for clarity
    colnames(Ribo1) <- c("count", "chr", "position", "strand")
    Ribo1
  })
  # Assign sample names to the list elements for easy reference
  names(Ribo_data_list) <- SampleNames
  return(Ribo_data_list)
}

# Helper Function: Assign Frames with Extended CDS Ranges
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
  if (length(extended_cds_ranges) > 0) {
    # Order exons based on strand and genomic positions
    if (strand == "+") {
      exons <- sort(extended_cds_ranges, decreasing = FALSE)
    } else {
      exons <- sort(extended_cds_ranges, decreasing = TRUE)
    }

    # Initialize vectors to store positions and transcript positions
    positions <- integer(0)
    tx_positions <- integer(0)
    cum_len <- 0

    # Iterate over exons to map genomic positions to transcript positions
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

    # Build data frame mapping positions to transcript positions
    position_df <- data.frame(position = positions, tx_pos = tx_positions)

    # Find the transcript position corresponding to the main ORF start
    if (strand == "+") {
      main_orf_start <- min(start(main_cds_ranges))
    } else {
      main_orf_start <- max(end(main_cds_ranges))
    }
    main_orf_tx_pos <- position_df$tx_pos[position_df$position == main_orf_start][1]

    # Calculate frames relative to the main ORF start
    position_df$frame <- (position_df$tx_pos - main_orf_tx_pos) %% 3
    position_df$frame <- factor(position_df$frame, levels = c(0,1,2))

    # Merge frame information with Ribo-seq data based on position
    Ribo_data <- merge(Ribo_data, position_df[, c("position", "frame")], by = "position", all.x = TRUE)
  } else {
    # If no extended CDS ranges, assign NA to frame
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
  if (length(ranges) > 0) {
    # Order exons based on strand and genomic positions
    if (strand == "+") {
      exons <- sort(ranges, decreasing = FALSE)
    } else {
      exons <- sort(ranges, decreasing = TRUE)
    }

    # Initialize vectors to store positions and transcript positions
    positions <- integer(0)
    tx_positions <- integer(0)
    cum_len <- 0

    # Iterate over exons to map genomic positions to transcript positions
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

    # Build data frame mapping positions to transcript positions
    position_df <- data.frame(position = positions, tx_pos = tx_positions)

    # Calculate frames based on transcript positions
    position_df$frame <- (position_df$tx_pos - 1) %% 3
    position_df$frame <- factor(position_df$frame, levels = c(0,1,2))

    # Merge frame information with Ribo-seq data based on position
    Ribo_data <- merge(Ribo_data, position_df[, c("position", "frame")], by = "position", all.x = TRUE)
  } else {
    # If no ranges provided, assign NA to frame
    Ribo_data$frame <- factor(NA, levels = c(0,1,2))
  }
  return(Ribo_data)
}

# Corrected helper function to assign frames with extension into UTRs
#' @title Assign Frames with Extension into UTRs
#' @description
#' Assigns reading frames to Ribo-seq reads with extension into UTRs.
#'
#' @param Ribo_data Data frame containing Ribo-seq reads.
#' @param cds_ranges CDS ranges as a \code{GRanges} object.
#' @param exons Exon ranges as a \code{GRanges} object.
#' @param fExtend Number of nucleotides to extend into the 5' UTR.
#' @param tExtend Number of nucleotides to extend into the 3' UTR.
#' @param strand Strand information ("+" or "-").
#' @return Data frame with an added \code{frame} column indicating the reading frame.
assign_frames_with_extension <- function(Ribo_data, cds_ranges, exons, fExtend, tExtend, strand) {
  if (length(exons) > 0 && length(cds_ranges) > 0) {
    # Order exons based on strand and genomic positions
    if (strand == "+") {
      exons <- sort(exons, decreasing = FALSE)
    } else {
      exons <- sort(exons, decreasing = TRUE)
    }

    # Build transcript coordinate mapping
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

    # Map CDS ranges to transcript positions
    cds_positions <- unlist(lapply(seq_along(cds_ranges), function(idx) {
      seq(start(cds_ranges[idx]), end(cds_ranges[idx]))
    }))
    if (strand == "-") {
      cds_positions <- rev(cds_positions)
    }
    cds_tx_pos <- position_df$tx_pos[position_df$position %in% cds_positions]
    cds_start_tx <- min(cds_tx_pos)
    cds_end_tx <- max(cds_tx_pos)

    # Extend the CDS in transcript coordinates
    extended_start_tx <- max(1, cds_start_tx - fExtend)
    extended_end_tx <- min(max(position_df$tx_pos), cds_end_tx + tExtend)
    extended_tx_pos <- seq(extended_start_tx, extended_end_tx)

    # Map extended transcript positions back to genomic positions
    extended_positions <- position_df$position[position_df$tx_pos %in% extended_tx_pos]
    frames <- ((extended_tx_pos - cds_start_tx) %% 3)
    frames <- factor(frames, levels = c(0,1,2))

    # Build data frame mapping positions to frames
    frame_df <- data.frame(tx_pos = extended_tx_pos, position = extended_positions, frame = frames)

    # Remove duplicates to ensure unique position-frame mapping
    frame_df <- frame_df[!duplicated(frame_df$position), ]

    # Merge frame information with Ribo-seq data
    Ribo_data <- merge(Ribo_data, frame_df[, c("position", "frame")], by = "position", all.x = TRUE)
  } else {
    # If no exons or CDS ranges provided, assign NA to frame
    Ribo_data$frame <- factor(NA, levels = c(0,1,2))
  }
  return(Ribo_data)
}

# Function to exclude Ribo-seq reads that overlap eORF regions
# Helper function to exclude eORF reads from main Ribo-seq data
#' @title Exclude eORF Reads from Main Ribo-seq Data
#' @description
#' Excludes Ribo-seq reads that overlap with eORF regions from the main Ribo-seq data.
#'
#' @param Ribo_data Data frame containing Ribo-seq reads.
#' @param eORFTxInfo An \code{eORF_info} object containing eORF information.
#' @param strand Strand information ("+" or "-").
#' @return Filtered Ribo-seq data frame.
exclude_eORF_reads <- function(Ribo_data, eORFTxInfo, strand) {
  if (is.null(eORFTxInfo) || length(eORFTxInfo$xlim.eORF) == 0) {
    # If no eORF information provided, return Ribo_data as is
    return(Ribo_data)
  }

  # Create a GRanges object for Ribo_data based on chromosome, position, and strand
  Ribo_gr <- GRanges(
    seqnames = Ribo_data$chr,
    ranges = IRanges(Ribo_data$position, Ribo_data$position),
    strand = Ribo_data$strand
  )

  # Create a GRangesList of eORF regions
  eORF_grl <- GRangesList(eORFTxInfo$xlim.eORF)

  # Find overlaps between Ribo-seq reads and eORF regions
  overlaps <- findOverlaps(Ribo_gr, unlist(eORF_grl))

  # Exclude reads that overlap eORF regions
  if (length(overlaps) > 0) {
    Ribo_data <- Ribo_data[-queryHits(overlaps), ]
  }

  return(Ribo_data)
}

# Function to plot DNA and Amino Acid sequences
#' @title Plot DNA and Amino Acid Sequences
#' @description
#' Generates a plot showing DNA nucleotides and their corresponding amino acids.
#'
#' @param GeneTxInfo A \code{Gene_info} object containing gene-specific information.
#' @param plot_range Optional vector specifying the genomic range to plot.
#' @param FASTA A \code{BSgenome} object containing the reference genome sequences.
#' @return A \code{ggplot2} object representing the DNA and amino acid sequences.
plotDNAandAA <- function(GeneTxInfo, plot_range = NULL, FASTA = NULL) {
  # Extract necessary information from GeneTxInfo
  genelim <- c(GeneTxInfo$range_left, GeneTxInfo$range_right)
  strand <- GeneTxInfo$strand

  # Adjust genelim to be within chromosome bounds
  chrom_length <- seqlengths(FASTA)[GeneTxInfo$chr]
  genelim_adj <- pmax(pmin(genelim, chrom_length), 1)

  # If plot_range is provided, sort and adjust genelim_adj accordingly
  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)  # Sort plot_range from small to large
    genelim_adj <- pmax(pmin(plot_range, chrom_length), 1)
  }

  # Determine the range length
  range_length <- abs(diff(range(genelim_adj))) + 1  # +1 to include both ends

  # Initialize suppress_labels
  suppress_labels <- FALSE

  # Check if we need to suppress labels or print warning
  if (range_length > 201) {
    suppress_labels <- TRUE
  }

  # Determine if range is longer than 201 nt
  long_range <- range_length > 201

  # Extract DNA sequence
  seq_region <- GRanges(
    seqnames = GeneTxInfo$chr,
    ranges = IRanges(min(genelim_adj), max(genelim_adj)),
    strand = "+"
  )
  seqs <- getSeq(FASTA, seq_region)
  dna_seq <- as.character(seqs)

  # Reverse complement if strand is "-"
  if (GeneTxInfo$strand == "-") {
    dna_seq <- as.character(reverseComplement(DNAString(dna_seq)))
  }

  # Split sequence into individual nucleotides
  dna_chars <- unlist(strsplit(dna_seq, split = ""))

  # Adjust positions_seq based on strand
  positions_seq <- if (GeneTxInfo$strand == "+") {
    seq(min(genelim_adj), max(genelim_adj))
  } else {
    seq(max(genelim_adj), min(genelim_adj), by = -1)
  }

  # Create data frame for DNA sequence
  dna_df <- data.frame(position = positions_seq,
                       nucleotide = dna_chars,
                       stringsAsFactors = FALSE,
                       row.names = NULL)

  # Assign colors based on nucleotides
  nucleotide_colors <- c(
    "A" = "#00FF00",   # Green
    "T" = "#FF0200",   # Red
    "C" = "#4747FF",   # Blue
    "G" = "#FFA503",   # Orange
    "N" = "grey"
  )

  # Add 'fill_value' column for nucleotides
  dna_df$fill_value <- dna_df$nucleotide

  # Calculate font size based on the number of nucleotides and plotting range
  num_nucleotides <- length(dna_chars)
  plot_width <- abs(diff(range(genelim_adj)))
  # Adjust font size; you may need to fine-tune the scaling factor
  font_size <- (plot_width / num_nucleotides) * 1.2  # Scaling factor of 2 for visibility

  # Ensure font size is within reasonable limits
  font_size <- max(min(font_size, 5), 2)  # Limit between 2 and 5

  # Increase font size by 1.5 fold
  font_size <- font_size * 1.5

  # Generate 3-frame translations
  aa_sequences <- list()
  frames <- c(0, 1, 2)
  # Define the colors as per your request
  frame_colors <- c("Annotated" = "#F1F1F1", "+1" = "#E6E6E6", "+2" = "#C9C9C9")

  # Compute annotated frame based on CDS start position
  cds_ranges <- GeneTxInfo$xlimCds[[GeneTxInfo$tx_id]]  # CDS ranges for main transcript
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
    # No CDS information, default to frame 1
    annotated_frame <- 1
  }
  annotated_frame_zero_based <- as.integer(annotated_frame - 1)  # Convert to 0-based index and ensure it's integer

  # Ensure annotated_frame_zero_based is within 0 to 2
  if (is.na(annotated_frame_zero_based) || annotated_frame_zero_based < 0 || annotated_frame_zero_based > 2) {
    annotated_frame_zero_based <- 0
  }

  # Create frame order: annotated frame first, then +1, then +2
  other_frames <- frames[frames != annotated_frame_zero_based]
  frame_order <- c(annotated_frame_zero_based, sort(other_frames))

  # Adjust y positions based on whether DNA sequence is plotted
  if (long_range) {
    # DNA sequence is not plotted, adjust amino acid positions
    frame_y_positions <- c(0.8, 0.6, 0.4)
  } else {
    frame_y_positions <- c(0.6, 0.4, 0.2)
  }
  frame_y_map <- setNames(frame_y_positions, frame_order)

  # Assign frame labels for fill_value
  frame_labels <- c("Annotated", "+1", "+2")
  frame_label_map <- setNames(frame_labels, frame_order)

  # Assign frame colors as per your request
  frame_color_map <- setNames(frame_colors, frame_labels)

  for (frame in frames) {
    # Determine codon starts and middles
    codon_starts <- seq(frame + 1, nchar(dna_seq) - 2, by = 3)
    codon_middles <- codon_starts + 1
    # Check if there are any codons in this frame
    if(length(codon_starts) == 0) next  # No complete codons in this frame

    # Create DNAString for complete codons only
    dna_subseq <- substr(dna_seq, codon_starts[1], codon_starts[length(codon_starts)] + 2)
    dna_string <- DNAString(dna_subseq)

    # Translate to amino acid sequence with suppressed warnings
    aa_seq <- suppressWarnings(as.character(translate(dna_string, if.fuzzy.codon = "X")))

    # Split amino acids into characters
    aa_chars <- unlist(strsplit(aa_seq, split = ""))

    # Positions of amino acids (middle nucleotide positions)
    aa_positions <- positions_seq[codon_middles]

    # Create data frame for amino acids
    aa_df <- data.frame(
      position = aa_positions,
      amino_acid = aa_chars,
      y = frame_y_map[as.character(frame)],
      frame = frame,  # Frame number (0-based)
      stringsAsFactors = FALSE,
      row.names = NULL
    )

    aa_sequences[[frame + 1]] <- aa_df
  }

  # Combine amino acid data frames
  aa_df_combined <- do.call(rbind, aa_sequences)

  # Set font size for amino acids (same as nucleotides)
  aa_font_size <- font_size*1.25

  # Assign fill_value based on amino acid and frame
  aa_df_combined$frame_label <- frame_label_map[as.character(aa_df_combined$frame)]
  aa_df_combined$fill_value <- aa_df_combined$frame_label

  # Override fill_value for 'M' and '*' amino acids
  aa_df_combined$fill_value[aa_df_combined$amino_acid == 'M'] <- 'Start'
  aa_df_combined$fill_value[aa_df_combined$amino_acid == '*'] <- 'Stop'

  # Create data frame for labeled amino acids (only 'M' and '*')
  aa_label_df <- aa_df_combined[aa_df_combined$amino_acid %in% c('M', '*'), ]

  # Define fill colors
  fill_colors <- c(
    nucleotide_colors,
    "Start" = "green",
    "Stop" = "red",
    frame_color_map
  )

  # Determine which breaks are present in the data
  available_breaks <- intersect(unique(aa_df_combined$fill_value), c("Start", "Stop"))

  # Split amino acid data into regular vs. start/stop for separate plotting
  aa_start_stop_df <- aa_df_combined[aa_df_combined$fill_value %in% c("Start","Stop"),]
  aa_regular_df <- aa_df_combined[!aa_df_combined$fill_value %in% c("Start","Stop"),]

  # Create the DNA and amino acid plot
  p_dna_aa <- ggplot()

  # Plot DNA sequence as colored tiles if long_range is FALSE
  if (!long_range) {
    p_dna_aa <- p_dna_aa +
      suppressWarnings(
        geom_tile(
          data = dna_df,
          aes(x = position, y = 0.8, fill = fill_value),
          width = 1,
          height = 0.2,
          color = "darkgrey",   # Thin dark grey border
          linewidth = 0.2,
          show.legend = FALSE,
          na.rm = TRUE
        )
      )

    # Only add nucleotide letters if suppress_labels is FALSE
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

  # Plot amino acid sequences (regular ones first)
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
        linewidth = 0.5, # Increased line width for start/stop codons
        show.legend = TRUE,
        na.rm = TRUE
      )
    )

  # Only add amino acid letters if suppress_labels is FALSE
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

  # Define fill scale with name set to NULL to remove legend title
  # Override the legend aesthetics to remove grey line around the legend keys
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

  # Adjust x-axis scale based on strand direction
  if (GeneTxInfo$strand == "-") {
    p_dna_aa <- p_dna_aa + scale_x_reverse(limits = c(max(genelim_adj), min(genelim_adj)))
  } else {
    p_dna_aa <- p_dna_aa + scale_x_continuous(limits = c(min(genelim_adj), max(genelim_adj)))
  }

  # Adjust y-axis limits and margins to prevent plot cutoff
  if (long_range) {
    y_limits <- c(0.2, 1)  # Adjust y-axis limits when DNA sequence is not plotted
  } else {
    y_limits <- c(0, 1)
  }
  p_dna_aa <- p_dna_aa +
    scale_y_continuous(limits = y_limits) +
    theme_void() +
    theme(
      plot.margin = unit(c(-1, 0.2, -0.5, 0.2), "lines")
    )

  # Return the DNA and amino acid plot
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

  isoforms <- GeneTxInfo$num_isoforms
  genelim <- c(GeneTxInfo$range_left, GeneTxInfo$range_right)
  tx_names <- GeneTxInfo$tx_names
  tx_id <- GeneTxInfo$tx_id
  strand <- GeneTxInfo$strand

  plot_data_list <- list()
  line_data_list <- list()
  idx <- 1

  other_tx_names <- setdiff(tx_names, tx_id)
  sorted_tx_names <- c(tx_id, sort(other_tx_names))

  y_step <- 0.3
  y_positions <- seq(1, by=y_step, length.out=length(sorted_tx_names))
  isoform_positions <- data.frame(
    isoform = sorted_tx_names,
    y = rev(y_positions),
    stringsAsFactors = FALSE
  )

  isoform_y_map <- setNames(isoform_positions$y, isoform_positions$isoform)

  for (isoform in isoform_positions$isoform) {
    y_value <- isoform_y_map[isoform]
    isoform_data_list <- list()
    isoform_idx <- 1

    exons_gr <- GeneTxInfo$exonByYFGtx[[isoform]]
    if (length(exons_gr)==0) {
      warning(paste("Exons for isoform", isoform, "not found in exonByYFGtx"))
      next
    }

    if (!is.null(plot_range)) {
      segment_gr <- GRanges(seqnames=GeneTxInfo$chr,
                            ranges=IRanges(plot_range[1], plot_range[2]),
                            strand=GeneTxInfo$strand)
      exons_gr <- pintersect(exons_gr, segment_gr)
      exons_gr <- exons_gr[width(exons_gr)>0]
      if (length(exons_gr)==0) {
        next
      }
    }

    transcript_start <- min(start(exons_gr))
    transcript_end <- max(end(exons_gr))

    cds_ranges <- GeneTxInfo$xlimCds[[isoform]]
    if (!is.null(cds_ranges) && length(cds_ranges)>0) {
      if (!is.null(plot_range)) {
        segment_gr <- GRanges(seqnames=GeneTxInfo$chr,
                              ranges=IRanges(plot_range[1], plot_range[2]),
                              strand=GeneTxInfo$strand)
        cds_ranges <- pintersect(cds_ranges, segment_gr)
        cds_ranges <- cds_ranges[width(cds_ranges)>0]
      }

      if (length(cds_ranges)>0) {
        cds_df <- data.frame(
          start=start(cds_ranges),
          end=end(cds_ranges),
          y=y_value,
          feature="CDS",
          isoform=isoform,
          height_factor=1,
          orf_id=NA,
          row.names=NULL
        )
        isoform_data_list[[isoform_idx]] <- cds_df
        isoform_idx <- isoform_idx+1
      }
    }

    no_cds_for_this_iso <- (length(cds_ranges)==0)

    if (!no_cds_for_this_iso) {
      if (isoform %in% names(GeneTxInfo$fiveUTRByYFGtx)) {
        fiveUTR_gr <- unlist(GeneTxInfo$fiveUTRByYFGtx[isoform])
        if (length(fiveUTR_gr)>0) {
          if (strand == "+") {
            end(fiveUTR_gr) <- end(fiveUTR_gr)+1
          } else {
            start(fiveUTR_gr) <- start(fiveUTR_gr)-1
          }
          if (!is.null(plot_range)) {
            segment_gr <- GRanges(seqnames=GeneTxInfo$chr,
                                  ranges=IRanges(plot_range[1],plot_range[2]),
                                  strand=GeneTxInfo$strand)
            fiveUTR_gr <- pintersect(fiveUTR_gr,segment_gr)
            fiveUTR_gr <- fiveUTR_gr[width(fiveUTR_gr)>0]
          }
          if (length(fiveUTR_gr)>0) {
            fiveUTR_df <- data.frame(
              start=start(fiveUTR_gr),
              end=end(fiveUTR_gr),
              y=y_value,
              feature="5' UTR",
              isoform=isoform,
              height_factor=1,
              orf_id=NA,
              row.names=NULL
            )
            isoform_data_list[[isoform_idx]] <- fiveUTR_df
            isoform_idx <- isoform_idx+1
          }
        }
      }

      if (isoform %in% names(GeneTxInfo$threeUTRByYFGtx)) {
        threeUTR_gr <- unlist(GeneTxInfo$threeUTRByYFGtx[isoform])
        if (length(threeUTR_gr)>0) {
          if (strand == "+") {
            start(threeUTR_gr) <- start(threeUTR_gr)-1
          } else {
            end(threeUTR_gr) <- end(threeUTR_gr)+1
          }
          if (!is.null(plot_range)) {
            segment_gr <- GRanges(seqnames=GeneTxInfo$chr,
                                  ranges=IRanges(plot_range[1], plot_range[2]),
                                  strand=GeneTxInfo$strand)
            threeUTR_gr <- pintersect(threeUTR_gr, segment_gr)
            threeUTR_gr <- threeUTR_gr[width(threeUTR_gr)>0]
          }
          if (length(threeUTR_gr)>0) {
            threeUTR_df <- data.frame(
              start=start(threeUTR_gr),
              end=end(threeUTR_gr),
              y=y_value,
              feature="3' UTR",
              isoform=isoform,
              height_factor=1,
              orf_id=NA,
              row.names=NULL
            )
            isoform_data_list[[isoform_idx]] <- threeUTR_df
            isoform_idx <- isoform_idx+1
          }
        }
      }

    } else {
      # no CDS => ncRNA
      exons_plot <- exons_gr
      if (!is.null(plot_range)) {
        segment_gr <- GRanges(seqnames=GeneTxInfo$chr,
                              ranges=IRanges(plot_range[1],plot_range[2]),
                              strand=GeneTxInfo$strand)
        exons_plot <- pintersect(exons_plot, segment_gr)
        exons_plot <- exons_plot[width(exons_plot)>0]
      }
      if (length(exons_plot)>0) {
        ncRNA_df <- data.frame(
          start=start(exons_plot),
          end=end(exons_plot),
          y=y_value,
          feature="ncRNA",
          isoform=isoform,
          height_factor=1,
          orf_id=NA,
          row.names=NULL
        )
        isoform_data_list[[isoform_idx]] <- ncRNA_df
        isoform_idx <- isoform_idx+1
      }
    }

    intron_df_list <- list()
    if (length(exons_gr)>1) {
      exons_sorted <- exons_gr[order(start(exons_gr))]
      for (ii in seq_len(length(exons_sorted)-1)) {
        intron_start <- end(exons_sorted[ii])
        intron_end <- start(exons_sorted[ii+1])
        intron_df <- data.frame(
          xstart=intron_start,
          xend=intron_end,
          y=y_value,
          isoform=isoform,
          row.names=NULL
        )
        intron_df_list[[length(intron_df_list)+1]] <- intron_df
      }
    }

    if (!is.null(plot_range)) {
      segment_left <- plot_range[1]
      segment_right <- plot_range[2]
      exons_sorted <- exons_gr[order(start(exons_gr))]

      if (start(exons_sorted[1])>segment_left) {
        intron_df <- data.frame(
          xstart=segment_left,
          xend=start(exons_sorted[1]),
          y=y_value,
          isoform=isoform,
          row.names=NULL
        )
        intron_df_list[[length(intron_df_list)+1]] <- intron_df
      }

      if (end(exons_sorted[length(exons_sorted)])<segment_right) {
        intron_df <- data.frame(
          xstart=end(exons_sorted[length(exons_sorted)]),
          xend=segment_right,
          y=y_value,
          isoform=isoform,
          row.names=NULL
        )
        intron_df_list[[length(intron_df_list)+1]] <- intron_df
      }
    }

    if (length(intron_df_list)>0) {
      intron_data <- do.call(rbind,intron_df_list)
      line_data_list[[length(line_data_list)+1]] <- intron_data
    }

    if (length(isoform_data_list)>0) {
      isoform_df <- do.call(rbind, isoform_data_list)
      plot_data_list[[idx]] <- isoform_df
      idx <- idx+1
    }
  }

  if (length(plot_data_list)>0) {
    plot_data <- do.call(rbind, plot_data_list)
  } else {
    stop("No valid exons or features found to plot.")
  }

  if (length(line_data_list)>0) {
    line_data <- do.call(rbind,line_data_list)
  } else {
    line_data <- data.frame()
  }

  if (!"orf_id" %in% names(plot_data)) {
    plot_data$orf_id <- NA
  }
  plot_data$orf_id <- as.character(plot_data$orf_id)

  feature_order <- c("uORF","ouORF","nORF","odORF","dORF","5' UTR","CDS","3' UTR","ncRNA")
  plot_data$feature <- factor(plot_data$feature, levels=feature_order)

  plot_data$height <- 0.08 * plot_data$height_factor
  plot_data$ymin <- plot_data$y - plot_data$height
  plot_data$ymax <- plot_data$y + plot_data$height

  unique_features <- levels(plot_data$feature)[levels(plot_data$feature) %in% plot_data$feature]

  feature_colors <- c(
    "uORF"="yellow",
    "ouORF"="#FFD700",
    "nORF"="orange",
    "odORF"="#FFD700",
    "dORF"="yellow",
    "5' UTR"="lightgrey",
    "CDS"="black",
    "3' UTR"="white",
    "ncRNA"="#FFB6C1"
  )

  feature_colors <- feature_colors[unique_features]

  num_legend_items <- length(unique_features)
  legend_text_size <- 8
  base_key_size <-1
  if (num_legend_items>3) {
    key_size <- base_key_size * 3 / num_legend_items
  } else {
    key_size <- base_key_size
  }
  key_size <- max(0.8,key_size)

  p_gene <- ggplot()

  if (nrow(line_data)>0) {
    p_gene <- p_gene +
      geom_segment(data=line_data,aes(x=xstart,xend=xend,y=y,yend=y),color="black",inherit.aes=FALSE)
  }

  p_gene <- p_gene +
    geom_rect(data=plot_data[is.na(plot_data$orf_id), ],
              aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=feature),
              color="black", inherit.aes=FALSE)

  if (nrow(plot_data[!is.na(plot_data$orf_id), ])>0) {
    p_gene <- p_gene +
      geom_rect(data=plot_data[!is.na(plot_data$orf_id), ],
                aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=orf_id),
                color="black", inherit.aes=FALSE)
  }

  p_gene <- p_gene +
    scale_fill_manual(
      name="Feature",
      values=feature_colors,
      breaks=unique_features,
      labels=unique_features,
      guide=guide_legend(
        override.aes = list(fill=feature_colors),
        ncol=1,
        keyheight=unit(key_size,"lines"),
        keywidth=unit(1,"lines")
      )
    ) +
    theme_bw() +
    theme(
      legend.position="right",
      legend.title=element_blank(),
      legend.text=element_text(size=legend_text_size),
      legend.key.size=unit(key_size,"lines"),
      axis.text.y=element_text(size=transcript_label_font_size),
      axis.ticks.y=element_blank(),
      axis.title.y=element_text(size=12),
      plot.margin=unit(c(-0.2,0.2,0,1.5),"lines"),
      panel.grid=element_blank(),
      panel.border=element_blank()
    ) +
    xlab("Genomic Position") +
    ylab("") +
    coord_cartesian(clip="off")

  if (strand=="-") {
    p_gene <- p_gene + scale_x_reverse(limits=c(max(genelim), min(genelim)))
  } else {
    p_gene <- p_gene + scale_x_continuous(limits=c(min(genelim), max(genelim)))
  }

  labels <- sapply(isoform_positions$isoform, function(x) {
    if (x == tx_id) {
      paste0("bold('", x, "')")
    } else {
      paste0("'", x, "'")
    }
  })
  labels <- parse(text=labels)

  padding <- 0.1

  p_gene <- p_gene + scale_y_continuous(
    breaks=isoform_positions$y,
    labels=labels,
    limits=c(min(isoform_positions$y)-padding,
             max(isoform_positions$y)+padding)
  )

  return(p_gene)
}

#' Plot RNA-seq and Ribo-seq coverage for a gene
#'
#' The `ggRibo` function creates a comprehensive plot that displays RNA-seq coverage and Ribo-seq read counts for a specified gene and transcript. It also includes options to display extended ORFs (eORFs), genomic sequences, and gene models.
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
#' @param data_types Vector of sample names corresponding to data types for different SNR data (e.g., Ribo-seq (default). PARE-seq, CAGE-seq)
#' @param GRangeInfo Genomic range information, typically an object like `Txome_Range`.
#' @param RNAseqBamPaired Vector indicating whether each RNA-seq BAM file is 'paired' or 'single'-end.
#' @param Y_scale Character string, either "all" or "each", specifying how to scale the Y-axis for RNA-seq coverage. Default is "all".
#' @param Ribo_fix_height Numeric value to fix the maximum height of Ribo-seq counts in the plot.
#' @param plot_ORF_ranges Logical indicating whether to plot ORF ranges in the gene model. Default is FALSE.
#' @param oORF_coloring Character string specifying the method for coloring overlapping ORFs. Options are "oORF_colors" or "extend_mORF".
#' @param frame_colors Named vector of colors for the reading frames. Default is c("0"="#FF0000", "1"="#3366FF", "2"="#009900").
#' @param plot_range Optional. Numeric vector of length two specifying a custom genomic range to plot.
#' @param sample_color Vector specifying colors for each sample or "color" to use default coloring.
#' @param show_seq Logical indicating whether to display the DNA and amino acid sequence. Default is FALSE.
#' @param FASTA FASTA file containing genomic sequences.
#' @param dna_aa_height_ratio Numeric value to adjust the height of the DNA and amino acid sequence plot. Default is 0.5.
#' @param gene_model_height_ratio Numeric value to adjust the height of the gene model plot. If NULL, it adjusts automatically based on the number of transcripts.
#' @param transcript_label_font_size Optional numeric value to control the font size of the transcript ID labels in the gene model plot.
#'
#' @return A combined ggplot object displaying RNA-seq coverage, Ribo-seq counts, gene models, and optionally genomic sequences.
#' @examples
#' # Example usage:
#' # plot <- ggRibo(gene_id = "GENE1", tx_id = "TRANSCRIPT1", RNAseq = list("RNAseq.bam"), ...)
#' # print(plot)
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
                   oORF_coloring = NULL,
                   frame_colors = c("0"="#FF0000", "1"="#3366FF", "2"="#009900"),
                   plot_range = NULL,
                   sample_color = rep("color", length(Riboseq)),
                   show_seq = FALSE,
                   FASTA = NULL,
                   dna_aa_height_ratio = 0.5,
                   gene_model_height_ratio = NULL,
                   transcript_label_font_size = 10,
                   data_types = rep("Ribo-seq", length(SampleNames))) {

  # Validate data_types length
  if (length(data_types) != length(SampleNames)) {
    stop("The length of data_types must match the number of samples.")
  }

  # Validate Y_scale parameter
  if (!(Y_scale %in% c("all", "each"))) {
    stop("Invalid Y_scale value. Please choose either 'all' or 'each'.")
  }

  # Ensure that required annotations are loaded
  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }

  # Handle eORF annotation
  has_overlapping_ORF <- FALSE
  if (!is.null(eORF.tx_id)) {
    if (is.null(eORFRangeInfo)) {
      if (exists("eORF_Range", envir = .GlobalEnv)) {
        eORFRangeInfo <- get("eORF_Range", envir = .GlobalEnv)
      } else {
        stop("eORFRangeInfo (e.g., eORF_Range) must be provided when eORF.tx_id is specified.")
      }
    }
    missing_tx_ids <- setdiff(eORF.tx_id, names(eORFRangeInfo$eORFByTx))
    if (length(missing_tx_ids) > 0) {
      stop(paste("eORF Transcript IDs", paste(missing_tx_ids, collapse = ", "), "not found in eORFRangeInfo$eORFByTx."))
    }
  }

  # Create a Gene_info object for the gene of interest
  txByYFG <- GRangeInfo$txByGene[gene_id]

  if (length(txByYFG) == 0 || length(txByYFG[[1]]) == 0) {
    stop(paste("No transcripts found for gene ID", gene_id))
  }

  num_isoforms <- length(txByYFG[[1]])
  if (!"tx_name" %in% names(mcols(txByYFG[[1]]))) {
    stop("Transcript names ('tx_name') not found in GRangeInfo$txByGene. Please ensure 'tx_name' is a metadata column.")
  }
  tx_names <- txByYFG[[1]]$tx_name

  # Check which tx_names are in GRangeInfo$cdsByTx
  tx_names_in_cdsByTx <- intersect(tx_names, names(GRangeInfo$cdsByTx))

  # If all isoforms are noncoding
  if (length(tx_names_in_cdsByTx) == 0) {
    # all transcripts considered, likely no CDS => ncRNA scenario
    message("This is a noncoding gene (no annotated CDS for all isoforms). The frame is calculated from the 1st position of the transcript. Please provide the CDS ranges in a gtf if you have evidence that an ORF is translated.")
  } else {
    if (!(tx_id %in% tx_names_in_cdsByTx)) {
      # Main transcript no ORF
      message(paste("The transcript", tx_id, "does not have an annotated ORF. The frame is calculated from the 1st position of the transcript. Please provide the CDS ranges in a gtf if you have evidence that an ORF is translated."), call.=FALSE)
      # Ensure tx_id is included even if no CDS
      if (!(tx_id %in% tx_names)) {
        stop(paste("Transcript ID", tx_id, "not found in gene."))
      } else {
        tx_names <- unique(c(tx_id, tx_names))
      }
    }
  }

  strand_info <- as.character(strand(unlist(txByYFG)))[1]
  chr <- as.character(seqnames(unlist(txByYFG)))[1]

  # Sort transcripts
  other_tx_names <- setdiff(tx_names, tx_id)
  tx_names <- c(tx_id, sort(other_tx_names))

  # Extract CDS and exon information
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

  # Adjust CDS ranges (Keep as GRanges)
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

  # Define genomic range with extension or plot_range
  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)
    range_left <- plot_range[1]
    range_right <- plot_range[2]
    gene_ranges <- GRanges(seqnames=chr,
                           ranges=IRanges(range_left, range_right),
                           strand=strand_info)
  } else {
    gene_ranges <- reduce(unlist(txByYFG))

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

  # Filter Ribo-seq data
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

  main_cds <- xlimCds[[tx_id]]
  if (length(main_cds) > 0) {
    cds_left <- min(start(main_cds))
    cds_right <- max(end(main_cds))
  } else {
    cds_left <- NA
    cds_right <- NA
    # No second warning here, just rely on the one we printed if needed.
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

  # Handle eORF if provided
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

    RNAseq_list <- lapply(seq_len(length(RNAseq)), function(i) {
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
            Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][start(GeneTxInfo$generangesplus):end(GeneTxInfo$generangesplus)])
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
            Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][start(GeneTxInfo$generangesplus):end(GeneTxInfo$generangesplus)])
          }
        }
        Gtx1
      } else {
        stop(paste("Invalid value in RNAseqBamPaired for sample", SampleNames[i],
                   "- expected 'paired' or 'single', got", RNAseqBamPaired[i]))
      }
    })
  }

  if (length(RNAseq_list)>0) {
    max_Y_global <- max(unlist(RNAseq_list), na.rm=TRUE)
  } else {
    max_Y_global <- 0
  }

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

  if (!is.null(RNAseq)) {
    positions <- seq(start(GeneTxInfo$generangesplus), end(GeneTxInfo$generangesplus))

    for (i in seq_len(length(RNAseq))) {
      RNAseq_counts <- RNAseq_list[[i]]
      RNAseq_df <- data.frame(position=positions, count=RNAseq_counts, row.names=NULL)
      RNAseq_df <- RNAseq_df[!is.na(RNAseq_df$count), ]
      RNAseq_df$isoform <- tx_id

      if (length(Riboseq_list)>0) {
        RiboRslt <- Riboseq_list[[i]]
      } else {
        RiboRslt <- data.frame()
      }

      if (!is.null(Ribo_fix_height)) {
        current_max_Y <- max(RNAseq_counts,na.rm=TRUE)
        if (current_max_Y==0) {
          scale_factor_Ribo <-1
        } else {
          scale_factor_Ribo <- current_max_Y / Ribo_fix_height
        }
        y_limits <- c(0, current_max_Y*1.1)
      } else if (Y_scale=="all") {
        current_max_Y <- max_Y_global
        current_max_P <- max_P_global
        if (current_max_P>0) {
          scale_factor_Ribo <- max_Y_global / max_P_global
        } else {
          scale_factor_Ribo <-1
        }
        y_limits <- c(0,current_max_Y*1.1)
      } else if (Y_scale=="each") {
        current_max_Y <- max(RNAseq_counts, na.rm=TRUE)
        if (length(RiboRslt)>0 && nrow(RiboRslt)>0) {
          current_max_P <- max(RiboRslt$count,na.rm=TRUE)
          if (current_max_P>0) {
            scale_factor_Ribo <- current_max_Y / current_max_P
          } else {
            scale_factor_Ribo <-1
          }
        } else {
          scale_factor_Ribo <-1
        }
        y_limits <- c(0,current_max_Y*1.1)
      }

      if (length(RiboRslt)>0 && !is.null(scale_factor_Ribo)) {
        RiboRslt$count_scaled <- RiboRslt$count*scale_factor_Ribo
      }

      sample_color_i <- sample_color[i]

      p <- ggplot() +
        geom_col(data=RNAseq_df, aes(x=position, y=count), fill=RNAbackground, color=RNAbackground, na.rm=TRUE) +
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

      if (GeneTxInfo$strand=="-") {
        p <- p + scale_x_reverse(limits=c(GeneTxInfo$range_right,GeneTxInfo$range_left))
        x_limits <- c(GeneTxInfo$range_right,GeneTxInfo$range_left)
      } else {
        p <- p + scale_x_continuous(limits=c(GeneTxInfo$range_left,GeneTxInfo$range_right))
        x_limits <- c(GeneTxInfo$range_left,GeneTxInfo$range_right)
      }

      p <- p + xlab("")

      main_has_cds <- length(GeneTxInfo$xlimCds[[tx_id]])>0

      # Previously condition was (main_has_cds && plot_ORF_ranges),
      # minimal fix: remove && plot_ORF_ranges for always plotting vertical lines:
      if (main_has_cds) {
        if (!is.null(eORFTxInfo)) {
          x_min <- min(x_limits)
          x_max <- max(x_limits)
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
            line_color <- if (overlaps_CDS) "orange" else "green"

            if (GeneTxInfo$strand=="+") {
              start_pos <- eORF_left_pos
              end_pos <- eORF_right_pos
            } else {
              start_pos <- eORF_right_pos
              end_pos <- eORF_left_pos
            }

            if (!is.null(start_pos) && !is.na(start_pos) && start_pos>=x_min && start_pos<=x_max) {
              p <- p + geom_vline(xintercept=start_pos, linetype="solid", color=line_color, alpha=0.5)
            }
            if (!is.null(end_pos) && !is.na(end_pos) && end_pos>=x_min && end_pos<=x_max) {
              p <- p + geom_vline(xintercept=end_pos, linetype="dashed", color=line_color, alpha=0.5)
            }
          }
        }

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

      if (length(RiboRslt)>0 && nrow(RiboRslt)>0) {
        cds_ranges <- GeneTxInfo$cdsByYFGtx[[tx_id]]
        exons <- GeneTxInfo$exonByYFGtx[[tx_id]]

        if (length(cds_ranges)==0) {
          # Noncoding transcript: assign frames from first nucleotide
          if (GeneTxInfo$strand=="+") {
            exons_sorted <- sort(exons, decreasing=FALSE)
          } else {
            exons_sorted <- sort(exons, decreasing=TRUE)
          }

          positions <- integer(0)
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
            positions <- c(positions, pos)
            tx_positions <- c(tx_positions, tx_pos)
            cum_len <- cum_len + len
          }
          position_df <- data.frame(position=positions, tx_pos=tx_positions)
          position_df$frame <- factor((position_df$tx_pos - 1) %% 3, levels=c(0,1,2))

          RiboRslt <- merge(RiboRslt, position_df[, c("position","frame")], by="position", all.x=TRUE)

          if (!is.null(Ribo_fix_height)) {
            RiboRslt$count <- pmin(RiboRslt$count, Ribo_fix_height)
            RiboRslt$count_scaled <- RiboRslt$count * scale_factor_Ribo
          }

          if (sample_color_i=="color") {
            p <- p +
              geom_segment(data=RiboRslt,
                           aes(x=position, xend=position, y=0, yend=count_scaled, color=frame))
            p <- p + scale_color_manual(values=frame_colors, na.value="grey")
          } else {
            p <- p +
              geom_segment(data=RiboRslt,
                           aes(x=position, xend=position, y=0, yend=count_scaled),
                           color=sample_color_i)
          }

        } else {
          # Coding transcripts logic below (unchanged)
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
                if (length(overlaps_CDS) > 0) {
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
                  if (nrow(eORF_Riboseq) > 0) {
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
                  eORF_Riboseq <- eORF_Riboseq_list[[i]][[j]]

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

          } else if (fExtend>0 || tExtend>0) {
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

      p <- p + scale_y_continuous(
        limits=y_limits,
        name="RNA-seq \ncoverage",
        sec.axis=sec_axis(~ . / scale_factor_Ribo, name = paste0(data_types[i], "\n counts"))
      )

      p <- p + xlab("")

      delta_x <- 0
      x_label <- if (GeneTxInfo$strand=="-") {
        GeneTxInfo$range_right - delta_x
      } else {
        GeneTxInfo$range_left + delta_x
      }
      hjust_label <- 0

      y_label <- if (!is.null(Ribo_fix_height) || Y_scale=="all") {
        current_max_Y
      } else {
        current_max_Y + 0.01*current_max_Y
      }

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

      plot_list[[i]] <- p
    }
  }

  if (show_seq && !is.null(FASTA)) {
    dna_aa_plot <- plotDNAandAA(
      GeneTxInfo=GeneTxInfo,
      plot_range=plot_range,
      FASTA=FASTA
    )
  } else {
    dna_aa_plot <- NULL
  }

  gene_model_plot <- plotGeneTxModel(
    GeneTxInfo = GeneTxInfo,
    eORFTxInfo = eORFTxInfo,
    plot_ORF_ranges = plot_ORF_ranges,
    plot_range = plot_range,
    transcript_label_font_size = transcript_label_font_size
  )

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

#' Generate RNA-seq coverage plots for a given gene and transcript
#'
#' @param gene_id Character. The gene ID.
#' @param tx_id Character. The transcript ID.
#' @param Extend Numeric or numeric vector of length 2. Number of bases to extend beyond the gene region.
#'   If a single number is provided, it is used for both left and right extensions.
#'   If a vector of two numbers is provided, the first is used for the left extension and the second for the right.
#' @param NAME Character. An optional name for the gene.
#' @param RNAcoverline Character. Color for the RNA-seq coverage line. Default is "grey".
#' @param RNAbackground Character. Color for the RNA-seq coverage background. Default is "#FEFEAE".
#' @param RNAseq Character vector. Paths to RNA-seq BAM files.
#' @param SampleNames Character vector. Names of the samples corresponding to the RNA-seq data.
#' @param GRangeInfo List. Genomic ranges information, typically containing annotations such as `cdsByTx`, `txByGene`, etc.
#' @param RNAseqBamPaired Character vector. Specifies whether each RNA-seq BAM file is "paired" or "single" end.
#' @param Y_scale Character. Scaling method for the y-axis. Either "all" to use the same scale across all samples,
#'   or "each" to scale each sample individually. Default is "all".
#' @param plot_ORF_ranges Logical. Whether to plot the ORF ranges. Default is FALSE.
#' @param frame_colors Named character vector. Colors for different reading frames.
#' @param plot_range Numeric vector of length 2. Specifies a custom range to plot instead of the gene's range plus extensions.
#' @param sample_color Character vector. Colors for the samples. Default is "color" for all samples.
#' @param show_seq Logical. Whether to show the DNA and amino acid sequence. Default is FALSE.
#' @param FASTA Character. Path to the FASTA file containing genomic sequences.
#' @param dna_aa_height_ratio Numeric. Height ratio for the DNA and amino acid sequence plot. Default is 0.5.
#' @param gene_model_height_ratio Numeric. Height ratio for the gene model plot. If NULL, it is calculated based on the number of transcripts.
#' @param transcript_label_font_size Numeric. Font size for the transcript labels. Default is 10.
#'
#' @return A ggplot object combining RNA-seq coverage plots and gene model.
#'
#' @examples
#' # Example usage:
#' plot <- ggRNA(
#'   gene_id = "GENE1",
#'   tx_id = "TRANSCRIPT1",
#'   RNAseq = c("sample1.bam", "sample2.bam"),
#'   SampleNames = c("Sample1", "Sample2"),
#'   RNAseqBamPaired = c("paired", "single"),
#'   GRangeInfo = Txome_Range,
#'   Y_scale = "all",
#'   plot_ORF_ranges = FALSE
#' )
#' print(plot)
ggRNA <- function(gene_id, tx_id,
                  Extend = 100, NAME = "",
                  RNAcoverline = "grey", RNAbackground = "#FEFEAE",
                  RNAseq = RNAseqData,
                  SampleNames = Samples,
                  GRangeInfo = Txome_Range,
                  RNAseqBamPaired = RNAseqBamPairorSingle,
                  Y_scale = "all",
                  plot_ORF_ranges = FALSE,
                  frame_colors = c("0"="#FF0000", "1"="#3366FF", "2"="#009900"),
                  plot_range = NULL,
                  sample_color = rep("color", length(RNAseq)),
                  show_seq = FALSE,
                  FASTA = NULL,
                  dna_aa_height_ratio = 0.5,
                  gene_model_height_ratio = NULL,
                  transcript_label_font_size = 10) {

  # Validate Y_scale parameter
  if (!(Y_scale %in% c("all", "each"))) {
    stop("Invalid Y_scale value. Please choose either 'all' or 'each'.")
  }

  # Ensure that required annotations are loaded
  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }

  # Ensure the transcript ID exists
  if (!(tx_id %in% names(GRangeInfo$cdsByTx))) {
    stop(paste("Transcript ID", tx_id, "not found in GRangeInfo$cdsByTx."))
  }

  # Create Gene_info object for the gene of interest
  txByYFG <- GRangeInfo$txByGene[gene_id]

  if (length(txByYFG) == 0 || length(txByYFG[[1]]) == 0) {
    stop(paste("No transcripts found for gene ID", gene_id))
  }

  # Number of isoforms and their names
  num_isoforms <- length(txByYFG[[1]])
  if (!"tx_name" %in% names(mcols(txByYFG[[1]]))) {
    stop("Transcript names ('tx_name') not found in GRangeInfo$txByGene. Please ensure 'tx_name' is a metadata column.")
  }
  tx_names <- txByYFG[[1]]$tx_name

  # Strand and chromosome information
  strand_info <- as.character(strand(unlist(txByYFG)))[1]
  chr <- as.character(seqnames(unlist(txByYFG)))[1]

  # Sort transcript names and prioritize the main transcript
  other_tx_names <- setdiff(tx_names, tx_id)
  tx_names <- c(tx_id, sort(other_tx_names))

  # Extract CDS and exon information
  cdsByYFGtx <- GRangeInfo$cdsByTx[tx_names]
  exonByYFGtx <- GRangeInfo$exonsByTx[tx_names]

  # Adjust CDS ranges (Keep as GRanges)
  xlimCds <- list()
  for (i in seq_along(tx_names)) {
    cds <- cdsByYFGtx[[tx_names[i]]]
    if (length(cds) > 0) {
      xlimCds[[i]] <- cds  # Keep as GRanges
    } else {
      xlimCds[[i]] <- GRanges()
    }
  }
  names(xlimCds) <- tx_names

  # Define genomic range with extension or plot_range
  if (!is.null(plot_range)) {
    plot_range <- sort(plot_range)
    # Use the provided plot_range
    range_left <- plot_range[1]
    range_right <- plot_range[2]
    gene_ranges <- GRanges(seqnames = chr,
                           ranges = IRanges(range_left, range_right),
                           strand = strand_info)
  } else {
    # Use gene ranges, extended by 'Extend'
    gene_ranges <- reduce(unlist(txByYFG))

    # Determine Extend_left and Extend_right
    if (length(Extend) == 1) {
      Extend_left <- Extend
      Extend_right <- Extend
    } else if (length(Extend) == 2) {
      Extend_left <- Extend[1]
      Extend_right <- Extend[2]
    } else {
      stop("Extend must be a numeric value or a vector of two numeric values.")
    }

    # Adjust range_left and range_right based on strand
    if (strand_info == "+") {
      range_left <- min(start(gene_ranges)) - Extend_left
      range_right <- max(end(gene_ranges)) + Extend_right
    } else if (strand_info == "-") {
      range_left <- min(start(gene_ranges)) - Extend_right
      range_right <- max(end(gene_ranges)) + Extend_left
    } else {
      stop("Invalid strand information.")
    }

    gene_ranges <- GRanges(seqnames = chr,
                           ranges = IRanges(range_left, range_right),
                           strand = strand_info)
  }

  # Extract UTR information
  isoforms_w_3UTR <- tx_names[tx_names %in% names(GRangeInfo$threeUTR)]
  threeUTRByYFGtx <- GRangeInfo$threeUTR[isoforms_w_3UTR]

  isoforms_w_5UTR <- tx_names[tx_names %in% names(GRangeInfo$fiveUTR)]
  fiveUTRByYFGtx <- GRangeInfo$fiveUTR[isoforms_w_5UTR]

  # Assign cds_left and cds_right
  if (length(xlimCds[[tx_id]]) > 0) {
    cds_left <- min(start(xlimCds[[tx_id]]))
    cds_right <- max(end(xlimCds[[tx_id]]))
  } else {
    cds_left <- NA
    cds_right <- NA
  }

  # Create Gene_info object
  GeneTxInfo <- Gene_info$new(
    gene_id = gene_id,
    tx_id = tx_id,
    txByGene = txByYFG,
    cdsByYFGtx = cdsByYFGtx,
    chr = chr,
    generanges = gene_ranges,
    generangesplus = gene_ranges,  # Use gene_ranges directly
    range_left = range_left,
    range_right = range_right,
    num_isoforms = num_isoforms,
    tx_names = tx_names,
    isoforms.w.3UTR = isoforms_w_3UTR,
    isoforms.w.5UTR = isoforms_w_5UTR,
    threeUTRByYFGtx = threeUTRByYFGtx,
    fiveUTRByYFGtx = fiveUTRByYFGtx,
    exonByYFGtx = exonByYFGtx,
    Extend = Extend,
    strand = strand_info,
    xlimCds = xlimCds,
    Riboseq_list = NULL,  # Set Riboseq_list to NULL to fix the error
    cds_left = cds_left,
    cds_right = cds_right
  )

  # Set up BAM parameters for reading RNA-seq data
  if (!is.null(RNAseq)) {
    what1 <- c("rname", "strand", "pos", "qwidth", "seq")
    param <- ScanBamParam(which = GeneTxInfo$generangesplus, what = what1)
  }

  # Read RNA-seq data
  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    # Validate RNAseqBamPaired
    if (is.null(RNAseqBamPaired)) {
      stop("RNAseqBamPaired must be provided when RNAseq data is supplied.")
    }

    if (length(RNAseqBamPaired) != length(RNAseq)) {
      stop("Length of RNAseqBamPaired must match the number of RNAseq samples.")
    }

    # Read each RNA-seq BAM file and compute coverage
    RNAseq_list <- lapply(seq_len(length(RNAseq)), function(i) {
      if (RNAseqBamPaired[i] == "paired") {
        # Read paired-end alignments
        readPairs1 <- suppressWarnings(
          readGAlignmentPairs(RNAseq[i], param = param, strandMode = 2)
        )
        if (length(readPairs1) == 0) {
          warning(paste("No paired-end reads found in", RNAseq[i]))
          Gtx1 <- numeric(length = width(GeneTxInfo$generangesplus))
        } else {
          # Filter reads by strand
          readPairs1 <- readPairs1[strand(readPairs1) == GeneTxInfo$strand]
          if (length(readPairs1) == 0) {
            warning(paste("No reads matching strand", GeneTxInfo$strand, "found in", RNAseq[i]))
            Gtx1 <- numeric(length = width(GeneTxInfo$generangesplus))
          } else {
            # Compute coverage over the gene region
            cvg1 <- coverage(readPairs1)
            Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][start(GeneTxInfo$generangesplus):end(GeneTxInfo$generangesplus)])
          }
        }
        Gtx1
      } else if (RNAseqBamPaired[i] == "single") {
        # Read single-end alignments
        alignments <- suppressWarnings(readGAlignments(RNAseq[i], param = param))
        if (length(alignments) == 0) {
          warning(paste("No single-end reads found in", RNAseq[i]))
          Gtx1 <- numeric(length = width(GeneTxInfo$generangesplus))
        } else {
          # Filter reads by strand
          alignments <- alignments[strand(alignments) == GeneTxInfo$strand]
          if (length(alignments) == 0) {
            warning(paste("No reads matching strand", GeneTxInfo$strand, "found in", RNAseq[i]))
            Gtx1 <- numeric(length = width(GeneTxInfo$generangesplus))
          } else {
            # Compute coverage over the gene region
            cvg1 <- coverage(alignments)
            Gtx1 <- as.numeric(cvg1[[GeneTxInfo$chr]][start(GeneTxInfo$generangesplus):end(GeneTxInfo$generangesplus)])
          }
        }
        Gtx1
      } else {
        stop(paste("Invalid value in RNAseqBamPaired for sample", SampleNames[i],
                   "- expected 'paired' or 'single', got", RNAseqBamPaired[i]))
      }
    })
  }

  # Determine maximum RNA-seq coverage across all samples
  if (length(RNAseq_list) > 0) {
    max_Y_global <- max(unlist(RNAseq_list), na.rm = TRUE)
  } else {
    max_Y_global <- 0
  }

  # Initialize list to store individual plots
  plot_list <- list()

  # Generate plots for RNA-seq data
  if (!is.null(RNAseq)) {
    # Define positions
    positions <- seq(start(GeneTxInfo$generangesplus), end(GeneTxInfo$generangesplus))

    for (i in seq_len(length(RNAseq))) {
      RNAseq_counts <- RNAseq_list[[i]]
      RNAseq_df <- data.frame(position = positions, count = RNAseq_counts, row.names = NULL)

      # Remove rows with NA counts
      RNAseq_df <- RNAseq_df[!is.na(RNAseq_df$count), ]

      # Prepare data for plotting
      RNAseq_df$isoform <- tx_id  # Assign the isoform to the data

      # Determine scaling based on Y_scale parameter
      if (Y_scale == "all") {
        current_max_Y <- max_Y_global
        y_limits <- c(0, current_max_Y * 1.1)
      } else if (Y_scale == "each") {
        current_max_Y <- max(RNAseq_counts, na.rm = TRUE)
        y_limits <- c(0, current_max_Y * 1.1)
      }

      # Get sample color for current sample
      sample_color_i <- sample_color[i]

      # Create the plot
      p <- ggplot() +
        # RNA-seq data: plot as filled columns and line
        geom_col(data = RNAseq_df, aes(x = position, y = count), fill = RNAbackground, color = RNAbackground, na.rm = TRUE) +
        geom_step(data = RNAseq_df, aes(x = position, y = count), color = RNAcoverline, na.rm = TRUE) +
        theme_bw() +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",  # Remove legend to maintain consistent width
          plot.margin = unit(c(0, 0.2, 0, 0.2), "lines"),
          panel.grid.major.x = element_blank(),  # Remove vertical grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.3),  # Light horizontal grid lines
          axis.title.y = element_text(size = 10),
          panel.background = element_rect(fill = "white", color = NA)  # Set background to white
        )

      # Set x-axis scales based on strand and limits
      if (GeneTxInfo$strand == "-") {
        p <- p + scale_x_reverse(limits = c(GeneTxInfo$range_right, GeneTxInfo$range_left))
        x_limits <- c(GeneTxInfo$range_right, GeneTxInfo$range_left)
      } else {
        p <- p + scale_x_continuous(limits = c(GeneTxInfo$range_left, GeneTxInfo$range_right))
        x_limits <- c(GeneTxInfo$range_left, GeneTxInfo$range_right)
      }

      # Remove x-axis label to save space
      p <- p + xlab("")

      x_min <- min(x_limits)
      x_max <- max(x_limits)

      # Add vertical lines for main ORF start (black dashed) and stop (grey dashed)
      # Determine main ORF start and stop positions based on strand
      main_orf_start <- if (GeneTxInfo$strand == "+") GeneTxInfo$cds_left else GeneTxInfo$cds_right
      main_orf_stop <- if (GeneTxInfo$strand == "+") GeneTxInfo$cds_right else GeneTxInfo$cds_left

      # Add vertical lines with specified colors, checking for NA and plotting range
      if (!is.na(main_orf_start) && main_orf_start >= x_min && main_orf_start <= x_max) {
        p <- p + geom_vline(xintercept = main_orf_start, linetype = "dashed", color = "black")
      }
      if (!is.na(main_orf_stop) && main_orf_stop >= x_min && main_orf_stop <= x_max) {
        p <- p + geom_vline(xintercept = main_orf_stop, linetype = "dashed", color = "darkgrey")
      }

      # Set y-axis limits and labels with adjusted limits to accommodate sample names
      p <- p + scale_y_continuous(
        limits = y_limits,  # Set limits based on Y_scale
        name = "RNA-seq \ncoverage"
      )

      # Remove x-axis label to save space
      p <- p + xlab("")

      # Adjust x position and hjust for sample name annotation
      delta_x <- 0  # Position sample names at the extreme left/right based on strand
      x_label <- if (GeneTxInfo$strand == "-") {
        GeneTxInfo$range_right - delta_x
      } else {
        GeneTxInfo$range_left + delta_x
      }
      hjust_label <- 0  # Left aligned

      # Adjust y position to move the label further up based on Y_scale
      y_label <- if (Y_scale == "all") {
        current_max_Y
      } else {
        current_max_Y + 0.01 * current_max_Y
      }

      # Add sample name inside the plot
      p <- p + annotate("text",
                        x = x_label,
                        y = y_label,
                        label = SampleNames[i],
                        hjust = hjust_label,
                        vjust = 0,  # Align text at the bottom
                        size = 3,
                        fontface = "bold")

      # Adjust theme to reduce empty space and align labels
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",  # Suppress legend
        plot.margin = unit(c(0, 0.2, -0.8, 0.2), "lines"),
        panel.grid.major.x = element_blank(),  # Remove vertical grid lines
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.3),
        axis.title.y = element_text(size = 10)
      )

      # Store the ggplot object in plot_list
      plot_list[[i]] <- p
    }
  }

  # Create the DNA and amino acid plot if show_seq is TRUE
  if (show_seq && !is.null(FASTA)) {
    dna_aa_plot <- plotDNAandAA(
      GeneTxInfo = GeneTxInfo,
      plot_range = plot_range,
      FASTA = FASTA
    )
  } else {
    dna_aa_plot <- NULL
  }

  # Create the gene model plot
  gene_model_plot <- plotGeneTxModel(
    GeneTxInfo = GeneTxInfo,
    plot_ORF_ranges = plot_ORF_ranges,
    plot_range = plot_range,
    transcript_label_font_size = transcript_label_font_size  # Pass the new parameter
  )

  # Adjust the height of the gene model plot dynamically based on the number of transcripts
  num_transcripts <- GeneTxInfo$num_isoforms
  num_datasets <- ifelse(!is.null(RNAseq), length(RNAseq), 0)

  # Define relative heights for the title, RNA-seq plots, gene model plot, and DNA/aa plot
  title_height <- 0.2  # Title height
  rna_height <- 0.8  # RNA-seq plot height

  # Adjust gene model height calculation based on number of transcripts
  if (is.null(gene_model_height_ratio)) {
    # Adjust automatically based on number of transcripts
    gene_model_height_ratio <- 0.2 + (num_transcripts) * 0.1
  }
  gene_model_height <- gene_model_height_ratio  # Use the computed or user-specified ratio

  # Adjust for DNA sequence if show_seq is TRUE
  if (show_seq && !is.null(FASTA)) {
    dna_aa_height <- dna_aa_height_ratio  # Height for DNA/aa plot, adjusted by dna_aa_height_ratio
  } else {
    dna_aa_height <- 0
  }

  # Total height units is the sum of heights of all plots
  total_height_units <- title_height + (num_datasets * rna_height) + dna_aa_height + gene_model_height

  # Adjust relative heights to proportionally allocate space
  rel_heights <- c(
    title_height,
    rep(rna_height, num_datasets),
    dna_aa_height,
    gene_model_height
  ) / total_height_units

  # Create a title plot using ggplot2
  title_plot <- ggplot() +
    theme_void() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "lines")
    ) +
    annotate("text",
             x = 0.5, y = 0.5,
             label = paste(gene_id, "  ", NAME),
             hjust = 0.5, vjust = 0.5,
             fontface = "italic", size = 5)

  # Combine the plots using cowplot::plot_grid
  combined_plot <- cowplot::plot_grid(
    title_plot,
    plotlist = c(plot_list, list(dna_aa_plot), list(gene_model_plot)),
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
