# -----------------------------------
# ggRibo_decon function
# -----------------------------------
#'
#' `ggRibo_decon` creates a combined visualization of RNA-Seq coverage and frame-specific Ribo-Seq counts for a specified gene and transcript. It generates three separate plots corresponding to the three reading frames (0, 1, and 2) of Ribo-Seq data, optionally including reads that do not fall within any annotated ORF regions.
#'
#' @param gene_id Character. The identifier for the gene of interest.
#' @param tx_id Character. The transcript identifier within the gene for which the main ORF is annotated.
#' @param eORF.tx_id Character vector, optional. Transcript identifiers for extended ORFs (eORFs) associated with the gene. Defaults to `NULL`.
#' @param eORFRangeInfo List, optional. A list containing eORF ranges, typically obtained from annotation data. Required if `eORF.tx_id` is specified. Defaults to `NULL`.
#' @param Extend Numeric or numeric vector of length 2. Extends the plotting range upstream and downstream of the gene. If a single value is provided, it is applied to both ends. Defaults to `100`.
#' @param NAME Character. An optional name to annotate the plot title. Defaults to an empty string `""`.
#' @param RNAcoverline Character. Color for the RNA-Seq coverage line. Defaults to `"grey"`.
#' @param RNAbackground Character or character vector. Fill color(s) for the RNA-Seq coverage bars. If a single color is provided, it is replicated for all samples. Defaults to `"#FEFEAE"`.
#' @param fExtend Numeric. Extends the ORF upstream by the specified number of bases. Defaults to `0`.
#' @param tExtend Numeric. Extends the ORF downstream by the specified number of bases. Defaults to `0`.
#' @param RNAseq Character vector. Paths to RNA-Seq BAM files. Defaults to `RNAseqData`.
#' @param Riboseq Character vector. Paths to Ribo-Seq BAM files. Defaults to `RiboseqData`.
#' @param SampleNames Character vector. Names of the samples corresponding to the Ribo-Seq data. Defaults to `Samples`.
#' @param GRangeInfo List. A list containing genomic range information, such as transcripts, CDS, exons, UTRs, etc. Typically obtained from annotation data. Defaults to `Txome_Range`.
#' @param RNAseqBamPaired Character vector. Indicates whether each RNA-Seq BAM file is paired-end (`"paired"`) or single-end (`"single"`). Must match the length of `RNAseq`. Defaults to `RNAseqBamPairorSingle`.
#' @param Y_scale Character. Determines the Y-axis scaling for Ribo-Seq counts. `"all"` uses a global scale across all samples, while `"each"` scales each sample individually. Defaults to `"all"`.
#' @param Ribo_fix_height Numeric, optional. If provided, Ribo-Seq counts are capped at this value, and `Y_scale` is ignored. Defaults to `NULL`.
#' @param plot_ORF_ranges Logical. If `TRUE`, annotated ORF ranges are highlighted on the gene model plot. Defaults to `FALSE`.
#' @param oORF_coloring Character, optional. Determines the coloring scheme for overlapping ORFs. Options include `"oORF_colors"`, `"extend_mORF"`, or `NULL` for default coloring. Defaults to `NULL`.
#' @param frame_colors Named character vector. Specifies colors for the three reading frames. Must have names `"0"`, `"1"`, and `"2"`. Defaults to `c("0"="#FF0000", "1"="#3366FF", "2"="#009900")`.
#' @param plot_range Numeric vector of length 2, optional. Specifies the exact genomic range to plot. If provided, `Extend` is ignored. Defaults to `NULL`.
#' @param sample_color Character. Color for the Ribo-Seq counts. If set to `"color"`, frame-specific colors are used. Otherwise, a single color is applied to all frames. Defaults to `"color"`.
#' @param show_seq Logical. If `TRUE`, DNA and amino acid sequences are displayed below the coverage plots. Defaults to `FALSE`.
#' @param FASTA Character, optional. Path to a FASTA file containing the genomic sequences. Required if `show_seq` is `TRUE`. Defaults to `NULL`.
#' @param dna_aa_height_ratio Numeric. Adjusts the height ratio of the DNA and amino acid sequence plot relative to the gene model plot. Defaults to `0.5`.
#' @param gene_model_height_ratio Numeric, optional. Adjusts the height ratio of the gene model plot. If `NULL`, it is automatically adjusted based on the number of transcripts. Defaults to `NULL`.
#' @param transcript_label_font_size Numeric. Controls the font size of the transcript ID labels in the gene model plot. Defaults to `10`.
#' @param plot_genomic_direction Logical. If `TRUE`, an arrow indicating the genomic strand direction is added to the top RNA-Seq plot. Defaults to `FALSE`.
#' @param data_types Character vector. Describes the type of data being plotted (e.g., `"Ribo-seq"`). Must match the length of `SampleNames`. Defaults to `"Ribo-seq"`.
#' @param plot_unassigned_reads Logical. If `TRUE`, Ribo-Seq reads that do not fall within any annotated or extended ORF regions are plotted as grey segments. Defaults to `FALSE`.
#'
#' @return A combined `ggplot` object displaying RNA-Seq coverage, frame-specific Ribo-Seq counts, gene models, and optionally genomic sequences.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' plot <- ggRibo_decon(
#'   gene_id = "GENE1",
#'   tx_id = "TRANSCRIPT1",
#'   RNAseq = list("RNAseq_sample1.bam"),
#'   Riboseq = list("Ribo_sample1.bam"),
#'   SampleNames = c("Sample1"),
#'   plot_unassigned_reads = TRUE
#' )
#' print(plot)
#' }
#'
#' @export

ggRibo_decon <- function(gene_id, tx_id, eORF.tx_id = NULL,
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
                       plot_unassigned_reads = TRUE) {
  # Validate parameters
  if (length(data_types) != length(SampleNames)) {
    stop("The length of data_types must match the number of samples.")
  }

  if (!(Y_scale %in% c("all", "each"))) {
    stop("Invalid Y_scale value. Please choose either 'all' or 'each'.")
  }

  if (length(RNAbackground) == 1) {
    RNAbackground <- rep(RNAbackground, length(SampleNames))
  } else if (length(RNAbackground) != length(SampleNames)) {
    stop("RNAbackground must be either a single color or a vector of colors with the same length as 'Samples'.")
  }

  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }

  # Check if multiple samples are provided
  if (length(SampleNames) > 1) {
    stop("ggRibo_decon only supports one sample at a time. Please provide a single RNA-seq and a single Ribo-seq sample.")
  }

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

  txByYFG <- GRangeInfo$txByGene[gene_id]

  if (length(txByYFG) == 0 || length(txByYFG[[1]]) == 0) {
    stop(paste("No transcripts found for gene ID", gene_id))
  }

  num_isoforms <- length(txByYFG[[1]])
  if (!"tx_name" %in% names(mcols(txByYFG[[1]]))) {
    stop("Transcript names ('tx_name') not found in GRangeInfo$txByGene.")
  }
  tx_names <- txByYFG[[1]]$tx_name

  tx_names_in_cdsByTx <- intersect(tx_names, names(GRangeInfo$cdsByTx))

  if (length(tx_names_in_cdsByTx) == 0) {
    message("This is a noncoding gene. The frame is calculated from the 1st position of the transcript.")
  } else {
    if (!(tx_id %in% tx_names_in_cdsByTx)) {
      message(paste("The transcript", tx_id, "does not have an annotated ORF. Frame is from 1st transcript position."))
      if (!(tx_id %in% tx_names)) {
        stop(paste("Transcript ID", tx_id, "not found in gene."))
      } else {
        tx_names <- unique(c(tx_id, tx_names))
      }
    }
  }

  strand_info <- as.character(strand(unlist(txByYFG)))[1]
  chr <- as.character(seqnames(unlist(txByYFG)))[1]

  other_tx_names <- setdiff(tx_names, tx_id)
  tx_names <- c(tx_id, sort(other_tx_names))

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
      stop("Extend must be numeric of length 1 or 2.")
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

    gene_ranges <- GRanges(seqnames=chr,
                           ranges=IRanges(range_left, range_right),
                           strand=strand_info)
  }

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
      stop("Length of RNAseqBamPaired must match the number of RNAseq samples.")
    }

    global_start <- min(start(GeneTxInfo$generangesplus))
    global_end <- max(end(GeneTxInfo$generangesplus))

    i <- 1
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

  i <- 1
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

    sample_color_i <- sample_color
    x_limits <- if (GeneTxInfo$strand=="-") {
      c(GeneTxInfo$range_right,GeneTxInfo$range_left)
    } else {
      c(GeneTxInfo$range_left,GeneTxInfo$range_right)
    }

    main_has_cds <- length(GeneTxInfo$xlimCds[[tx_id]])>0
    cds_ranges <- GeneTxInfo$cdsByYFGtx[[tx_id]]
    exons <- GeneTxInfo$exonByYFGtx[[tx_id]]

    make_frame_plot <- function(RNAseq_df, Ribo_df, frame_color, frame_label, y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned) {
      p <- ggplot() +
        geom_col(data=RNAseq_df, aes(x=position, y=count), fill=RNAbackground[1], color=RNAbackground[1], na.rm=TRUE) +
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

      if (plot_unassigned && nrow(na_data)>0) {
        p <- p + geom_segment(data=na_data, aes(x=position, xend=position, y=0, yend=count_scaled), color="grey", na.rm=TRUE)
      }

      if (nrow(Ribo_df)>0) {
        p <- p + geom_segment(data=Ribo_df, aes(x=position, xend=position, y=0, yend=count_scaled), color=frame_color, na.rm=TRUE)
      }

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

          line_color <- if (overlaps_CDS) "orange" else "green"

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
        sec.axis=sec_axis(~ . / scale_factor_Ribo, name = paste0(data_types[1], "\n counts"))
      )

      return(p)
    }

    if (length(cds_ranges)==0) {
      # Noncoding frame assignment
      positions_all <- integer(0)
      tx_positions <- integer(0)
      if (GeneTxInfo$strand=="+") {
        exons_sorted <- sort(exons,decreasing=FALSE)
      } else {
        exons_sorted <- sort(exons,decreasing=TRUE)
      }
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
      # Coding
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

    frame0_data <- RiboRslt[RiboRslt$frame=="0", ]
    frame1_data <- RiboRslt[RiboRslt$frame=="1", ]
    frame2_data <- RiboRslt[RiboRslt$frame=="2", ]
    na_data <- RiboRslt[is.na(RiboRslt$frame),]

    p0 <- make_frame_plot(RNAseq_df, frame0_data, frame_colors["0"], "Frame0", y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned_reads)
    p1 <- make_frame_plot(RNAseq_df, frame1_data, frame_colors["1"], "Frame1", y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned_reads)
    p2 <- make_frame_plot(RNAseq_df, frame2_data, frame_colors["2"], "Frame2", y_limits, scale_factor_Ribo, GeneTxInfo, main_has_cds, eORFTxInfo, fExtend, tExtend, na_data, plot_unassigned_reads)

    delta_x <- 0
    x_label <- if (GeneTxInfo$strand=="-") {
      GeneTxInfo$range_right - delta_x
    } else {
      GeneTxInfo$range_left + delta_x
    }
    y_label <- current_max_Y + 0.01*current_max_Y

    # Annotate each plot with its frame instead of sample name
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

    p2 <- p2 + theme(plot.margin=unit(c(0,0.2,-0.8,0.2),"lines"))

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
# Modified ggRibo function to ggRNA
# -----------------------------------

#' Plot RNA-seq coverage for a gene
#'
#' The `ggRNA` function creates a comprehensive plot that displays RNA-seq coverage for a specified gene and transcript.
#' It optionally includes genomic sequences and gene models.
#'
#' @param gene_id Character string specifying the gene ID of interest.
#' @param tx_id Character string specifying the transcript ID to be used as the main isoform.
#' @param Extend Numeric value specifying the number of base pairs to extend the plot beyond the gene range. Default is 100.
#' @param NAME Optional. Character string for an additional name or title to display in the plot.
#' @param RNAcoverline Color for the RNA-seq coverage line. Default is "grey".
#' @param RNAbackground Color for the RNA-seq coverage background. Default is "#FEFEAE".
#' @param RNAseq List of file paths to RNA-seq BAM files.
#' @param SampleNames Vector of sample names corresponding to the RNAseq data.
#' @param GRangeInfo Genomic range information, typically an object like `Txome_Range`.
#' @param RNAseqBamPaired Vector indicating whether each RNA-seq BAM file is 'paired' or 'single'-end.
#' @param Y_scale Character string, either "all" or "each", specifying how to scale the Y-axis for RNA-seq coverage. Default is "all".
#' @param plot_ORF_ranges Logical indicating whether to plot ORF ranges in the gene model. Default is FALSE.
#' @param plot_range Optional. Numeric vector of length two specifying a custom genomic range to plot.
#' @param show_seq Logical indicating whether to display the DNA and amino acid sequence. Default is FALSE.
#' @param FASTA FASTA file containing genomic sequences.
#' @param plot_genomic_direction Plot the original direction of the gene in the genome according to the annotation.
#' @param dna_aa_height_ratio Numeric value to adjust the height of the DNA and amino acid sequence plot. Default is 0.5.
#' @param gene_model_height_ratio Numeric value to adjust the height of the gene model plot. If NULL, it adjusts automatically.
#' @param transcript_label_font_size Optional numeric value to control the font size of the transcript ID labels in the gene model plot.
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
                  plot_genomic_direction = FALSE) {

  if (!(Y_scale %in% c("all", "each"))) {
    stop("Invalid Y_scale value. Please choose either 'all' or 'each'.")
  }

  if (length(RNAbackground) == 1) {
    RNAbackground <- rep(RNAbackground, length(SampleNames))
  } else if (length(RNAbackground) != length(SampleNames)) {
    stop("RNAbackground must be either a single color or a vector of colors with the same length as 'Samples'.")
  }

  if (is.null(GRangeInfo)) {
    stop("GRangeInfo (e.g., Txome_Range) must be provided.")
  }

  txByYFG <- GRangeInfo$txByGene[gene_id]

  if (length(txByYFG) == 0 || length(txByYFG[[1]]) == 0) {
    stop(paste("No transcripts found for gene ID", gene_id))
  }

  num_isoforms <- length(txByYFG[[1]])
  if (!"tx_name" %in% names(mcols(txByYFG[[1]]))) {
    stop("Transcript names ('tx_name') not found in GRangeInfo$txByGene. Please ensure 'tx_name' is a metadata column.")
  }
  tx_names <- txByYFG[[1]]$tx_name

  tx_names_in_cdsByTx <- intersect(tx_names, names(GRangeInfo$cdsByTx))
  if (length(tx_names_in_cdsByTx) == 0) {
    message("This is a noncoding gene (no annotated CDS for all isoforms).")
  } else {
    if (!(tx_id %in% tx_names_in_cdsByTx)) {
      message(paste("The transcript", tx_id, "does not have an annotated ORF."))
      if (!(tx_id %in% tx_names)) {
        stop(paste("Transcript ID", tx_id, "not found in gene."))
      } else {
        tx_names <- unique(c(tx_id, tx_names))
      }
    }
  }

  strand_info <- as.character(strand(unlist(txByYFG)))[1]
  chr <- as.character(seqnames(unlist(txByYFG)))[1]
  other_tx_names <- setdiff(tx_names, tx_id)
  tx_names <- c(tx_id, sort(other_tx_names))

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

  # Only RNA-seq from here
  if (!is.null(RNAseq)) {
    what1 <- c("rname","strand","pos","qwidth","seq")
    param <- ScanBamParam(which=gene_ranges, what=what1)
  }

  RNAseq_list <- list()
  if (!is.null(RNAseq)) {
    if (is.null(RNAseqBamPaired)) {
      stop("RNAseqBamPaired must be provided when RNAseq data is supplied.")
    }

    if (length(RNAseqBamPaired) != length(RNAseq)) {
      stop("Length of RNAseqBamPaired must match the number of RNAseq samples.")
    }

    global_start <- min(start(gene_ranges))
    global_end <- max(end(gene_ranges))

    RNAseq_list <- lapply(seq_len(length(RNAseq)), function(i) {
      if (RNAseqBamPaired[i] == "paired") {
        readPairs1 <- suppressWarnings(
          readGAlignmentPairs(RNAseq[i], param=param, strandMode=2)
        )
        if (length(readPairs1)==0) {
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

  plot_list <- list()

  if (!is.null(RNAseq)) {
    global_start <- min(start(gene_ranges))
    global_end <- max(end(gene_ranges))
    positions <- seq(global_start, global_end)

    for (i in seq_len(length(RNAseq))) {
      RNAseq_counts <- RNAseq_list[[i]]
      RNAseq_df <- data.frame(position=positions, count=RNAseq_counts, row.names=NULL)
      RNAseq_df <- RNAseq_df[!is.na(RNAseq_df$count), ]
      RNAseq_df$isoform <- tx_id

      current_max_Y <- if (Y_scale=="all") max_Y_global else max(RNAseq_counts, na.rm=TRUE)
      if (current_max_Y == 0) {
        scale_factor_Ribo <- 1
      } else {
        # No Ribo data, so no scaling needed. Just set scale_factor_Ribo=1.
        scale_factor_Ribo <- 1
      }
      y_limits <- c(0,current_max_Y*1.1)

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

      if (strand_info=="-") {
        p <- p + scale_x_reverse(limits=c(range_right,range_left))
        x_limits <- c(range_right,range_left)
      } else {
        p <- p + scale_x_continuous(limits=c(range_left,range_right))
        x_limits <- c(range_left,range_right)
      }

      p <- p + xlab("")
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

      p <- p + scale_y_continuous(
        limits=y_limits,
        name="RNA-seq \ncoverage"
      )

      p <- p + xlab("")

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

      plot_list[[i]] <- p
    }
  }

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
