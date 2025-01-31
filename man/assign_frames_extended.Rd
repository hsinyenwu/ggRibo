% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/code.R
\name{assign_frames_extended}
\alias{assign_frames_extended}
\title{Assign Frames with Extended CDS Ranges}
\usage{
assign_frames_extended(Ribo_data, extended_cds_ranges, strand, main_cds_ranges)
}
\arguments{
\item{Ribo_data}{Data frame containing Ribo-seq reads.}

\item{extended_cds_ranges}{Extended CDS ranges as a \code{GRanges} object.}

\item{strand}{Strand information ("+" or "-").}

\item{main_cds_ranges}{Main CDS ranges as a \code{GRanges} object.}
}
\value{
Data frame with an added \code{frame} column indicating the reading frame.
}
\description{
Assigns reading frames to Ribo-seq reads based on extended CDS ranges, including overlapping eORFs.
}
