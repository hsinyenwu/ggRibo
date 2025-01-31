% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/code.R
\name{plotDNAandAA}
\alias{plotDNAandAA}
\title{Plot DNA and Amino Acid Sequences}
\usage{
plotDNAandAA(GeneTxInfo, plot_range = NULL, FASTA = NULL)
}
\arguments{
\item{GeneTxInfo}{A \code{Gene_info} object containing gene-specific information.}

\item{plot_range}{Optional vector specifying the genomic range to plot.}

\item{FASTA}{A \code{BSgenome} object containing the reference genome sequences.}
}
\value{
A \code{ggplot2} object representing the DNA and amino acid sequences.
}
\description{
Generates a plot showing DNA nucleotides and their corresponding amino acids for the genomic region of interest.
}
