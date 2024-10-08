% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{make_se_customized}
\alias{make_se_customized}
\title{Data.frame to SummarizedExperiment object
conversion using an experimental design}
\usage{
make_se_customized(
  proteins_unique,
  columns,
  expdesign,
  log2transform = F,
  exp = "LFQ",
  lfq_type = NULL,
  level = NULL,
  exp_type = NULL
)
}
\arguments{
\item{proteins_unique}{Data.frame,
Protein table with unique names annotated in the 'name' column
(output from \code{\link{make_unique}()}).}

\item{columns}{Integer vector,
Column numbers indicating the columns containing the assay data.}

\item{expdesign}{Data.frame,
Experimental design with 'label', 'condition' and 'replicate' information.}

\item{exp}{quantification method i.e. LFQ, TMT, or DIA}

\item{level}{which level of the quantification table summarized at. For example, protein or peptide}
}
\value{
A SummarizedExperiment object
with log2-transformed values.
}
\description{
\code{make_se_customized} creates a SummarizedExperiment object
based on two data.frames: the protein table and experimental design.
}
