% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/e_cellType_gene_max.R
\name{e_cellType_gene_max}
\alias{e_cellType_gene_max}
\title{e_cellType_gene_max}
\usage{
e_cellType_gene_max(anno_interactions_allSamples, verbose = FALSE)
}
\arguments{
\item{anno_interactions_allSamples}{a list of data frames containing annotation information for interactions in all samples.}

\item{verbose}{a logical value indicating whether to print progress messages (default: FALSE).}
}
\value{
a list of data frames (numeric \link{0,Inf}) with additional columns for the maximum mean expression for ligand in sending cell type, and for the maximum active fraction for receptor in receiving cell type.
}
\description{
calculates max mean expression of each gene in all cell types in all samples.
}
