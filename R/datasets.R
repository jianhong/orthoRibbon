#' Human and Drosophila Ortholog Groups
#'
#' Ortholog groups for human (*Homo sapiens*) and fruit fly
#' (*Drosophila melanogaster*) identified using OrthoFinder.
#'
#' @format The dataset contains three objects:
#' \describe{
#'   \item{chrom_infos}{A list of data frames for
#'   \code{hsapiens} and \code{dmelanogaster}. Each data frame contains
#'   the columns \code{name}, \code{length}, \code{coord_system},
#'   \code{synonyms}, \code{toplevel}, \code{non_ref}, and \code{circular}.}
#'   \item{homologs_df}{A data frame containing the columns
#''   \code{gene_id}, \code{ortholog_group}, \code{symbol},
#'   \code{species}, \code{seq}, and \code{start}.}
#'   \item{genes_gr}{A \code{GRanges} object containing 40,528 ranges
#'   with two metadata columns: \code{gene_name} and \code{species}.}
#' }
#' @name human_fly
#' @docType data
#' @keywords datasets
#' @rdname ortholog_datasets
NULL


#' Ensembl Homologs for Zebrafish, Human, and Mouse
#'
#' Preprocessed homolog information for zebrafish, human, and mouse
#' retrieved from Ensembl. The original data were filtered and reduced
#' to minimize the dataset size.
#'
#' @rdname ortholog_datasets
#' @format The dataset contains three objects:
#' \describe{
#'   \item{chrom_infos}{A list of data frames for
#'   \code{hsapiens}, \code{mmusculus}, and \code{drerio}. Each data
#'   frame contains the columns \code{name}, \code{length},
#'   \code{coord_system}, \code{synonyms}, \code{toplevel},
#'   \code{non_ref}, and \code{circular}.}
#'   \item{homologs_df}{A data frame containing two columns:
#'   \code{gene_id1} and \code{gene_id2}.}
#'   \item{genes_gr}{A \code{GRanges} object containing 44,392 ranges
#'   with two metadata columns: \code{gene_name} and \code{species}.}
#' }
#'
#' @source Data retrieved from Ensembl using BioMart on March 13, 2026.
#' @name ensembl_hs_mm_dr
#' @docType data
#' @keywords datasets
NULL

#' @rdname ortholog_datasets
#' @name chrom_infos
NULL

#' @rdname ortholog_datasets
#' @name homologs_df
NULL

#' @rdname ortholog_datasets
#' @name genes_gr
NULL
