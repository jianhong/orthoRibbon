#' Build the plot data
#' @description
#' Build the plot data from homolog data.frame, gene annotations and chromosome
#' informations.
#' @param com_name Species abbreviations eg. "hsapiens", "mmusculus", "drerio"
#' @param homolog_df The data.frame for homologs with column names
#' "gene_id1", "gene_id2".
#' @param genes_gr A GRanges object for genes. It must contain the information
#' for all the homolog ids.
#' 'gene_name', 'species' must be metadata column name of 'genes_gr'.
#' @param chrom_infos A list of chromosome information.
#' The Chromosome information for each species must be a data.frame with columns
#' "name" and "length".
#' @param sp_min_chr_size The minimal chromosome size.
#' @param chr_orders Chromosome orders for plot.
#' @param filterByCoordSystem Filter the chromosomes or not.
#' @param max_links A numeric. The maiximal link number to show in the plot.
#' @param chromosome_order_method a character string with the name of the
#' seriation method or 'chr'.
#'  If using 'chr', it will try to find the best order by the chromosome
#'  name. Otherwise see \link[seriation:seriate]{seriate}. It will be worth to
#'  try 'TSP' first.
#' @return A list with elements "homolog_df_list", "chrom_bars_df",
#' "chrom_label_df" , "symbol_list_top" and, "symbol_list_bottom" for plot.
#' @export
#' @examples
#' # example code
#'
buildPlotData <- function(com_name, homolog_df, genes_gr, chrom_infos,
                          sp_min_chr_size=10000000, chr_orders=NULL,
                          filterByCoordSystem = TRUE, max_links=Inf,
                          chromosome_order_method='chr'){
  # Step1 check input
  stopifnot(is.data.frame(homolog_df))
  df_from_flag <- NA
  if(all(c('gene_id', 'ortholog_group') %in% colnames(homolog_df))){
    if(!all(c('species', 'seq', 'start') %in% colnames(homolog_df))){
      stopifnot(is(genes_gr, 'GRanges'))
      stopifnot("'gene_name', 'species' must be metadata column name of 'genes_gr'"=
                  all(c('gene_name', 'species') %in% colnames(mcols(genes_gr))))
      df_from_flag <- 'ortholog_group_only'
    }else{
      df_from_flag <- 'ortholog_group_with_gene_info'
    }
  }else{
    stopifnot(
      'The input homolog_df must be a data.frame with column names
    "gene_id1", "gene_id2"'=
        all(c('gene_id1', 'gene_id2') %in% colnames(homolog_df)))
    if(!all(c('symbol1', 'symbol2',
              'species1', 'seq1', 'start1',
              'species2', 'seq2', 'start2') %in% colnames(homolog_df))){
      stopifnot(is(genes_gr, 'GRanges'))
      stopifnot("'gene_name', 'species' must be metadata column name of 'genes_gr'"=
                  all(c('gene_name', 'species') %in% colnames(mcols(genes_gr))))
      df_from_flag <- 'ortholog_pair_only'
    }else{
      df_from_flag <- 'ortholog_pair_with_gene_info'
    }
  }

  null <- lapply(chrom_infos, function(chrInfo){
    stopifnot('chromosome info must have columns "name" and "length"'=
                all(c('name', 'length') %in% colnames(chrInfo)))
  })
  stopifnot('The names of "com_name" and "chrom_infos" must be identical'=
              identical(names(com_name), names(chrom_infos)))
  stopifnot(is.numeric(sp_min_chr_size))
  stopifnot(is.character(chromosome_order_method))

  # Step2 add gene information to homolog data.frame
  if(isTRUE(df_from_flag!='ortholog_pair_with_gene_info')){
    homolog_df <- addGeneInfo(homolog_df = homolog_df,
                              genes_gr = genes_gr,
                              type = df_from_flag)
  }

  # Step3 filter the homolog data.frame by chromosome names
  if(isTRUE(filterByCoordSystem)){
    chrom_infos <- lapply(chrom_infos, filterChrom,
                          sp_min_chr_size=sp_min_chr_size)
  }
  homolog_df <- subsetHomologsByChrom(homolog_df, chrom_infos,
                                      max_links=max_links)

  # Step4 order the chromosome names if user does not provided
  if(!all(com_name %in% names(chr_orders))){
    chr_orders <- getChrOrders(homolog_df, chrom_infos,
                               method = chromosome_order_method)
  }

  # Step5 Build chromosome plotting dataframes
  chrom_df <- buildChromDF(chrom_infos, chr_orders)

  # Step6 Build homolog links plotting dataframes
  homolog_df_list <- buildHomologLinksDF(homolog_df, chrom_df)

  # Step7 Create chromosome bar data
  chrom_bars_df <- buildChromBarDF(chrom_df)
  chrom_label_df <- buildChromLabelDF(chrom_bars_df)

  # Step8 # create symbol labels
  symbol_list_top <- buildSymbolList(homolog_df_list[[1]], y=0)
  symbol_list_bottom <-
    buildSymbolList(homolog_df_list[[length(homolog_df_list)]],
                    y=length(homolog_df_list), topX=FALSE)

  return(list(
    homolog_df_list = homolog_df_list,
    chrom_bars_df = chrom_bars_df,
    chrom_label_df = chrom_label_df,
    symbol_list_top = symbol_list_top,
    symbol_list_bottom = symbol_list_bottom
  ))
}

