#' Add gene coordinates to homolog data
#' @description
#' Add chromosome, start information to the homolog data.
#' @param homolog_df A data frame with two columns which contains the homolog
#' pairs.
#' @param genes_gr A GRanges object for genes. It must contain the information
#' for all the homolog ids.
#' @param symbol_colunm_name The column name of gene symbol in 'genes_gr'.
#' @return A data frame with annotated homologs.
#' @importFrom S4Vectors mcols
add_coord <- function(homolog_df, genes_gr, symbol_colunm_name){
  stopifnot(is.matrix(homolog_df) || is.data.frame(homolog_df))
  stopifnot(ncol(homolog_df)>=2)
  stopifnot(is(genes_gr, 'GRanges'))
  stopifnot(is.character(symbol_colunm_name))
  stopifnot("'symbol_colunm_name' must be metadata column name of 'genes_gr'"=
              symbol_colunm_name %in% colnames(mcols(genes_gr)))
  stopifnot(all(homolog_df[, 1] %in% names(genes_gr)))
  stopifnot(all(homolog_df[, 2] %in% names(genes_gr)))
  ## make sure the first two columns of homolog_df be gene_id1 and gene_id2
  homolog_df <- as.data.frame(homolog_df)
  colnames(homolog_df)[c(1, 2)] <- c('gene_id1', 'gene_id2')

  homolog_df$symbol1 <-
    mcols(genes_gr[homolog_df$gene_id1])[[symbol_colunm_name]]
  homolog_df$symbol2 <-
    mcols(genes_gr[homolog_df$gene_id2])[[symbol_colunm_name]]

  homolog_df$seq1<- as.character(seqnames(genes_gr[homolog_df$gene_id1]))
  homolog_df$start1 <- start(genes_gr[homolog_df$gene_id1])

  homolog_df$seq2<- as.character(seqnames(genes_gr[homolog_df$gene_id2]))
  homolog_df$start2 <- start(genes_gr[homolog_df$gene_id2])

  return(homolog_df)
}

#' Retrieve all the ensembl gene IDs by given species abbreviations
#' @param com_name species abbreviations eg. "hsapiens", "mmusculus", "drerio"
#' @param marts A named list with the Mart object for each species
#' @return A list of character with ensembl gene IDs.
#' @export
#' @importFrom geneClusterPattern guessSpecies
#' @importFrom biomaRt getBM
#' @importFrom methods is
#' @importFrom AnnotationDbi keys
getGeneIDs <- function(com_name, marts){
  stopifnot(is.character(com_name))
  if(missing(marts)){
    full_name <- guessSpecies(com_name, output = "scientific name")
    if(requireNamespace("ChIPpeakAnno")){
      orgList <- lapply(full_name, ChIPpeakAnno::egOrgMap)
      names(orgList) <- com_name
      ids <- lapply(orgList, function(org){
        if(requireNamespace(org)){
          org <- get(org)
          ens_ids <- keys(org, keytype = "ENSEMBL")
        }else{
          stop(org, ' is required. Please install via BiocManager::install("',
               org, '")')
        }
      })
    }else{
      stop('"ChIPpeakAnno" is required.
           Please install via BiocManager::install("ChIPpeakAnno")')
    }

  }else{
    stopifnot(identical(com_name, names(marts)))
    null <- lapply(marts, function(.ele) {
      stopifnot("marts should be a list of Mart objects" =
                  is(.ele, 'Mart'))
    })
    ids <- lapply(marts, function(mart){
      getBM(attributes='ensembl_gene_id', mart=mart)
    })
  }
  return(ids)
}

#' Filter chromosomes
#' @param chrInfo A data.frame retrived from ensembl by function
#' getChromInfoFromEnsembl.
#' @param sp_min_chr_size The minimal size for a chromosome.
filterChrom <- function(chrInfo, sp_min_chr_size=1e7){
  stopifnot(is.data.frame(chrInfo))
  stopifnot('chromosome info must have columns "name" and "length"'=
              all(c('name', 'length') %in% colnames(chrInfo)))
  stopifnot(is.numeric(sp_min_chr_size))
  if('toplevel' %in% colnames(chrInfo)){
    chrInfo <- chrInfo[chrInfo$toplevel, , drop=FALSE]
  }else{
    warning('toplevel is not a colname in chromosome info.')
  }
  if('coord_system' %in% colnames(chrInfo)){
    if('chromosome' %in% chrInfo$coord_system){
      chrInfo <- chrInfo[chrInfo$coord_system %in% 'chromosome', , drop=FALSE]
    }else{
      warning('chromosome is not a selection in coord_system.')
    }
  }else{
    warning('coord_system is not a colname in chromosome info.')
  }
  if('length' %in% colnames(chrInfo)){
    chrInfo <- chrInfo[chrInfo$length>=sp_min_chr_size, , drop=FALSE]
  }else{
    stop('length is not a colname in chromosome info.')
  }
  return(chrInfo)
}

#' Retrieve the homolog pairs
#' @param ids A named list with the gene ids for each species
#' @param com_name species abbreviations eg. "hsapiens", "mmusculus", "drerio"
#' @param marts A named list with the Mart object for each species
#' @return A list with homolog GRanges.
#' @export
#' @importFrom geneClusterPattern getHomologGeneList
getHomologGRs <- function(ids, com_name, marts){
  stopifnot(length(names(ids))==length(ids))
  stopifnot(all(names(ids) %in% com_name))
  stopifnot(identical(names(ids), names(marts)))
  null <- lapply(marts, function(.ele) {
    stopifnot("marts should be a list of Mart objects" =
                is(.ele, 'Mart'))
  })
  target_species <- lapply(names(ids), function(name){
    com_name[com_name!=name]
  })
  homologs <- mapply(getHomologGeneList,
                     target_species, marts, ids,
                     SIMPLIFY = FALSE)
  names(homologs) <- names(ids)
  return(homologs)
}
#' Retrieve the homolog pairs
#' @param homologs A named list with the homolog GRanges. The output from
#' function \link{getHomologGRs}.
#' @return A data.frame with homologs.
#' @export
#' @importFrom S4Vectors mcols
#' @importFrom utils combn
getHomologIDs <- function(homologs){
  stopifnot(is.list(homologs))
  # merge the homolog ids
  homolog_ids <- lapply(homologs, function(hl){
    mc <- lapply(hl, mcols)
    mc <- lapply(mc, FUN=function(.ele){
      .ele$ensembl_gene_ids <- rownames(.ele)
      .ele[, c('ensembl_gene_ids', 'homolog_ensembl_gene_ids')]
    })
    ## bridged homologs
    cmb <- combn(seq_along(mc), 2, simplify = FALSE)
    mc_pair <- lapply(cmb, function(.ele){
      ## only keep the bridged homologs
      merge(mc[[.ele[1]]], mc[[.ele[2]]], by='ensembl_gene_ids')[, c(2, 3)]
    })
    ## fix the column names
    mc <- lapply(c(mc, mc_pair), function(.ele){
      colnames(.ele) <- c('x', 'y')
      .ele
    })
    do.call(rbind, mc)
  })

  homolog_ids <- do.call(rbind, homolog_ids)
  homolog_ids <- apply(homolog_ids, 1, sort)
  homolog_ids <- unique(t(homolog_ids))
  homolog_df <- as.data.frame(homolog_ids)
  colnames(homolog_df) <- c('gene_id1', 'gene_id2')
  return(homolog_df)
}

#' Retrieve all the gene ranges
#' @param ids A named list with the gene ids for each species
#' @param marts A named list with the Mart object for each species
#' @param homologs A list of homologs. It must be the output from function
#' \link[orthoRibbon:getHomologGRs]{getHomologGRs}l
#' @return A list with gene information GRanges.
#' @export
#' @importFrom geneClusterPattern grangesFromEnsemblIDs
#' @importFrom GenomicRanges GRangesList
getGeneGRs <- function(ids, marts, homologs){
  stopifnot(length(names(ids))==length(ids))
  stopifnot(identical(names(ids), names(marts)))
  null <- lapply(marts, function(.ele) {
    stopifnot("marts should be a list of Mart objects" =
                is(.ele, 'Mart'))
  })
  genes <- mapply(grangesFromEnsemblIDs, marts, ids, SIMPLIFY = FALSE)
  genes <- mapply(genes, names(genes), FUN=function(.ele, .species){
    .ele$species <- .species
    .ele
  })
  genes_gr <- lapply(homologs, function(hls){
    hls <- mapply(hls, names(hls), FUN=function(hl, species){
      names(hl) <- hl$homolog_ensembl_gene_ids
      hl$species <- species
      hl$homolog_ensembl_gene_ids <- NULL
      hl
    }, SIMPLIFY = FALSE)
    unlist(GRangesList(hls), use.names = FALSE)
  })
  genes_gr <- unlist(GRangesList(c(genes_gr, genes)), use.names=FALSE)
}

#' Add gene information to homolog data.frame
#' @param homolog_df The data.frame with homologs. The output from function
#' \link{getHomologIDs}
#' @param genes_gr The GRanges object with gene info. The output from function
#' \link{getGeneGRs}
#' @importFrom methods is
#' @importFrom Seqinfo seqnames
#' @importFrom BiocGenerics start
addGeneInfo <- function(homolog_df, genes_gr){
  stopifnot(is.data.frame(homolog_df))
  stopifnot(all(c('gene_id1', 'gene_id2') %in% colnames(homolog_df)))
  stopifnot(is(genes_gr, 'GRanges'))
  stopifnot(all(c('gene_name', 'species') %in% colnames(mcols(genes_gr))))

  homolog_df$symbol1 <- genes_gr[homolog_df$gene_id1]$gene_name
  homolog_df$symbol2 <- genes_gr[homolog_df$gene_id2]$gene_name
  homolog_df$species1 <- genes_gr[homolog_df$gene_id1]$species
  homolog_df$seq1<- as.character(seqnames(genes_gr[homolog_df$gene_id1]))
  homolog_df$start1 <- start(genes_gr[homolog_df$gene_id1])
  homolog_df$species2 <- genes_gr[homolog_df$gene_id2]$species
  homolog_df$seq2<- as.character(seqnames(genes_gr[homolog_df$gene_id2]))
  homolog_df$start2 <- start(genes_gr[homolog_df$gene_id2])

  return(homolog_df)
}

#' Filter the homolog data.frame by chromosome names
#' @param homolog_df The data.frame with homologs. The output from function
#' \link{addGeneInfo}
#' @param chrom_infos A list with the chromosome information data.frame.
#' @param max_links A numeric. The maiximal link number to show in the plot.
#' @return A data.frame with filtered homologs
subsetHomologsByChrom <- function(homolog_df, chrom_infos, max_links=5000){
  stopifnot(is.list(chrom_infos))
  stopifnot(is.data.frame(homolog_df))
  stopifnot(all(c('species1', 'seq1', 'species2', 'seq2') %in%
                  colnames(homolog_df)))
  used_seqs <- lapply(chrom_infos, function(.ele){
    sortSeqlevels(.ele$name)
  })

  used_seqs_ <- mapply(paste, rep(names(chrom_infos), lengths(used_seqs)),
                       unlist(used_seqs))

  ## subset the homolog_df by the used chromosomes
  homolog_df$chr_sp1 <- paste(homolog_df$species1, homolog_df$seq1)
  homolog_df$chr_sp2 <- paste(homolog_df$species2, homolog_df$seq2)
  keep1 <- homolog_df$chr_sp1 %in% used_seqs_
  keep2 <- homolog_df$chr_sp2 %in% used_seqs_
  homolog_df <- homolog_df[keep1 & keep2, , drop=FALSE]
  if(max_links<nrow(homolog_df)){
    homolog_df <- homolog_df[sample.int(nrow(homolog_df), size = max_links), ,
                             drop=FALSE]
  }
  return(homolog_df)
}

#' Get the best order of the chromosomes
#' @description
#' The best order of the chromosomes decided by the distance matrix.
#' The distance matrix among chromosomes are the count of the homologs
#' @param homolog_df The data.frame with homologs. The output from function
#' \link{subsetHomologsByChrom}
#' @param chrom_infos A list with the chromosome information data.frame.
#' @param method a character string with the name of the seriation method or
#'  'chr'. If using 'chr', it will try to find the best order by the chromosome
#'  name. Otherwise see \link[seriation:seriate]{seriate}. It will be worth to
#'  try 'TSP' first.
#' @return A list with the ordered chromosome names
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom seriation get_order seriate
#' @importFrom stats hclust as.dist dist cor
getChrOrders <- function(homolog_df, chrom_infos, method = 'TSP'){
  stopifnot(is.list(chrom_infos))
  stopifnot(is.data.frame(homolog_df))
  stopifnot(all(c('chr_sp1', 'chr_sp2') %in%
                  colnames(homolog_df)))
  chr_dist <- table(homolog_df[, c('chr_sp1', 'chr_sp2')])
  used_seqs <- lapply(chrom_infos, function(.ele){
    sortSeqlevels(.ele$name)
  })
  ## the mapping chain is the order of com_name
  sp_pairs <- lapply(seq_along(used_seqs)[-length(used_seqs)],
                     function(x) names(used_seqs)[c(x, x+1)])
  ## Finds an order that minimizes crossings
  chr_pairs <- lapply(sp_pairs,
                      FUN=function(p){
                        chr_a <- paste(p[1], used_seqs[[p[1]]])
                        chr_b <- paste(p[2], used_seqs[[p[2]]])
                        if(sum(chr_a %in% rownames(chr_dist)) >=
                           sum(chr_a %in% colnames(chr_dist)) &
                           sum(chr_b %in% colnames(chr_dist)) >=
                           sum(chr_b %in% rownames(chr_dist))){
                          cur_dist <-
                            chr_dist[chr_a[chr_a %in% rownames(chr_dist)],
                                     chr_b[chr_b %in% colnames(chr_dist)]]
                        }else{
                          cur_dist <-
                            t(chr_dist[chr_b[chr_b %in% rownames(chr_dist)],
                                       chr_a[chr_a %in% colnames(chr_dist)]])
                        }
                        # adjust by weight
                        rs <- rowSums(cur_dist)
                        cs <- colSums(cur_dist)
                        cur_dist <- cur_dist / (outer(rs, cs, "+"))
                        if(method=='chr'){
                          ## always set the second ordered by chromosome name
                          hc <- hclust(dist(cur_dist), method = "complete")
                          ords <- list(
                            sub(paste0(p[1], ' '), '',
                                rownames(cur_dist)[hc$order]),
                            sub(paste0(p[2], ' '), '', chr_b)
                          )
                        }else{
                          row_order <- seriate(
                            as.dist(1 - cor(t(cur_dist), method = "spearman")),
                            method = method)
                          col_order <- seriate(
                            as.dist(1 - cor(cur_dist, method = "spearman")),
                            method = method)
                          ords <- list(
                            sub(paste0(p[1], ' '), '',
                                rownames(cur_dist)[get_order(row_order)]),
                            sub(paste0(p[2], ' '), '',
                                colnames(cur_dist)[get_order(col_order)]))
                        }
                        names(ords) <- p
                        return(ords)
                      })
  ## sort the chromosome by previous one
  chr_orders <- chr_pairs[[1]]
  for(j in seq_along(chr_pairs)[-1]){
    cur_chr_pairs <- chr_pairs[[j]]
    chr_orders[[names(cur_chr_pairs)[2]]] <-
      cur_chr_pairs[[2]][match(chr_orders[[names(cur_chr_pairs)[1]]],
                               cur_chr_pairs[[1]])]
  }
  ## add missing chromosomes
  chr_orders <- mapply(chr_orders, used_seqs, FUN=function(.chr_ord, .all_chr){
    c(.chr_ord, .all_chr[!.all_chr %in% .chr_ord])
  }, SIMPLIFY = FALSE)
}

#' Build chromosome plotting dataframes
#' @param chrom_infos A list with the chromosome information data.frame.
#' @param chr_orders A list with the ordered chromosome names
#' @return A list with the ordered chromosome plot information
buildChromDF <- function(chrom_infos, chr_orders){
  stopifnot(identical(names(chrom_infos), names(chr_orders)))
  chrom_infos <- mapply(chrom_infos, chr_orders,
                          FUN=function(chr_size, chr_order){
                            chr_size[match(chr_order, chr_size$name), ,
                                     drop=FALSE]
                          }, SIMPLIFY = FALSE)


  chrom_df <- lapply(chrom_infos, function(.ele){
    stopifnot('chromosome info must have columns "name" and "length"'=
                all(c('name', 'length') %in% colnames(.ele)))
    chromdf <- data.frame(
      chrom = .ele$name,
      chromix = seq.int(nrow(.ele)),
      chrsize = .ele$length
    )
    num_chroms <- nrow(.ele)
    min_chr_size <- 10
    max_chr_size <- 50
    min_gap_percent <- 0.12
    max_gap_percent <- 0.5

    if (num_chroms < min_chr_size) {
      percent_chrom_as_spaces <- min_gap_percent
    } else if (num_chroms > max_chr_size) {
      percent_chrom_as_spaces <- max_gap_percent
    } else {
      percent_chrom_as_spaces <- min_gap_percent +
        (max_gap_percent - min_gap_percent) *
        (num_chroms - min_chr_size) / (max_chr_size - min_chr_size)
    }

    percent_chrom_as_chroms <- 1 - percent_chrom_as_spaces
    space_between_chrom <- percent_chrom_as_spaces / (num_chroms - 1)
    # Calculate plot positions (chromosome coordinates)
    total_chrom_len <- sum(chromdf$chrsize)
    chromdf$chrPlotPercent <- (chromdf$chrsize / total_chrom_len) *
      percent_chrom_as_chroms
    chromdf$chrPlotOffset <-
      c(0, cumsum(chromdf$chrPlotPercent)[-nrow(chromdf)]) +
      (space_between_chrom * chromdf$chromix)
    chromdf
  })
}

#' Build the homolog links dataframe
#' @param homolog_df The data.frame with homologs. The output from function
#' \link{subsetHomologsByChrom}
#' @param chrom_df A list with the ordered chromosome plot information. The
#' output of function
#' @return A list with the homolog links plot information.
#' @importFrom stats ave
buildHomologLinksDF <- function(homolog_df, chrom_df){
  used_seqs <- lapply(chrom_df, function(.ele){
    sort(.ele$name)
  })
  ## the mapping chain is the order of com_name
  sp_pairs <- lapply(seq_along(used_seqs)[-length(used_seqs)],
                     function(x) names(used_seqs)[c(x, x+1)])

  homolog_df_list <- split(homolog_df,
                           paste(homolog_df$species1,
                                 homolog_df$species2))
  homolog_df_list_pairs <- strsplit(names(homolog_df_list), ' ')

  homolog_df_list <- lapply(sp_pairs, function(.ele){
    id <- lapply(homolog_df_list_pairs, function(x) all(.ele %in% x))
    id <- which(unlist(id))
    if(length(id)<1) stop('no homologs found in ', .ele[1], ' and ', .ele[2])
    x <- homolog_df_list[[id[1]]]
    cn <- colnames(x)
    if(all(.ele==homolog_df_list_pairs[[id[1]]])){
      cn <- sub('1$', '_top', cn)
      cn <- sub('2$', '_bottom', cn)
    }else{
      # switch the 1 and 2 in column names
      cn <- sub('1$', '_bottom', cn)
      cn <- sub('2$', '_top', cn)
    }
    colnames(x) <- cn
    x$sp_top <- .ele[1]
    x$sp_bottom <- .ele[2]
    # Calculate positions for index plotting (top species)
    chromdf <- chrom_df[[.ele[1]]]
    x$top_ChromSize <- chromdf[match(x$seq_top, chromdf$chrom), 'chrsize']
    x$top_ChromPercent <-
      chromdf[match(x$seq_top, chromdf$chrom), 'chrPlotPercent']
    x$top_ChromOffset <-
      chromdf[match(x$seq_top, chromdf$chrom), 'chrPlotOffset']
    x$topIx <- ave(x$start_top, x$seq_top, FUN = function(i)
      rank(i, ties.method = "first")) ##
    x$topIx_Size <- ave(x$start_top, x$seq_top, FUN = length)
    x$topIx_geneOffset <- (x$topIx / x$topIx_Size) * x$top_ChromPercent
    x$topIx_finalOffset <- x$top_ChromOffset + x$topIx_geneOffset
    # Calculate positions for chromosome plotting (top species)
    x$topChr_geneOffset <- (x$start_top / x$top_ChromSize) * x$top_ChromPercent
    x$topChr_finalOffset <- x$top_ChromOffset + x$topChr_geneOffset

    # Calculate positions for index plotting (bottom species)
    chromdf <- chrom_df[[.ele[2]]]
    x$bottom_ChromSize <-
      chromdf[match(x$seq_bottom, chromdf$chrom), 'chrsize']
    x$bottom_ChromPercent <-
      chromdf[match(x$seq_bottom, chromdf$chrom), 'chrPlotPercent']
    x$bottom_ChromOffset <-
      chromdf[match(x$seq_bottom, chromdf$chrom), 'chrPlotOffset']
    x$bottomIx <- ave(x$start_bottom, x$seq_bottom,
                      FUN = function(i) rank(i, ties.method = "first"))
    x$bottomIx_Size <- ave(x$start_bottom, x$seq_bottom, FUN = length)
    x$bottomIx_geneOffset <-
      (x$bottomIx / x$bottomIx_Size) * x$bottom_ChromPercent
    x$bottomIx_finalOffset <- x$bottom_ChromOffset + x$bottomIx_geneOffset
    # Calculate positions for chromosome plotting (bottom species)
    x$bottomChr_geneOffset <-
      (x$start_bottom / x$bottom_ChromSize) * x$bottom_ChromPercent
    x$bottomChr_finalOffset <- x$bottom_ChromOffset + x$bottomChr_geneOffset

    # Calculate positions for minimimal plotting
    # The idea is that sort all the seq_top and seq_bottom to make the plot
    # by blocks
    x$topMini <- as.integer(ave(x$seq_bottom, x$seq_top, FUN = function(i){
      i <- as.numeric(factor(i, levels = chrom_df[[.ele[2]]]$chrom))
      rank(i, ties.method = 'first')
    }))
    x$topMini_Size <- as.integer(ave(x$seq_bottom, x$seq_top, FUN = length))
    x$topMini_geneOffset <-
      (x$topMini/x$topMini_Size) * x$top_ChromPercent
    x$topMini_finalOffset <- x$top_ChromOffset + x$topMini_geneOffset
    x$bottomMini <- as.integer(ave(x$seq_top, x$seq_bottom, FUN = function(i){
      i <- as.numeric(factor(i, levels = chrom_df[[.ele[1]]]$chrom))
      rank(i, ties.method = 'first')
    }))
    x$bottomMini_Size <- as.integer(ave(x$seq_top, x$seq_bottom, FUN = length))
    x$bottomMini_geneOffset <-
      (x$bottomMini/x$bottomMini_Size) * x$bottom_ChromPercent
    x$bottomMini_finalOffset <- x$bottom_ChromOffset + x$bottomMini_geneOffset

    x
  })
}

create_bezier_points <- function(finalOffset1, finalOffset2, i,
                                 col1, col2, cl1=4){
  stopifnot(cl1<=4&&cl1>=0)
  data.frame(x = seq(finalOffset1, finalOffset2, length.out=4),
             y = c(i-1, i, i-1, i),
             col = rep(c(col1, col2), c(cl1, 4-cl1)))
}

create_bezier_matrix <- function(data, i, colname1, colname2,
                                 col1='seq_top', col2='seq_bottom', cl1=4){
  x <- apply(data[, c(colname1, colname2, col1, col2)],
             1, function(.ele){
               create_bezier_points(.ele[1], .ele[2], i,
                                    .ele[3], .ele[4],
                                    cl1=cl1)
             }, simplify = FALSE)
  id <- rep(paste(i, seq_along(x), sep = '_'), each=4)
  x <- do.call(rbind, x)
  x$id <- id
  x
}

#' Build the data.frame for Bezier curve
#' @param homolog_df_list A list with the data.frame of plot data for homologs.
#' @param colname1,colname2 The column names for top and bottom final positions.
#' @param col1,col2 The column names for color.
#' @param cl1 A nummeric. The Bezier curve is created with two control points.
#' And cl1 should be no more than 4 and no less than 0.
#' If it is set to 4, all points in Bezier curve will be set to col1.
#' If it is set to 3, the top point, two control points will be set to col1 and
#' the bottom point will be set to col2.
#' If it is set to 2, the top point, top control point will be set to col1 and
#' the bottom control point, bottom point will be set to col2.
#' If it is set to 1, the top point will be set to col1 and others col2.
#' If it is set to 0, all points will be set to col2.
#' @export
#' @examples
#' # example code
#'
buildBezierDF <- function(homolog_df_list, colname1, colname2,
                          col1='seq_top', col2='seq_bottom',
                          cl1=4){
  bezier_df <- lapply(seq_along(homolog_df_list), function(i){
    create_bezier_matrix(homolog_df_list[[i]],
                         i, colname1, colname2, col1, col2, cl1=cl1)
  })
  bezier_df <- do.call(rbind, bezier_df)
}

buildChromBarDF <- function(chrom_df){
  chrom_bars <- mapply(chrom_df, seq_along(chrom_df)-1, names(chrom_df),
                       FUN=function(chromdf, i, sp){
                         data.frame(x=c(chromdf$chrPlotOffset,
                                        chromdf$chrPlotOffset +
                                          chromdf$chrPlotPercent),
                                    y=i,
                                    chrom = chromdf$chrom,
                                    label = chromdf$chrom,
                                    sp = sp)
                       }, SIMPLIFY = FALSE)
  chrom_bars_df <- do.call(rbind, chrom_bars)
  return(chrom_bars_df)
}

buildSymbolList <- function(homolog_df, y, topX=TRUE){
  if(isTRUE(topX)){
    xChr <- homolog_df$topChr_finalOffset
    xIx <- homolog_df$topIx_finalOffset
    xMini <- homolog_df$topMini_finalOffset
  }  else{
    xChr <- homolog_df$bottomChr_finalOffset
    xIx <- homolog_df$bottomIx_finalOffset
    xMini <- homolog_df$bottomMini_finalOffset
  }
  data.frame(xChr = xChr,
             xIx = xIx,
             xMini = xMini,
             y=y,
             top_label=homolog_df$symbol_top,
             bottom_label=homolog_df$symbol_bottom,
             seq_top=homolog_df$seq_top,
             seq_bottom=homolog_df$seq_bottom)
}
#' @importFrom dplyr %>% group_by mutate ungroup
#' @importFrom rlang .data
buildChromLabelDF <- function(chrom_bars_df){
  chrom_label_df <- chrom_bars_df %>% group_by(.data$sp, .data$chrom) %>%
    mutate(
      x = mean(.data$x, na.rm = TRUE),
      y = mean(.data$y, na.rm = TRUE)
    ) %>%
    ungroup() %>% unique()
  return(chrom_label_df)
}
