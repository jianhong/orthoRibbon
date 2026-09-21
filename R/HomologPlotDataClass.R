#' HomologPlotData S4 class
#'
#' Wraps the collection of data frames/lists used for an orthoRibbon plot into
#' a single S4 object, with list-like `$`, `[[`, and `[` access.
#'
#' @slot homolog_df_list list of data.frame with top and bottom plot information.
#' @slot chrom_bars_df data.frame of chromosome plot information.
#' @slot chrom_label_df data.frame of chromosome labels information.
#' @slot symbol_df_top,symbol_df_bottom data.frame of symbols for top and
#' bottom layer.
#' @slot common_name species abbreviations e.g. "hsapiens", "mmusculus", "drerio"
#' @param x A HomologPlotData object.
setClass(
  "HomologPlotData",
  representation(
    homolog_df_list = "list",
    chrom_bars_df = "data.frame",
    chrom_label_df = "data.frame",
    symbol_df_top = "data.frame",
    symbol_df_bottom = "data.frame",
    common_name = "character"
  )
)

#' @rdname HomologPlotData-class
#' @param homolog_df_list list of data.frame with top and bottom plot information.
#' @param chrom_bars_df data.frame of chromosome plot information.
#' @param chrom_label_df data.frame of chromosome labels information.
#' @param symbol_df_top,symbol_df_bottom data.frame of symbols for top and
#' bottom layer.
#' @param common_name species abbreviations e.g. "hsapiens", "mmusculus", "drerio"
#' @importFrom methods new
#' @export
HomologPlotData <- function(homolog_df_list,
                            chrom_bars_df,
                            chrom_label_df,
                            symbol_df_top,
                            symbol_df_bottom,
                            common_name) {
  new(
    "HomologPlotData",
    homolog_df_list    = homolog_df_list,
    chrom_bars_df       = chrom_bars_df,
    chrom_label_df       = chrom_label_df,
    symbol_df_top      = symbol_df_top,
    symbol_df_bottom   = symbol_df_bottom,
    common_name = common_name
  )
}

setValidity("HomologPlotData", function(object) {
  errors <- character()
  if (!is.list(object@homolog_df_list)) {
    errors <- c(errors, "homolog_df_list must be a list")
  }else{
    if(length(object@homolog_df_list)>0){
      ## check the columns required.
      isDF <- vapply(object@homolog_df_list, is.data.frame, logical(1L))
      if(!all(isDF)){
        errors <- c(errors, "homolog_df_list must be a list of data.frame")
      }else{
        cn_errors <- lapply(object@homolog_df_list, function(.ele){
          if(!all(c("gene_id_top", "gene_id_bottom", "symbol_top", "symbol_bottom", "species_top", "seq_top", "start_top", "species_bottom", "seq_bottom", "start_bottom", "chr_sp_top", "chr_sp_bottom", "sp_top", "sp_bottom", "top_ChromSize", "top_ChromPercent", "top_ChromOffset", "topIx", "topIx_Size", "topIx_geneOffset", "topIx_finalOffset", "topChr_geneOffset", "topChr_finalOffset", "bottom_ChromSize", "bottom_ChromPercent", "bottom_ChromOffset", "bottomIx", "bottomIx_Size", "bottomIx_geneOffset", "bottomIx_finalOffset", "bottomChr_geneOffset", "bottomChr_finalOffset", "topMini", "topMini_Size", "topMini_geneOffset", "topMini_finalOffset", "bottomMini", "bottomMini_Size", "bottomMini_geneOffset", "bottomMini_finalOffset") %in% colnames(.ele))){
            "The dataframe in homolog_df_list must have columns 'gene_id_top', 'gene_id_bottom', 'symbol_top', 'symbol_bottom', 'species_top', 'seq_top', 'start_top', 'species_bottom', 'seq_bottom', 'start_bottom', 'chr_sp_top', 'chr_sp_bottom', 'sp_top', 'sp_bottom', 'top_ChromSize', 'top_ChromPercent', 'top_ChromOffset', 'topIx', 'topIx_Size', 'topIx_geneOffset', 'topIx_finalOffset', 'topChr_geneOffset', 'topChr_finalOffset', 'bottom_ChromSize', 'bottom_ChromPercent', 'bottom_ChromOffset', 'bottomIx', 'bottomIx_Size', 'bottomIx_geneOffset', 'bottomIx_finalOffset', 'bottomChr_geneOffset', 'bottomChr_finalOffset', 'topMini', 'topMini_Size', 'topMini_geneOffset', 'topMini_finalOffset', 'bottomMini', 'bottomMini_Size', 'bottomMini_geneOffset', and 'bottomMini_finalOffset'"
          }else{
            NULL
          }
        })
        cn_errors <- unique(unlist(cn_errors))
        errors <- c(errors, cn_errors)
      }
    }
  }
  if(length(object@common_name)){
    if(length(object@common_name)<2){
      errors <- c(errors, "common_name must be more than one element")
    }
  }
  if (!is.data.frame(object@chrom_bars_df)) {
    errors <- c(errors, "chrom_bars_df must be a data.frame")
  }
  if(length(object@chrom_bars_df)){
    if(!all(c('x', 'y', 'sp', 'chrom', 'label') %in% colnames(object@chrom_bars_df))){
      errors <-
        c(errors,
          "chrom_bars_df must have column names 'x', 'y', 'sp', 'label' and 'chrom'")
    }
    if(!length(object@common_name)){
      errors <-
        c(errors, "common_name is missing")
    }
    if(!all(object@common_name %in% object@chrom_bars_df$sp)){
      errors <- c(
        errors, "common_name must be elements in chrom_bars_df$sp"
      )
    }
  }
  if (!is.data.frame(object@chrom_label_df)) {
    errors <- c(errors, "chrom_label_df must be a data.frame")
  }
  if(length(object@chrom_label_df)){
    if(!all(c('x', 'y', 'sp', 'label', 'chrom') %in% colnames(object@chrom_label_df))){
      errors <-
        c(errors,
          "chrom_bars_df must have column names 'x', 'y', 'sp', 'chrom', and 'label'")
    }
  }

  if (length(errors) == 0) TRUE else errors
})

#' @rdname HomologPlotData-class
#' @param name,i slot name of HomologPlotData.
#' @method $ HomologPlotData
#' @aliases $,HomologPlotData-method
#' @importFrom methods slot
#' @export
setMethod("$", "HomologPlotData", function(x, name) {
  slot(x, name)
})

#' @rdname HomologPlotData-class
#' @param value value to be set
#' @method $<- HomologPlotData
#' @aliases $<-,HomologPlotData-method
#' @importFrom methods "slot<-"
#' @export
setMethod("$<-", "HomologPlotData", function(x, name, value){
  slot(x, name, check=TRUE) <- value
  x
})
#' @rdname HomologPlotData-class
#' @method [[ HomologPlotData
#' @aliases "[[",HomologPlotData,ANY,ANY-method
#' @param j,drop,... not used.
#' @export
setMethod("[[", "HomologPlotData", function(x, i, ...) {
  slot(x, i)
})

#' @rdname HomologPlotData-class
#' @method [ HomologPlotData
#' @aliases "[",HomologPlotData,ANY,ANY,ANY-method
#' @importFrom methods slotNames
#' @export
setMethod("[", "HomologPlotData", function(x, i, ...) {
  nms <- slotNames(x)
  if (is.character(i)) {
    keep <- i
  } else {
    keep <- nms[i]
  }
  if (!all(keep %in% nms)) {
    stop("Unknown slot(s): ", paste(setdiff(keep, nms), collapse = ", "))
  }
  stats::setNames(lapply(keep, function(nm) slot(x, nm)), keep)
})

#' @rdname HomologPlotData-class
#' @method names HomologPlotData
#' @aliases "names",HomologPlotData-method
#' @export
setMethod("names", "HomologPlotData", function(x) {
  slotNames(x)
})

#' @rdname HomologPlotData-class
#' @importFrom utils .DollarNames
#' @method .DollarNames HomologPlotData
#' @aliases .DollarNames,HomologPlotData,character-method
#' @param pattern A regular expression. Only matching names are returned.
#' @export
.DollarNames.HomologPlotData <- function(x, pattern = "") {
  grep(pattern, slotNames(x), value = TRUE)
}

#' @rdname HomologPlotData-class
#' @method subset HomologPlotData
#' @aliases subset,HomologPlotData,list-method
#' @param subset is a named list, e.g.
#'   list(speciesA = c("chr1", "chr2"), speciesB = c("chr1", "chrX"))
#' Names must match values in chrom_bars_df$sp; each element's values must
#' match either chrom_bars_df$chrom or chrom_bars_df$label for that species.
#' The first name is treated as the "top" species, the last as the "bottom"
#' species, and homolog_df_list / symbol_df_top / symbol_df_bottom are
#' filtered to match.
#' @export
setMethod("subset", "HomologPlotData", function(x, subset, ...){
  if(!missing(subset)){
    stopifnot('subset must be a named list'=is.list(subset))
    stopifnot('subset must be a named list'=length(names(subset))==length(subset))
    n <- names(subset)
    if(!all(n %in% x$chrom_bars_df$sp)){
      stop('not all the names of subset exist in the plotData as species')
    }
    null <- lapply(subset, function(i){
      stopifnot('the elements of subset must be labels or chrom exists in chrom_bars_df' =
                  all(i %in% x$chrom_bars_df$chrom | i %in% x$chrom_bars_df$label))
    })

    subset <- data.frame(sp=rep(names(subset), lengths(subset)),
                         values=unlist(subset, use.names=FALSE))
    subsetN <- paste(subset$sp, subset$values)
    keep <- (!x$chrom_bars_df$sp %in% subset$sp) |
      (paste(x$chrom_bars_df$sp,x$chrom_bars_df$label) %in% subsetN) |
      (paste(x$chrom_bars_df$sp,x$chrom_bars_df$chrom) %in% subsetN)
    x$chrom_bars_df <- x$chrom_bars_df[keep, , drop=FALSE]
    keep <- (!x$chrom_label_df$sp %in% subset$sp) |
      (paste(x$chrom_label_df$sp,x$chrom_label_df$label) %in% subsetN) |
      (paste(x$chrom_label_df$sp,x$chrom_label_df$chrom) %in% subsetN)
    x$chrom_label_df <- x$chrom_label_df[keep, , drop=FALSE]
    for(sp in x$common_name){
      seqs <- unique(x$chrom_bars_df$chrom[x$chrom_bars_df$sp==sp])
      for(j in seq_along(x$homolog_df_list)){
        if(isTRUE(length(x$homolog_df_list[[j]]) &&
                  sp %in% subset$sp)){
          if(sp %in% x$homolog_df_list[[j]]$species_top){
            target_column <- '_top'
          }else if(sp %in% x$homolog_df_list[[j]]$species_bottom){
            target_column <- '_bottom'
          }else{
            next
          }
          tmp <- x$homolog_df_list[[j]]
          tmp <- tmp[tmp[[paste0('species', target_column)]]!= sp |
                       (tmp[[paste0('species', target_column)]] == sp &
                          tmp[[paste0('seq', target_column)]] %in% seqs),
                     , drop=FALSE]
          x$homolog_df_list[[j]] <- tmp
        }
      }
      filter_symbol_df <- function(sym_df, sps, seqs){
        if(all(sps %in% subset$sp)){
          ## both need to filter
          sym_df <- sym_df[sym_df$seq_top %in% seqs[[1]] &
                       sym_df$seq_bottom %in% seqs[[2]],
                     , drop=FALSE]
        }else{
          if(sps[1] %in% subset$sp){
            sym_df <- sym_df[sym_df$seq_top %in% seqs[[1]],
                       , drop=FALSE]
          }else{
            sym_df <- sym_df[sym_df$seq_bottom %in% seqs[[2]],
                       , drop=FALSE]
          }
        }
        sym_df
      }
      if(sp==x$common_name[1]){ # topSpeices
        spTop <- x$common_name[c(1, 2)]
        seqTop <- lapply(spTop, function(sp){
          x$chrom_bars_df[x$chrom_bars_df$sp==sp, 'chrom', drop=TRUE]
        })
        if(isTRUE(length(x$symbol_df_top) &&
                  any(spTop %in% subset$sp))){
          x$symbol_df_top <- filter_symbol_df(x$symbol_df_top, spTop, seqTop)
        }
      }else if(sp==x$common_name[length(x$common_name)]){# bottomSpecies
        spBottom <- x$common_name[c(length(x$common_name)-1,
                                    length(x$common_name))]
        seqBottom <- lapply(spBottom, function(sp){
          x$chrom_bars_df[x$chrom_bars_df$sp==sp, 'chrom', drop=TRUE]
        })
        if(isTRUE(length(x$symbol_df_bottom) &&
                  any(spBottom %in% subset$sp))){
          x$symbol_df_bottom <- filter_symbol_df(x$symbol_df_bottom,
                                                 spBottom, seqBottom)
        }
      }
    }
  }
  return(x)
})

setMethod("show", "HomologPlotData", function(object) {
  cat("<HomologPlotData>\n")
  for (nm in slotNames(object)) {
    val <- slot(object, nm)
    cat(sprintf(
      "  %-20s: %s [%s]\n",
      nm,
      class(val)[1],
      if (is.data.frame(val)) paste(dim(val), collapse = " x ")
      else paste0("length ", length(val))
    ))
  }
  invisible(object)
})
