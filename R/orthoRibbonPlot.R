#' Plot the data
#' @description
#' Plot data by ribbon plot.
#'
#' @param com_name Species abbreviations eg. "hsapiens", "mmusculus", "drerio".
#' @param bezier_df A data.frame. The coordinates for Bezier curves.
#' It must have columns "x", "y", "id", and "col".
#' The "id" column saves group id for each pair of homologs.
#' @param chrom_bars_df A data.frame. The coordinates for chromosome line.
#' It must have columns "x", "y", "sp" and "chrom". Here "sp" means species and
#' "chrom" indicates the chromosome.
#' @param chrom_label_df A data.frame. The coordinates for chromosome names.
#' It must have columns "x", "y" and "labels"
#' @param symbol_list_top,symbol_list_bottom The data.frame for the gene labels
#' @param link_lwd,chr_lwd The line width for Bezier curve, chromsome bar,
#' @param chr_lineend The line end for chromosome bar.
#' @param chr_size,label_size,symbol_size The size for chromosome label,
#' species, and symbols.
#' @param xlim Two numeric values, specifying the x limit of the scale.
#' @param show_symbol A logical value. Show the symbols or not.
#' @return A ggolot object
#' @importFrom ggplot2 ggplot aes geom_line geom_text scale_y_reverse xlim
#'  theme_minimal element_blank theme element_text coord_cartesian
#' @importFrom rlang .data
#' @importFrom ggforce geom_bezier2
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples
#' # example code
#'
orthoRibbonPlot <- function(com_name,
                            bezier_df, chrom_bars_df, chrom_label_df,
                            symbol_list_top, symbol_list_bottom,
                            link_lwd=0.25,
                            chr_lwd=3, chr_lineend='round',
                            chr_size = 6, label_size = 12, symbol_size = 2,
                            xlim=c(-0.02, 1.02),
                            show_symbol=FALSE){
  stopifnot(all(c('x', 'y', 'id', 'col') %in% colnames(bezier_df)))
  stopifnot(all(c('x', 'y', 'sp', 'chrom') %in% colnames(chrom_bars_df)))
  stopifnot(all(c('x', 'y', 'label') %in% colnames(chrom_label_df)))
  if(isTRUE(show_symbol)){
    stopifnot(all(c('x', 'y', 'top_label', 'bottom_label') %in%
                    colnames(symbol_list_top)))
    stopifnot(all(c('x', 'y', 'top_label', 'bottom_label') %in%
                    colnames(symbol_list_bottom)))
  }

  p <- ggplot() +
    ggforce::geom_bezier2(data = bezier_df,
                          aes(x = .data$x, y = .data$y,
                              group = .data$id, colour = .data$col),
                          linewidth = link_lwd, show.legend = FALSE)+
    geom_line(data = chrom_bars_df,
              aes(x = .data$x, y = .data$y,
                  group = interaction(.data$sp, .data$chrom),
                  colour = .data$chrom),
              linewidth = chr_lwd, lineend = chr_lineend,
              show.legend = FALSE) +
    geom_text(data = chrom_label_df,
              aes(x=.data$x, y=.data$y-0.03, label = .data$label),
              hjust = 0.5, vjust = 0, size = chr_size, show.legend = FALSE)
  if(isTRUE(show_symbol)){
    p <- p + scale_y_reverse(
      breaks = c(-.2, seq(0, length(com_name)-1), length(com_name)-.9),
      labels = c('symbol', com_name, 'symbol'),
      limits=c(-.3, length(com_name)-.8)) +
      geom_text_repel(data=symbol_list_top,
                      aes(x=.data$x, y=.data$y, label=.data$bottom_label),
                      force_pull = 0,
                      nudge_y=0.1,
                      direction = 'x',
                      angle=90,
                      hjust =0,
                      segment.size=0.2, max.iter=1e4,
                      max.time=1,
                      size=symbol_size)+
      geom_text_repel(data=symbol_list_top,
                      aes(x=.data$x, y=.data$y, label=.data$top_label),
                      force_pull = 0,
                      nudge_y=0.2,
                      direction = 'x',
                      angle=90,
                      hjust =0,
                      segment.size=0.2, max.iter=1e4,
                      max.time=1,
                      size=symbol_size) +
      geom_text_repel(data=symbol_list_bottom,
                      aes(x=.data$x, y=.data$y, label=.data$bottom_label),
                      force_pull = 0,
                      nudge_y=-0.1,
                      direction = 'x',
                      angle=90,
                      hjust =1,
                      segment.size=0.2, max.iter=1e4,
                      max.time=1,
                      size=symbol_size)
  }else{
    p <- p + scale_y_reverse(breaks = seq(0, length(com_name)-1),
                             labels = com_name)
  }
  p <- p + coord_cartesian(xlim=xlim) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = label_size))
  return(p)
}
