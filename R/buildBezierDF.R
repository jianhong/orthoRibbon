#' Build the data.frame for B\'ezier curve
#' @param homolog_df_list A list with the data.frame of plot data for homologs.
#' @param colname1,colname2 The column names for top and bottom final positions.
#' @param col1,col2 The column names for color.
#' @param alpha The alpha value for color.
#' @param cl1 A numeric. The B\'ezier curve is created with two control points.
#' And cl1 should be no more than 4 and no less than 0.
#' If it is set to 4, all points in B\'ezier curve will be set to col1.
#' If it is set to 3, the top point, two control points will be set to col1 and
#' the bottom point will be set to col2.
#' If it is set to 2, the top point, top control point will be set to col1 and
#' the bottom control point, bottom point will be set to col2.
#' If it is set to 1, the top point will be set to col1 and others col2.
#' If it is set to 0, all points will be set to col2.
#' @returns A data.frame with the coordinates, colors and group id of
#' B\'ezier curve.
#' @export
#' @examples
#' extdata <- system.file('extdata', 'human_fly', package='orthoRibbon')
#' chrom_infos <- readRDS(file.path(extdata, 'chrom_infos.rds'))
#' homologs_df <- readRDS(file.path(extdata, 'homologs_df.rds'))
#' common_name <- setNames(nm=c("hsapiens", "dmelanogaster"))
#' plotData <- buildPlotData(common_name, homologs_df,
#'                           chrom_infos = chrom_infos,
#'                           chromosome_order_method = 'max',
#'                           max_links=1000)
#' bezier_chr <- buildBezierDF(plotData$homolog_df_list,
#'                             colname1='topChr_finalOffset',
#'                             colname2='bottomChr_finalOffset')
#'
buildBezierDF <- function(homolog_df_list, colname1, colname2,
                          col1='seq_top', col2='seq_bottom',
                          cl1=4, alpha=1){
  bezier_df <- lapply(seq_along(homolog_df_list), function(i){
    create_bezier_matrix(homolog_df_list[[i]],
                         i, colname1, colname2, col1, col2, cl1=cl1)
  })
  bezier_df <- do.call(rbind, bezier_df)
}
