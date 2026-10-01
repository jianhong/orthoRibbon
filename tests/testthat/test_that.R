extdata <- system.file('extdata', package='orthoRibbon')
chrom_infos <- readRDS(file.path(extdata, 'ensembl', 'chrom_infos.rds'))
homologs_df <- readRDS(file.path(extdata, 'ensembl', 'homolog_df.human.mouse.zebrafish.rds'))
genes_gr <- readRDS(file.path(extdata, 'ensembl', 'gene.GRanges.obj.rds'))


# test_that("getGeneIDs", {
#   ids <- getGeneIDs('hsapiens')
#   expect_is(ids, 'list')
#   expect_equal(names(ids), 'hsapiens')
#   expect_true(all(grepl('^ENSG', ids[[1]])))
# })

test_that("filterChrom", {
  res <- filterChrom(chrom_infos[['mmusculus']], sp_min_chr_size = 1e8)
  expect_equal(nrow(res), sum(chrom_infos[['mmusculus']]$length>1e8))
})

test_that("getHomologIDs", {
  homologs <- list('a2b'=GRangesList(
   x=GRanges('seq1', IRanges(seq.int(5), width=1, name=letters[seq.int(5)]),
    homolog_ensembl_gene_ids=LETTERS[seq.int(5)]),
   y=GRanges('seq2', IRanges(seq.int(3), width=1, names=letters[seq.int(3)]),
    homolog_ensembl_gene_ids=c('m', 'n', 'k'))
  ))
  res <- getHomologIDs(homologs)
  expect_equal(colnames(res), c('gene_id1', 'gene_id2'))
  ## check if the bridge works
  expect_all_true(c('A m', 'B n', 'C k') %in% paste(res[, 1], res[, 2]))

  ## check single elements
  homologs <- list('a2b'=GRangesList(
    x=GRanges('seq1', IRanges(seq.int(5), width=1, name=letters[seq.int(5)]),
              homolog_ensembl_gene_ids=LETTERS[seq.int(5)])))
  res <- getHomologIDs(homologs)
  lapply(seq.int(5), function(i){
    expect_true(letters[i] %in% res[i, ])
    expect_true(LETTERS[i] %in% res[i, ])
  })
})

test_that('addGeneInfo', {
  homolog_df <- readRDS(file.path(extdata, 'human_fly/homologs_df.rds'))
  gene_gr <- readRDS(file.path(extdata, 'human_fly/genes_gr.rds'))
  res <- addGeneInfo(homolog_df = homolog_df,
                     genes_gr = gene_gr,
                     type = 'ortholog_group_with_gene_info')
  expect_equal(colnames(res), c(
    'gene_id1', 'gene_id2',
    'symbol1', 'symbol2',
    'species1', 'seq1', 'start1',
    'species2', 'seq2', 'start2'
  ))
  res <- addGeneInfo(homolog_df = homologs_df,
                     genes_gr = genes_gr,
                     type = 'ortholog_pair_only')
  expect_equal(colnames(res), c(
    'gene_id1', 'gene_id2',
    'symbol1', 'symbol2',
    'species1', 'seq1', 'start1',
    'species2', 'seq2', 'start2'
  ))
})

test_that('subsetHomologsByChrom', {
  homolog_df <- addGeneInfo(homolog_df = homologs_df,
                            genes_gr = genes_gr,
                            type = 'ortholog_pair_only')
  res <- subsetHomologsByChrom(homolog_df, chrom_infos, max_links=5000)
  expect_equal(nrow(res), 5000)

  res <- subsetHomologsByChrom(homolog_df,
                               chrom_infos[c('mmusculus', 'drerio')],
                               max_links=Inf)
  expect_true(all(res$species1!='hsapiens' & res$species2!='hsapiens'))
})

test_that('get_unique_max_rows', {
  m <- matrix(c(10, 1,
                1, 10),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("r1", "r2"), c("c1", "c2")))

  res <- get_unique_max_rows(m)
  expect_equal(res$column, c("c1", "c2"))
  expect_equal(res$row,    c("r1", "r2"))
  expect_equal(res$value,  c(10, 10))

  ## c1: r1=8, r2=1  -> colSum=9,  argmax=r1
  ## c2: r1=9, r2=8  -> colSum=17, argmax=r1
  ## Both want r1 -> conflict.
  ## Raw-count tiebreak would give r1 to c2 (9 > 8).
  ## Proportion tiebreak: c1 = 8/9 = 0.889, c2 = 9/17 = 0.529 -> c1 wins.
  m <- matrix(c(8, 9,
                1, 8),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("r1", "r2"), c("c1", "c2")))

  res <- get_unique_max_rows(m)
  res <- res[order(res$column), ]

  expect_equal(res$row, c("r1", "r2"))     # c1 -> r1, c2 -> r2
  expect_equal(res$value, c(8, 8))

  ## All three columns' argmax is r1 initially.
  ## c1: r1=5,r2=1,r3=1 -> colSum=7,  share=5/7=0.714
  ## c2: r1=6,r2=1,r3=1 -> colSum=8,  share=6/8=0.75   <- wins round 1
  ## c3: r1=4,r2=3,r3=1 -> colSum=8,  share=4/8=0.5
  ## Round 2 (rows r2,r3 remain; cols c1,c3 remain):
  ##   c1: r2=1,r3=1 -> tie, ties.method="first" -> r2
  ##   c3: r2=3,r3=1 -> argmax r2
  ##   Both want r2 -> conflict. shares: c1=1/7=0.143, c3=3/8=0.375 -> c3 wins.
  ## Round 3 (row r3 remains; col c1 remains): c1 -> r3, value=1.
  m <- matrix(c(5, 6, 4,
                1, 1, 3,
                1, 1, 1),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("r1", "r2", "r3"), c("c1", "c2", "c3")))

  res <- get_unique_max_rows(m)
  res <- res[order(res$column), ]

  expect_equal(res$row,   c("r3", "r1", "r2"))
  expect_equal(res$value, c(1, 6, 3))
  expect_equal(length(unique(res$row)), 3) # every row used exactly once

  m <- matrix(c(5, 3), nrow = 1,
              dimnames = list("r1", c("c1", "c2")))
  ## Only one row: both columns' argmax is r1 -> conflict.
  ## shares: c1 = 5/5 = 1.0, c2 = 3/3 = 1.0 -> tie -> first listed (c1) wins.
  ## c2 has no rows left afterwards -> NA.
  res <- get_unique_max_rows(m)

  expect_equal(res$column, c("c1", "c2"))
  expect_equal(res$row,    c("r1", NA_character_))
  expect_equal(res$value,  c(5, NA_real_))

  make_real_fixture <- function() {
    m <- rbind(
      `mmusculus 1`  = c(55,150,15,22,20,117,34,55,191,22,69,13, 43,12,38,17,43,14,5,92,7,88,34, 44,6),
      `mmusculus 2`  = c(25,27,30,28,165,138,97,92,128,47,72,6, 52,16,1,12,84,46,5,52,74,13,152, 26,65),
      `mmusculus 3`  = c(94,112,22,2,12,37,23,68,30,23,33,7, 26,58,41,114,17,29,83,13,7,48,36, 62,21),
      `mmusculus 4`  = c(43,86,52,4,59,73,71,84,1,28,73,7, 20,16,10,125,50,1,137,53,26,23,136, 17,4),
      `mmusculus 5`  = c(102,29,95,32,195,21,112,60,1,97,1,27, 12,66,38,18,34,16,15,58,60,5,19, 39,12),
      `mmusculus 6`  = c(31,25,47,167,49,66,23,51,11,69,54,13, 51,20,30,88,8,44,75,11,14,66,38, 18,74),
      `mmusculus 7`  = c(65,17,132,6,33,15,182,57,13,46,8,120, 28,7,127,115,25,140,55,2,48,6,16, 23,135),
      `mmusculus 8`  = c(152,50,50,13,35,5,140,38,24,41,39,4, 30,38,7,10,7,143,7,25,14,40,14, 6,82),
      `mmusculus 9`  = c(20,66,43,6,79,68,83,22,11,46,24,21, 31,3,102,85,10,124,34,18,36,28,21, 47,97),
      `mmusculus 10` = c(14,47,14,122,14,39,5,31,31,11,51,31, 67,5,3,22,56,33,8,124,8,75,69, 4,34),
      `mmusculus 11` = c(52,20,303,6,120,60,46,23,4,76,32,228, 39,75,134,23,16,3,30,9,129,12,15, 28,17),
      `mmusculus 12` = c(1,2,6,21,19,1,11,2,0,0,2,1, 79,0,10,15,237,2,15,194,0,1,6, 2,16),
      `mmusculus 13` = c(6,44,38,6,114,2,20,55,1,70,14,18, 13,39,0,44,5,4,56,14,71,45,10, 51,24),
      `mmusculus 14` = c(48,44,3,9,12,32,36,41,85,17,42,49, 64,8,11,11,51,5,16,37,5,10,18, 30,4),
      `mmusculus 15` = c(8,36,98,69,26,60,14,12,8,8,18,36, 4,2,1,76,2,17,80,9,26,23,72, 17,29),
      `mmusculus 16` = c(52,41,56,4,48,26,2,20,80,59,13,8, 4,1,44,1,2,10,3,7,14,37,2, 39,2),
      `mmusculus 17` = c(44,30,91,11,17,33,5,37,13,5,27,43, 95,9,15,29,54,6,60,61,22,48,20, 39,39),
      `mmusculus 18` = c(3,43,4,2,21,7,7,17,7,27,2,12, 0,62,4,22,0,1,27,13,77,5,3, 36,1),
      `mmusculus 19` = c(50,6,4,49,66,2,68,26,0,23,6,93, 103,55,6,5,38,5,5,5,51,15,1, 2,17),
      `mmusculus X`  = c(19,11,7,3,63,24,24,79,39,19,28,30, 4,135,1,8,2,3,5,4,52,6,76, 24,64),
      `mmusculus Y`  = c(1,0,0,0,3,2,0,53,2,1,0,0, 0,0,0,0,0,0,0,0,0,0,0, 3,76)
    )
    colnames(m) <- paste("drerio", 1:25)
    m
  }
  m <- make_real_fixture()
  res <- get_unique_max_rows(m)

  expect_equal(res$column, colnames(m))
  matched <- res[!is.na(res$row), ]
  expect_equal(length(unique(matched$row)), nrow(matched))
  expect_true(all(matched$row %in% rownames(m)))

  ## 21 rows, 25 columns -> exactly 4 columns must be left unmatched,
  ## regardless of which specific pairs the proportion rule picks
  ## (every round strictly consumes rows, so rows always deplete first)
  expect_equal(sum(is.na(res$row)), ncol(m) - nrow(m))

})

test_that('getChrOrders', {
  homolog_df <- addGeneInfo(homolog_df = homologs_df,
                               genes_gr = genes_gr,
                               type = 'ortholog_pair_only')
  homolog_df <- subsetHomologsByChrom(homolog_df, chrom_infos, max_links=Inf)
  ## only keep sex chromosome for human and mouse
  homolog_df <-
    homolog_df[(homolog_df$species1 %in% c('hsapiens', 'mmusculus') &
                  homolog_df$seq1 %in% c('X', 'Y')) |
                 (homolog_df$species2 %in% c('hsapiens', 'mmusculus') &
                    homolog_df$seq2 %in% c('X', 'Y')),  ]
  chrom_infos$hsapiens <-
    chrom_infos$hsapiens[chrom_infos$hsapiens$name %in% c('X', 'Y'), ]
  chrom_infos$mmusculus <-
    chrom_infos$mmusculus[chrom_infos$mmusculus$name %in% c('X', 'Y'), ]
  chr_orders <- lapply(c('max', 'spearman', 'TSP', 'GW', 'OLO', 'Spectral'),
                       function(chromosome_order_method){
    getChrOrders(homolog_df, chrom_infos,
                             method = chromosome_order_method)
    })
  names(chr_orders) <-
    c('max', 'spearman', 'TSP', 'GW', 'OLO', 'Spectral')
  chr_orders_hm <- lapply(chr_orders[-3], function(.ele){
    expect_all_true(.ele$hsapiens==.ele$mmusculus)
  })
  ## only keep X chromosome for human and mouse
  homolog_df <-
    homolog_df[(homolog_df$species1 %in% c('hsapiens', 'mmusculus') &
                  homolog_df$seq1 %in% c('X')) |
                 (homolog_df$species2 %in% c('hsapiens', 'mmusculus') &
                    homolog_df$seq2 %in% c('X')),  ]
  chrom_infos$hsapiens <-
    chrom_infos$hsapiens[chrom_infos$hsapiens$name %in% c('X'), ]
  chrom_infos$mmusculus <-
    chrom_infos$mmusculus[chrom_infos$mmusculus$name %in% c('X'), ]
  chr_orders <- lapply(c('max', 'spearman', 'TSP', 'GW', 'OLO', 'Spectral'),
                       function(chromosome_order_method){
                         getChrOrders(homolog_df, chrom_infos,
                                      method = chromosome_order_method)
                       })
})

test_that('buildChromDF & buildChromBarDF & buildChromLabelDF', {
  chrom_infos <- list(
    'hsapiens'=data.frame(
      'name'=c(1, 2, 3),
      'length'=rev(c(100, 200, 300))
    ),
    'mmusculus'=data.frame(
      'name'=c(1, 2, 3),
      'length'=rev(c(100, 200, 300))),
    'drerio'=data.frame(
      'name'=c(1, 2, 3, 4, 5),
      'length'=rev(c(100, 200, 300, 400, 500)))
  )
  chr_orders <- list(
    'hsapiens'=c(1, 2, 3),
    'mmusculus'=c(3, 2, 1),
    'drerio'=c(1, 2, 5, 4, 3)
  )
  chrDF <- buildChromDF(chrom_infos, chr_orders)
  mapply(function(ord, species){
    expect_equal(chrDF[[species]]$chrom, ord)
  }, chr_orders, names(chr_orders))
  lapply(chrDF, function(df){
    expect_equal(df$chrPlotPercent/min(df$chrPlotPercent),
                 df$chrsize/min(df$chrsize))
  })

  chrBar <- buildChromBarDF(chrDF)
  expect_equal(chrBar[match(names(chr_orders), chrBar$sp), 'y'], c(0, 1, 2))
  expect_all_true(c('x', 'y', 'chrom', 'label', 'sp') %in% colnames(chrBar))

  chrBarLabel <- buildChromLabelDF(chrBar)
  expect_all_true(c('x', 'y', 'chrom', 'label', 'sp') %in% colnames(chrBarLabel))
})

