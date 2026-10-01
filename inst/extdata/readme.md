# Pre-prepared data from Ensembl

OrthoRibbon supports flexible input formats including commonly used databases such as Ensembl.
OrthoRibbon can be readily adapted to using the homolog pairs like data retrieved from Ensembl database.

Here is the code to retrieve chromosome, homologs and gene coordinates via BioMart at Mar 13, 2026.

```r
# prepare mart for BioMart
marts <- geneClusterPattern::guessSpecies(
  common_name, output='mart')

# prepare all the Ensembl ids
ids <- orthoRibbon::getGeneIDs(common_name, marts)

# get chromosome info
chrom_infos<- lapply(common_name, GenomeInfoDb::getChromInfoFromEnsembl)

# retrieve all the homologs
homologs <- orthoRibbon::getHomologGRs(ids, common_name, marts)

# retrieve all the gene ranges
genes_gr <- orthoRibbon::getGeneGRs(ids, marts, homologs)

# reformat the homologs to data.frame
homologs_df <- getHomologIDs(homologs)
```

And modified at Oct 1st, 2026 to reduce the size.

```r
set.seed(42)
mouse_ids <- unique(as.character(unlist(homologs_df)))
mouse_ids <- mouse_ids[grepl('ENSMUSG', mouse_ids)]
idx <- sample(mouse_ids, 5000)
keep <- homologs_df[, 1] %in% idx | homologs_df[, 2] %in% idx
homologs_df <- homologs_df[keep, ]
keep <- names(genes_gr) %in% homologs_df[, 1] |
  names(genes_gr) %in% homologs_df[, 2]
genes_gr <- genes_gr[keep]
```

# human fly orthoFinder results

The ortholog groups were created via orthoFinder for human and fly.

