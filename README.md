# orthoRibbon

## Ribbon-based visualization of orthologous genes across species

OrthoRibbon is an R package designed to visualize orthologous gene relationships across species using ribbon diagrams. The package connects homologous genes between genomes through ribbons linking their genomic coordinates, enablingcomparison of conserved genes and genomic organization across species.

OrthoRibbon transforms ortholog tables into clear cross-species visualizations and produces customizable, publication-quality figures suitable for comparative genomics and evolutionary analysis.

## Features

Visualize orthologous gene relationships across species

Ribbon-based representation of cross-species genomic correspondence

Flexible support for multiple species and chromosomes

Customizable ribbon aesthetics and genomic layouts

Generate publication-quality figures

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("jianhong/orthoRibbon",
                        build_vignettes = TRUE,
                        dependencies = TRUE)
```

## Quick Start

Example workflow to visualize orthologous genes across species.

```r
library(orthoRibbon)

# example ortholog table
ortho <- data.frame(
  species1 = c("Human", "Human", "Human"),
  gene1    = c("TP53", "BRCA1", "MYC"),
  species2 = c("Mouse", "Mouse", "Mouse"),
  gene2    = c("Trp53", "Brca1", "Myc")
)

# plot ribbon diagram
orthoRibbon(ortho)
```

## Input Data

OrthoRibbon requires genomic coordinate information to position chromosomes and genes across species and to draw ribbons connecting orthologous genes.

### Method 1: Supply orthologous gene pairs along with their corresponding genomic information.

### 1. An ortholog table describing gene relationships across species.

Typical fields include:

| Column         | Description               |
| :------------- | :------------------------ |
| gene_id1       | gene identifier           |
| gene_id2       | gene identifier           |

Users can obtain ortholog information from databases or orthogroups searching tools such as:

* [Ensembl](http://www.ensembl.org/)

* [OrthoDB](https://www.orthodb.org/)

* [eggNOG](http://eggnog6.embl.de/)

* [OrthoFinder](https://github.com/OrthoFinder/OrthoFinder)

* [odp](https://github.com/conchoecia/odp)

### 2. Chromosome Length Information

Chromosome length information defines the genomic layout for each species. These data are used to determine the relative sizes and positions of chromosomes in the visualization.

The chromosome length information should be provided as a list, where each element corresponds to a species and is named by the species identifier.

Each element of the list should be a chromosome table containing the following columns:

| Column       | Description                       |
| :----------- | :-------------------------------- |
| `name`       | Chromosome or scaffold name       |
| `length`     | Chromosome length (in base pairs) |

### 3. Gene Information

Gene information specifies the genomic coordinates of genes used in the ortholog visualization. These coordinates determine where ribbons will attach on each chromosome.

Gene information should be provided as a `GRanges` object, with `gene_name` and `species` stored as metadata columns, and the gene identifier assigned as the name of each range.

### Method 2: Supply orthologous gene groups along with the corresponding genomic information.

OrthoRibbon expects an ortholog table describing gene relationships across species.

Typical fields include:

| Column         | Description               |
| :------------- | :------------------------ |
| species        | species name              |
| gene_id        | gene identifier           |
| seq            | chromosome name           |
| start          | genomic start coordinate  |
| ortholog_group | ortholog group identifier |

And the chromosome length information as above.

## Visualization

OrthoRibbon generates ribbon diagrams where:

* chromosomes or genomic regions are displayed as tracks

* orthologous genes are connected by ribbons

* ribbons represent cross-species gene relationships

These visualizations help reveal:

* conserved genomic regions

* ortholog clusters

* evolutionary relationships across species.

## Documentation

Full tutorial and examples are available in the vignette:

👉 https://jianhong.github.io/orthoRibbon/articles/orthoRibbon.html

To view the documentation locally:

```r
browseVignettes("orthoRibbon")
```

## Contributing

Contributions are welcome.

Typical workflow:

1. [Fork the repository](https://github.com/jianhong/orthoRibbon/fork)

2. Create a feature branch

3. [Submit a pull request](https://github.com/jianhong/orthoRibbon/compare)

Bug reports and feature requests can be submitted via GitHub [Issues](https://github.com/jianhong/orthoRibbon/issues/new/choose).



