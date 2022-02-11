# iPheGWAS: bringing intelligence to visualizing the genome-phenome wide association results.

<!-- badges: start -->
[![R-CMD-check](https://github.com/georgeg0/PheGWAS2/workflows/R-CMD-check/badge.svg)](https://github.com/georgeg0/PheGWAS2/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

First of all, thanks to the community for using [PheGWAS](https://www.biorxiv.org/content/10.1101/694794v2.article-info). Since it's release we have been working on a heuristic approach to give quick insights into the genetic relatedness for phenotypes using GWAS summary statistics file. We are adding this module to PheGWAS and reintroducing it as iPheGWAS (Intelligent PheGWAS). The traits in your landscape are no more ordered as random but will be ordered based on their genetic correlations. 

In addition to the heuristic approach that we developed, all the functionalities outlined in the ***PheGWAS*** are also available in ***iphegwas*** package. Considering our community's request to improve the speed of PheGWAS, the entire codebase is rewritten, and you will notice that the ***iphegwas*** package is significantly faster than the ***PheGWAS*** package. 

## Installation:
### Prerequisites
Please Install plotly, devtools and Biomart before installing the PheGWAS package.

Install plotly
```
install.packages("plotly")
```

Install devtools
```
install.packages("devtools")
```

Install BioMart
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```
## Installation of iphegwas package
```
install_github("georgeg0/iphegwas")
```
