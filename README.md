
<!-- README.md is generated from README.Rmd. Please edit that file -->

# blantyreESBL

[![DOI](https://zenodo.org/badge/325831306.svg)](https://zenodo.org/badge/latestdoi/325831306)

<!-- badges: start -->
<!-- badges: end -->

The blantyreESBL R package is a repository for data generated as part of
a study into carriage of resistant bacteria in Blantyre, Malawi. It also
includes reproducible analysis script for three manuscripts arising from
the study:

[**Colonization dynamics of extended-spectrum beta-lactamase-producing
Enterobacterales in the gut of Malawian
adults**](https://joelewis101.github.io/blantyreESBL/articles/analysis.html)

Joseph M Lewis<sup>1,2,3,4</sup>, , Madalitso Mphasa<sup>1</sup>, Rachel
Banda<sup>1</sup>, Matthew Beale<sup>4</sup>, Eva Heinz<sup>2</sup>,
Jane Mallewa<sup>5</sup>, Christopher Jewell<sup>6</sup>, Nicholas R
Thomson<sup>4,7</sup>, Nicholas A Feasey<sup>1,2</sup>

Paper available at Nature Microbiology
[here](https://doi.org/10.1038/s41564-022-01216-7)

[**Genomic analysis of extended-spectrum beta-lactamase (ESBL) producing
Escherichia coli colonising adults in Blantyre, Malawi reveals
previously undescribed
diversity**](https://joelewis101.github.io/blantyreESBL/articles/analysis-ecoli.html)

Joseph M Lewis<sup>1,2,3,4</sup>, , Madalitso Mphasa<sup>1</sup>, Rachel
Banda<sup>1</sup>, Matthew Beale<sup>4</sup>, Jane Mallewa<sup>5</sup>,
Catherine Anscombe<sup>1,2</sup>, Allan Zuza<sup>1</sup>, Adam P
Roberts<sup>2</sup>, Eva Heinz*<sup>2</sup>, Nicholas
Thomson*<sup>4,7</sup>, Nicholas A Feasey\*<sup>1,2</sup>

-   = contributed equally

Preprint available at BioRxiv
[here](https://doi.org/10.1101/2021.10.07.463523)

[**Genomic and antigenic diversity of colonising *Klebsiella pneumoniae*
isolates mirrors that of invasive isolates in Blantyre,
Malawi**](https://joelewis101.github.io/blantyreESBL/articles/analysis-kleb.html)

Joseph M Lewis<sup>1,2,3,4</sup>, , Madalitso Mphasa<sup>1</sup>, Rachel
Banda<sup>1</sup>, Matthew Beale<sup>4</sup>, Jane Mallewa<sup>5</sup>,
Eva Heinz<sup>2</sup>, Nicholas Thomson<sup>4,7</sup>, Nicholas A
Feasey<sup>1,2</sup>

Available at Microbial Genomics
[here](https://doi.org/10.1099/mgen.0.000778)

### Installing and accessing data

If you just want the data, then all the data to replicate the analysis
are bundled with the package. To install the package from GitHub:

``` r
install.packages("devtools")
devtools::install_github("https://github.com/joelewis101/blantyreESBL")
```

The various data objects are described in the
[pkgdown](https://joelewis101.github.io/blantyreESBL/) site for this
package
[here](https://joelewis101.github.io/blantyreESBL/reference/index.html),
and available via R in the usual way (i.e. `?btESBL_participants` brings
up the definitions for the `btESBL_participants` data. They are all lazy
loaded so will be available for use immediately; they all start with
`btESBL_` to make it easy to choose the one you want using autocomplete.

### Whole genome sequence accession numbers and metadata

Reads from all isolates sequenced as part of this study have been
depositied in the European Nucleotide Archive (ENA); accession numbers
are available in the `btESBL_sequence_sample_metadata` data frame,
available on installing the package as above.

### Analysis scripts

The analysis scripts to reproduce tables and figures for each manuscript
are available as package vignettes; these can be built when downloading
the package by running:

``` r
devtools::install_github(
  "https://github.com/joelewis101/blantyreESBL", 
  build_vignettes = TRUE, 
  dependencies = TRUE )
```

The `dependencies = TRUE` option will install all the packages necessary
to run the vignette. Building the vignettes may take some time - you
have been warned!

Alternatively the source code for the vignettes are `analysis.Rmd`
`analysis-ecoli.Rmd` and `analysis-kleb.Rmd`in the `vignettes/` folder
of the [GitHub](https://github.com/joelewis101/blantyreESBL) repo or the
[pkgdown](https://joelewis101.github.io/blantyreESBL/) site for this
package has a rendered version of each vignette.

### Running Stan models

The longitudinal modelling paper uses models fit with Stan, the
probabilistic programming language, via the *rstan* R package. Unlike
the rest if the vignettes (which run the analysis as they are built) the
Stan models are not fit as part of the package vignettes as they take a
long time to fit. The outputs of the models are available as data
objects and a
[vignette](https://joelewis101.github.io/blantyreESBL/articles/fitting-stan-models.html)
provides instructions on how to fit and simulate from the posterior. The
Stan code is available as .stan files in the package directory - see
vignette for how to locate it.

### Author affiliations

1.  Malawi Liverpool Wellcome Research Programme, Blantyre, Malawi
2.  Liverpool School of Tropical Medicine, Liverpool, UK
3.  University of Liverpool, Liverpool, UK
4.  Wellcome Sanger Institute, Hinxton, UK
5.  Kamuzu University of Health Sciences, Malawi
6.  University of Lancaster, Lancaster, UK
7.  London School of Tropical Medicine and Hygiene, London, UK
