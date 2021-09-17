
<!-- README.md is generated from README.Rmd. Please edit that file -->

# blantyreESBL

<!-- badges: start -->
<!-- badges: end -->

The blantyreESBL R package contains data and code to replicate the
analysis of the manuscript:

<br />

*Dynamics of gut mucosal colonisation with extended spectrum
beta-lactamase producing Enterobacterales in Malawi*

<br />

Joseph M Lewis<sup>1,2,3,4</sup>, , Madalitso Mphasa<sup>1</sup>, Rachel
Banda<sup>1</sup>, Matthew Beale<sup>4</sup>, Eva Heinz<sup>2</sup>,
Jane Mallewa<sup>5</sup>, Christopher Jewell<sup>6</sup>, Nicholas R
Thomson<sup>4</sup>, Nicholas A Feasey<sup>1,2</sup>

1.  Malawi Liverpool Wellcome Clinical Research Programme, Blantyre,
    Malawi
2.  Department of Clinical Sciences, Liverpool School of Tropical
    Medicine, Liverpool, UK
3.  Department of Clinical Infection, Microbiology and Immunology,
    University of Liverpool, Liverpool, UK
4.  Wellcome Sanger Institute, Hinxton, UK
5.  College of Medicine, University of Malawi, Malawi
6.  University of Lancaster, Lanaster, UK

### Installing and accessing data

If you just want the data, then all the data to replicate the analysis
are bundled with the package. To install the package from GitHub:

``` r
install.packages("devtools")
devtools::install_github("https://github.com/joelewis101/blantyreESBL)
```

The various data objects are described in the
[pkgdown](https://joelewis101.github.io/blantyreESBL/) site for this
package, and available via R in the usual way
(i.e. `?btESBL_participants` brings up the definitions for the
`btESBL_participants` data. They are all lazy loaded so will be
available for use immediately; they all start with `btESBL_` to make it
easy to choose the one you want using autocomplete.

The analysis is available as a package vignette; this can be built when
downloading the package by typing:

``` r
devtools::install_github("https://github.com/joelewis101/blantyreESBL", build_vignettes = TRUE, dependencies = TRUE )
```

The `dependencies = TRUE` option will install all the packages necessary
to run the vignette.

Alternatively the source code for the vignette is `analysis.Rmd` in the
`vignettes/` folder of the
[GitHub](https://github.com/joelewis101/blantyreESBL) repo or the
[pkgdown](https://joelewis101.github.io/blantyreESBL/) site for this
package has a rendered version.
