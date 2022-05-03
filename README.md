
<!-- README.md is generated from README.Rmd. Please edit that file -->

# leapfrog

<!-- badges: start -->

[![R-CMD-check](https://github.com/mrc-ide/leapfrog/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/leapfrog/actions)
[![Codecov test
coverage](https://codecov.io/gh/mrc-ide/leapfrog/branch/master/graph/badge.svg)](https://codecov.io/gh/mrc-ide/leapfrog?branch=master)
<!-- badges: end -->

Leapfrog is a multistate population projection model for demographic and
HIV epidemic estimation.

The name *leapfrog* is in honor of
[Professor](https://blogs.lshtm.ac.uk/alumni/2018/07/16/obituary-professor-basia-zaba/)
Basia
[Zaba](https://translate.google.co.uk/#view=home&op=translate&sl=pl&tl=en&text=%C5%BBaba).

*Note: the CCMPP model and implementation of Bayesian Population
Reconstruction in TMB previously here have moved to a separate
repository: <https://www.github.com/mrc-ide/ccmpp.tmb>*

## Simulation model

The simulation model is implemented in a header-only C++ library located
in [`inst/include/leapfrog-raw.h`](inst/include/leapfrog-raw.h). This
location allows the C++ code to be imported in other R packages via
specifying `LinkingTo: leapfrog` in the `DESCRIPTION` file.

The simulation model is callable in R via a wrapper function
`leafrogR()` created with [Rcpp](https://www.rcpp.org).

## Installation

Install the development version from
[GitHub](https://github.com/mrc-ide/leapfrog) via devtools:

``` r
# install.packages("devtools")
devtools::install_github("mrc-ide/leapfrog")
```

## Example

The file `pjnz/bwa2021_v6.13.pjnz` contains an example Spectrum file
constructed from default country data for Botswana in Spectrum v2.13
(December 2021).

Prepare model inputs.

``` r
library(leapfrog)

pjnz <- system.file("pjnz/bwa2021_v6.13.pjnz", package = "leapfrog")

demp <- prepare_leapfrog_demp(pjnz)
hivp <- prepare_leapfrog_projp(pjnz)
```

Simulate ‘full’ age group (single-year age) and ‘coarse’ age group
(collapsed age groups) models.

``` r
lsimF <- leapfrogR(demp, hivp, hiv_strat = "full")
lsimC <- leapfrogR(demp, hivp, hiv_strat = "coarse")
```

Compare the HIV prevalence age 15-49 years and AIDS deaths 50+ years.
Deaths 50+ years are to show some noticeable divergence between the
`"full"` and `"coarse"` age group simulations.

``` r
prevF <- colSums(lsimF$hivpop1[16:50,,],,2) / colSums(lsimF$totpop1[16:50,,],,2)
prevC <- colSums(lsimC$hivpop1[16:50,,],,2) / colSums(lsimC$totpop1[16:50,,],,2)

deathsF <- colSums(lsimF$hivdeaths[51:81,,],,2)
deathsC <- colSums(lsimC$hivdeaths[51:81,,],,2)

plot(1970:2030, prevF, type = "l", main = "Prevalence 15-49")
lines(1970:2030, prevC, col = 2)
```

<img src="man/figures/README-sim_prev-1.png" width="100%" />

``` r
plot(1970:2030, deathsF, type = "l", main = "AIDS Deaths 50+ years")
lines(1970:2030, deathsC, col = 2)
```

<img src="man/figures/README-sim_prev-2.png" width="100%" />

### Benchmarking

*Note: The function `devtools::load_all()` invokes
`pkgbuild::compile_dll()`, which executes with default arguemnt
`debug = TRUE`. This compiles the package with `-O0` compiler
optimisation flags. For benchmarking, make sure to install the package
or compile the DLL with `pkgbuild::compile_dll(debug = FALSE)`.*

``` r
library(bench)
library(eppasm)

fp <- eppasm::prepare_directincid(pjnz)

bench::mark(leapfrogR(demp, hivp, hiv_strat = "full"),
            leapfrogR(demp, hivp, hiv_strat = "coarse"),
            eppasm::simmod(fp),
            check = FALSE)
#> # A tibble: 3 × 6
#>   expression                                       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                                  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 leapfrogR(demp, hivp, hiv_strat = "full")     4.34ms   5.01ms      194.    4.33MB     45.2
#> 2 leapfrogR(demp, hivp, hiv_strat = "coarse")  761.4µs 967.22µs      971. 1006.65KB     34.8
#> 3 eppasm::simmod(fp)                            1.01ms   1.29ms      741.    1.43MB     32.8
```

## Code design

### Simulation model

The simulation model is implemented as templated C++ code in
`inst/include/leapfrog-raw.h`. This is the simulation model may be
developed as a standalone C++ library that can be called by other
software without requiring R-specific code features. The code uses
header-only open source libraries to maximize portability.

### R functions

The file `src/leapfrogR.cpp` contains R wrapper functions for the model
simulation via [Rcpp](http://dirk.eddelbuettel.com/code/rcpp.html) and
[RcppEigen](http://dirk.eddelbuettel.com/code/rcpp.eigen.html).

## Development notes

### Simulation model

-   The model was implemented using *Eigen::Tensor* containers. These
    were preferred for several reasons:

    -   Benchmarking found they were slighlty more efficient than
        *boost::multi_array*.
    -   Column-major indexing in the same order as R
    -   Other statistical packages (e.g. TMB, Stan) rely heavily on
        *Eigen* so using *Eigen* containers slims the dependencies.
