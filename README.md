
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Overview

`lorenzgini` is an [R](https://www.r-project.org) package that
calculates the Gini coefficient and plots Lorenz curves. It builds on
the `gini` function from package `DescTools`, and extends this function
to the generalized Gini as derived in:

> Bernasco, W. and W. Steenbeek (2017). More places than crimes:
> Implications for evaluating the law of crime concentration at place.
> *Journal of Quantitative Criminology*.
> <https://doi.org/10.1007/s10940-016-9324-7>

## Version

The most recent version of the package is:

> Steenbeek, W., Bernasco, W. (2018). *lorenzgini: generalized Gini for
> sparse data situations*. R package version 0.1.2. URL:
> <https://github.com/wsteenbeek/lorenzgini>

## Installation

You can install the package from this [GitHub
repository](https://github.com/wsteenbeek/lorenzgini). You first need to
install the [remotes](https://CRAN.R-project.org/package=remotes)
package.

``` r
install.packages("remotes")
```

Then install `lorenzgini` using the `install_github` function in the
[remotes](https://CRAN.R-project.org/package=remotes) package.

``` r
remotes::install_github("wsteenbeek/lorenzgini")
```

## How to start

After loading the package with `library(lorenzgini)`, see the
documentation of the two functions (these include examples) by typing:

  - `?gini`
  - `?lorenz`

## Vignette

The package includes a vignette explaining the background and the
functions in `lorenzgini` in more detail.

By far the easiest way to view the vignette is this direct link,
courtesy of the [GitHub HTML Preview
service](http://htmlpreview.github.io/):

  - [Introduction to
    lorenzgini](http://htmlpreview.github.io/?https://github.com/wsteenbeek/lorenzgini/blob/main/doc/lorenzgini.html)

If instead you want to access the vignettes from R itself you need to
take a few additional steps, because `remotes::install_github()` does
not build vignettes by default to save time and because it may require
additional packages.

1.  Install the rmarkdown package with `install.packages("rmarkdown")`

2.  [Install pandoc](http://johnmacfarlane.net/pandoc/installing.html)
    (and afterwards restart your computer)

3.  Then, install the package again but force building of the vignettes
    using `remotes::install_github("wsteenbeek/lorenzgini",
    build_vignettes = TRUE, force = TRUE)`. This will take a few
    minutes.

Afterwards, you should be able to view which vignettes are available
using:

``` r
browseVignettes("lorenzgini")
```

To directly read the vignettes rather than going through
`browseVignettes("lorenzgini")` you can use:

``` r
vignette("lorenzgini", package = "lorenzgini")
```

## License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3, as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the GNU General
Public License for more details.

A copy of the GNU General Public License, version 3, is available at
<https://www.r-project.org/Licenses/GPL-3>
