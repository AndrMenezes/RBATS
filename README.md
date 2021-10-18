
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RBATS

> The `RBATS` is a experimental R package for Bayesian Dynamic Models.
> The functionaly is based on the [Bayesian Forecasting and Dynamic
> Models](http://www2.stat.duke.edu/~mw/West&HarrisonBook/) from West
> and Harrison (1997). The currently version of the package implements
> the Dynamic Linear Models and Dynamic Generalized Exponential Growth
> Models.

## Development cycle

This package is early in development. Changes and additional
functionality may occur. The **development** version can be installed
from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("AndrMenezes/RBATS", build_vignettes = TRUE)
```

## Quick overview

The package is structured based on S3 classes objects. The following
functions created the two main objects of the package:

-   `dlm`: create and object of class `dlm` for Dynamic Linear Models.

-   `dgegm`: create and object of class `dgegm` for Dynamic Generalized
    Exponential Growth Models.

Each of this class objects have the following methods:

-   `update_moments`: update the moments of state parameters for `dlm`
    and `dgegm` objects.

-   `retrospective_moments`: smoothing moments of state parameters for
    `dlm` objects.

-   `fit`: perform the filtering and smoothing for `dlm` and `dgegm`
    objects.

-   `forecast`: perform the marginal forecast for `dlm.fit` objects.

-   `logLik`: return the predictive log-likelihood for `dlm.fit` and
    `dgegm.fit` objects.

Two vignettes showing the package functionality can be accessed through:

``` r
utils::browseVignettes(package = "RBATS")
```

## Contribute

The `RBATS` is still an experimental `R` package. I am inspired in the
great [`PyBATS`](https://lavinei.github.io/pybats/) python package from
[Isaac Lavine](https://www.linkedin.com/in/isaac-lavine-70495929/) and
[Andrew Cron](https://www.linkedin.com/in/andrewjcron/). I hope in
future `RBATS` became a robust tool as well as `PyBATS`.

Please, any suggestions and contributions are more than welcome. Feel
free to open issues, pull requests and/or forks. Thanks!

## Other R packages

Other R packages that can be used for Dynamic Linear Model or also known
State Space Models are:

-   [`dlm`](https://cran.r-project.org/web/packages/dlm/index.html) by
    Giovanni Petris and Wally Gilks.

-   [`bsts`](https://cran.r-project.org/web/packages/bsts/index.html) by
    Steven L. Scott.

-   [`KFAS`](https://cran.r-project.org/web/packages/KFAS/index.html) by
    Jouni Helske.

-   [`bssm`](https://cran.r-project.org/web/packages/bssm/index.html) by
    Jouni Helske and Matti Vihola.
