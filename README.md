# RBATS

The RBATS package is an R package that provides routines for sequential
inference within the class of univariate Normal Bayesian Dynamic Models.
Much of the code is based on the descriptions from
[West and Harrison (1997)](http://www2.stat.duke.edu/~mw/West&HarrisonBook/).

The main idea behind the package is that it implements Bayesian Dynamic Models,
maintaining Bayesian sequential inference while keeping a lightweight structure
with fewer dependencies.

In the current version, I translated the main functions,
[`forward_filter_dlm`](https://github.com/AndrMenezes/RBATS/blob/master/src/dlm.cpp#L138)
and [`backward_smoother_dlm`](https://github.com/AndrMenezes/RBATS/blob/master/src/dlm.cpp#L227), to `C++` using the
`Rcpp` and `RcppArmadillo` interfaces.

You can install the development version from GitHub with:
```
if (!require(remotes)) install.packages('remotes')
remotes::install_github("AndrMenezes/RBATS")
```

I mainly wrote RBATS for use in routine data analysis and for learning purposes.
I drew inspiration from the great [`PyBATS`](https://lavinei.github.io/pybats/)
python package by
[Isaac Lavine](https://www.linkedin.com/in/isaac-lavine-70495929/) and
[Andrew Cron](https://www.linkedin.com/in/andrewjcron/).


## Version Information

- v0.2.0: Translated the main code to `C++` and included new non-linear DLM models.

- v0.1.0: First release version written in `R` and used in [Migon et al. (2023)](https://onlinelibrary.wiley.com/doi/10.1002/asmb.2756).


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
