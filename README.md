
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nout

<!-- badges: start -->
<!-- badges: end -->

The goal of ***nout*** is to provide tools to quantify outliers through
closed testing procedures which use conformal *p*-values. According to
which local test is used, three closed testing procedures are available:
Simes version, Simes version with Storey estimator for the proportion of
true null hypotheses and Wilcoxon-Mann-Whitney version.

## Installation

You can install the development version of nout from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("chiaragaiamagnani/nout")
```

## Example

This is a basic example which shows the use of the three main procedures
contained in ***nout*** package:

``` r
library(nout)

set.seed(321)

# Generating score vectors
Sxy = sample(x=1:1000, size=100)
Sx = sample(Sxy, size=70)
Sy = setdiff(Sxy, Sx)

# Simes closed testing procedure
d_Simes(S_Y=Sy, S_X=Sx)
#> [1] 0
# Simes with Storey closed testing procedure
d_StoreySimes(S_Y=Sy, S_X=Sx)
#> [1] 0
# Wilcoxon-Mann-Whitney closed testing procedure
crit = critWMW(m=length(Sy), n=length(Sx))
d_mannwhitney(S_Y=Sy, S_X=Sx, crit=crit)
#> [1] 0
```

<!--You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
