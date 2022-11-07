# depthtestr: Two-sample testing with data depths

`depthtestr` is an R package for performing two-sample hypothesis tests using data depths. The main function is `depthTest`, which implements the Liu-Singh test (Liu and Singh, 1993). Different depth functions may be used in the hypothesis test; the `depthTest` function includes options for the following depths:

* Halfspace depth (`"halfspace"`) (Tukey, 1975)
* Mahalanobis depth (`"mahalanobis"`) (Zuo and Serfling, 2000)
* $L_p$ depth (`"lp"`) (Zuo and Serfling, 2000)
* Local depth from Paindaveine and Van Bever (`"pvb"`) (Paindaveine and Van Bever, 2013)
* Pairwise depth from Rendsburg and Garreau (`"pd"`) (Rendsburg and Garreau, 2021)
* Local community depth (`"lcd"`) (Evans and Berenhaut)
* Partitioned local depth (`"pald"`) (Berenhaut *et al.*, 2022)
* A user-defined custom depth function (`"custom"`)

## Installation

You can install `depthtestr` from GitHub using:

```r
# install.packages("devtools")
devtools::install_github("ciaran-evans/depthtestr")
```

## References

Berenhaut, K. S., Moore, K. E., and Melvin, R. L. (2022). A social perspective on perceived distances reveals deep community structure. *Proceedings of the National Academy of Sciences*, 119(4).

Evans, C. and Berenhaut, K.S. Two-sample testing with local community depth.

Liu, R. Y. and Singh, K. (1993). A quality index based on data depth and multivariate rank tests. *Journal of the American Statistical Association*, 88(421):252–260.

Paindaveine, D. and Van Bever, G. (2013). From depth to local depth: a focus on centrality. *Journal of the American Statistical Association*, 108(503):1105–1119.

Rendsburg, L. and Garreau, D. (2021). Comparison-based centrality measures. *International Journal of Data Science and Analytics*, 11(3):243–259.

Tukey, J. W. (1975). Mathematics and the picturing of data. In *Proceedings of the International Congress of Mathematicians, Vancouver, 1975*, volume 2, pages 523–531.

Zuo, Y. and Serfling, R. (2000). General notions of statistical depth function. *Annals of Statistics*, pages 461–482.
