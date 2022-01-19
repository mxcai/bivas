bivas
===

bivas is an efficient statistical approach for bi-level variable selection.

Installation
===========

To install the development version of bivas, it's easiest to use the 'devtools' package. Note that bivas depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("mxcai/bivas")
```

Usage
===========
[The 'bivas' vignette](https://github.com/mxcai/bivas/blob/master/vignettes/bivas_package.pdf) provides a quick start for the usage of the package. The following help page also provides quick reference and examples:

```
library(bivas)
package?bivas
```

Reproducibility
==========

All the simulation results can be reproduced by using the code at [sim-bivas](https://github.com/mxcai/sim-bivas).

References
==========

Mingxuan Cai, Mingwei Dai, Jingsi Ming, Heng Peng, Jin Liu and Can Yang. BIVAS: A Scalable Bayesian Method for Bi-Level Variable Selection With Applications. Journal of Computational and Graphical Statistics. https://doi.org/10.1080/10618600.2019.1624365.


Development
==========

This R package is developed by Mingxuan Cai and Can Yang (macyang@ust.hk)
