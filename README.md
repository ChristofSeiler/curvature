# R package `curvature`

This package implements a simulatio-based approach to our NIPS paper:

```
Positive Curvature and Hamiltonian Monte Carlo 
C. Seiler, S. Rubinstein-Salzedo, and S. Holmes 
NIPS, December, 2014, Montreal, Canada
```

and our long version:

```
Curvature and Concentration of Hamiltonian Monte Carlo in High Dimensions 
S. Holmes, S. Rubinstein-Salzedo, and C. Seiler 
https://arxiv.org/abs/1407.1114
```

## Installation

``` r
install.packages("devtools")
devtools::install_github("ChristofSeiler/curvature")
```

## Getting Started

``` r
library("curvature")
fit = rstan::sampling(...)
res = curvature(fit)
summary(res)
plot(res)
```
