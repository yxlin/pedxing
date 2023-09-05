############################################################################## 80
## Multi-phase decision making
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
pkg <- c("data.table", "ggplot2", "rstan", "rethinking", "dplyr", "pedxing")
sapply(pkg, require, character.only = TRUE)
win_path <- "D:/03_Projects/pedxing"
unix_path <- "/media/yslin/kea/03_Projects/pedxing"
wk <- ifelse(.Platform$OS.type == "windows",
  shortPathName(win_path), unix_path
)


setwd(wk)
options(digits = 4)
rm(list = ls())

load("tests/extdata/mid-stage/proportion_exp.rda")
data <- data.frame(dp_be)

m0 <- map2stan(
  alist(
    count ~ dbinom(total, p),
    logit(p) <- a[s],
    a[s] ~ dnorm(0, 10)
  ),
  data = data, chain = 4, iter = 3000, warmup = 1500, cores = 4
)

m1 <- map2stan(
  alist(
    count ~ dbinom(total, p),
    logit(p) <- a[s] + b1 * TTAint,
    a[s] ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10)
  ),
  data = data, iter = 3000, warmup = 1500, chain = 4, cores = 4
)

m2 <- map2stan(
  alist(
    count ~ dbinom(total, p),
    logit(p) <- a[s] + b2 * side,
    a[s] ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10)
  ),
  data = data, iter = 3000, warmup = 1500, chain = 4, cores = 4
)

m3 <- map2stan(
  alist(
    count ~ dbinom(total, p),
    logit(p) <- a[s] + b1 * TTAint + b2 * side,
    a[s] ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10)
  ),
  data = data, iter = 3000, warmup = 1500, chain = 4, cores = 4
)

m4 <- map2stan(
  alist(
    count ~ dbinom(total, p),
    logit(p) <- a[s] + b1 * TTAint + b2 * side + b3 * TTAint * side,
    a[s] ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10)
  ),
  data = data, iter = 3000, warmup = 1500, chain = 4, cores = 4
)

save(m0, m1, m2, m3, m4, file = "tests/extdata/mid-stage/precar_given_safe_exp.RData")
