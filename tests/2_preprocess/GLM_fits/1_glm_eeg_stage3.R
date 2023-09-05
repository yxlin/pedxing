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

# pre-car RT | Safe ---------------------------
load("tests/extdata/mid-stage/proportion_rt_eeg.rda")
data <- data.frame(dbe)

m0 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s],
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m1 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m2 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m3 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m4 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side + b3 * TTAint * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)


# pre-car RT | Safe
save(m0, m1, m2, m3, m4, file = "tests/extdata/mid-stage/safe_rt_given_precar_eeg.RData")

# postcar RT | Safe ---------------------------
data <- data.frame(daf)
m0 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s],
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m1 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m2 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m3 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m4 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side + b3 * TTAint * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

# post-car RT | Safe
save(m0, m1, m2, m3, m4, file = "tests/extdata/mid-stage/safe_rt_given_postcar_eeg.RData")

# pre-car RT | Unsafe ---------------------------
data <- data.frame(dbe_hit)
m0 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s],
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m1 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m2 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m3 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m4 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side + b3 * TTAint * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

save(m0, m1, m2, m3, m4, file = "tests/extdata/mid-stage/unsafe_rt_given_precar_eeg.RData")

# post-car RT | Unsafe ---------------------------
data <- data.frame(daf_hit)

m0 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s],
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m1 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m2 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m3 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)

m4 <- map2stan(
  alist(
    value ~ dnorm(mu, sigma),
    mu <- a[s] + b1 * TTAint + b2 * side + b3 * TTAint * side,
    a[s] ~ dnorm(amu, sigma_a),
    amu ~ dnorm(0, 10),
    b1 ~ dnorm(0, 10),
    b2 ~ dnorm(0, 10),
    b3 ~ dnorm(0, 10),
    c(sigma, sigma_a) ~ dcauchy(0, 2)
  ),
  data = data, iter = 5000, warmup = 2500, chains = 4, cores = 4
)



# post-car RT | Unsafe
save(m0, m1, m2, m3, m4, file = "tests/extdata/mid-stage/unsafe_rt_given_postcar_eeg.RData")
