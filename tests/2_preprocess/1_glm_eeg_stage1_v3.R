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
data(ped)


# Load data ---------------------------------------------------
data(ped) # or load("data/pedxing.RData")


## Logistic models ----------------
### Pre-car, and safe proportion ------------------
dbeh <- d[E == "eeg"]
unique(dbeh$s)
sort(unique(dbeh$sold))
table(dbeh$sold)
# 2  4  7  9  11 12 13 14 15 17 18 19 20 21 22
# 2  4  7  9  11 12 13 14 15 17 18 19 20 21 22 1  3  5  6  8  10 16 23 24
dbeh$s <- ifelse(dbeh$sold == 1, 1,
  ifelse(dbeh$sold == 2, 2,
    ifelse(dbeh$sold == 3, 3,
      ifelse(dbeh$sold == 4, 4,
        ifelse(dbeh$sold == 5, 5,
          ifelse(dbeh$sold == 6, 6,
            ifelse(dbeh$sold == 7, 7,
              ifelse(dbeh$sold == 8, 8,
                ifelse(dbeh$sold == 9, 9,
                  ifelse(dbeh$sold == 10, 10,
                    ifelse(dbeh$sold == 11, 11,
                      ifelse(dbeh$sold == 12, 12,
                        ifelse(dbeh$sold == 13, 13,
                          ifelse(dbeh$sold == 14, 14,
                            ifelse(dbeh$sold == 15, 15,
                              ifelse(dbeh$sold == 16, 16,
                                ifelse(dbeh$sold == 17, 17,
                                  ifelse(dbeh$sold == 18, 18,
                                    ifelse(dbeh$sold == 19, 19,
                                      ifelse(dbeh$sold == 20, 20,
                                        ifelse(dbeh$sold == 21, 21,
                                          ifelse(dbeh$sold == 22, 22,
                                            ifelse(dbeh$sold == 23, 23,
                                              ifelse(dbeh$sold == 24, 24, NA)
                                            )
                                          )
                                        )
                                      )
                                    )
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)

unique(dbeh$sold)
sort(unique(dbeh$s))

length(unique(dbeh$s))
length(unique(dbeh$sold))

ns <- length(unique(dbeh$s))
ns


pro <- dbeh[R == "safe", .(count = .N), .(TTAint, TTANme, side, G, s)]
pro[TTAint == 0 & G == "before"]
pro[, total := sum(count), .(TTAint, side, s)]
pro[, value := count / total]

# pre-car percentage | Safe responses
# post-car percentage | Safe responses
dp_be <- pro[G == "before"]
dp_af <- pro[G == "after"]

pro <- dbeh[G == "before", .(count = .N), .(TTAint, TTANme, side, R, s)]
pro[, total := sum(count), .(TTAint, side, s)]
pro[, value := count / total]

# Safe prop | pre-car
dsafe_be <- pro[R == "safe"]

## replace with data.frame(dsafe_af) to get safe proportion | after the car
pro <- dbeh[G == "after", .(count = .N), .(TTAint, TTANme, side, R, s)]
pro[, total := sum(count), .(TTAint, side, s)]
pro[, value := count / total]

# Safe prop | post-car
dsafe_af <- pro[R == "safe"]


tibble::as_tibble(dp_be)



### RT before and after ----------------
## RT before and after | Safe
rt <- dbeh[R == "safe", .(
  count = .N,
  value = median(RT)
), .(TTAint, TTANme, side, G, s)]

# pre-car and post-car  RT | Safe
dbe <- rt[G == "before"]
daf <- rt[G == "after"]

## RT before and after | Unsafe
rt <- dbeh[R == "hit", .(
  count = .N,
  value = median(RT)
), .(TTAint, TTANme, side, G, s)]

# pre-car and post-car RT | Unsafe
dbe_hit <- rt[G == "before"]
daf_hit <- rt[G == "after"]

save(dp_be, dp_af, dsafe_be, dsafe_af, dbe, daf, dbe_hit, daf_hit,
  file = "tests/extdata/mid-stage/proportion_rt_eeg.rda"
)
