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
cols <- c("TTA", "Side", "R", "G", "TTANme", "RT", "s", "sold", "E")
d[, ..cols]
table(d$s)
d[, .N, .(E, TTA, s)]
unique(d[, .N, .(E, TTA, s)]$N)

# Basic checks ------------
dhit <- d[R == "hit"] %>% count(TTA, G, .drop = FALSE)
dhit

prohit <- d[R == "hit", .(
  count = .N,
  MRT = mean(RT)
), .(TTA, Side, G, s)]
prohit[, total := sum(count), .(TTA, Side, s)]
prohit[, value := count / total]
setorder(prohit, TTA, G, Side, s)
tibble::as_tibble(prohit[TTA == 2.5 & G == "before"])

# 44.2/ c(16.00, 13.33, 11.43, 10.00)
## Road width: 4.2m
## Initial distance from Kerb: 0.5 m
## Car Length: 4.96 m
## Car Width: 1.9 m
## Car initial distance: 40 m


## Logistic models ----------------
### Pre-car, and safe proportion ------------------
dbeh <- d[E == "exp"]
unique(dbeh$s)
sort(unique(dbeh$sold))
# 2  4  7  9  11 12 13 14 15 17 18 19 20 21 22
dbeh$s <- ifelse(dbeh$sold == 2, 1,
  ifelse(dbeh$sold == 4, 2,
    ifelse(dbeh$sold == 7, 3,
      ifelse(dbeh$sold == 9, 4,
        ifelse(dbeh$sold == 11, 5,
          ifelse(dbeh$sold == 12, 6,
            ifelse(dbeh$sold == 13, 7,
              ifelse(dbeh$sold == 14, 8,
                ifelse(dbeh$sold == 15, 9,
                  ifelse(dbeh$sold == 17, 10,
                    ifelse(dbeh$sold == 18, 11,
                      ifelse(dbeh$sold == 19, 12,
                        ifelse(dbeh$sold == 20, 13,
                          ifelse(dbeh$sold == 21, 14,
                            ifelse(dbeh$sold == 22, 15, NA)
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


dp_be[TTAint == 0]
# TTAint TTANme side      G  s count total   value
# 1:      0    2.5    0 before  7     1    10 0.10000
# 2:      0    2.5    1 before  9     1    22 0.04545
# 3:      0    2.5    0 before 16     1    18 0.05556
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

dbe[s == 1]
dbe[s == 2]
dbe[s == 3]

daf[s == 1]

# save(dp_be, dp_af, dsafe_be, dsafe_af, file = "tests/extdata/mid-stage/proportion_exp.rda")
# save(dbe, daf, dbe_hit, daf_hit, file = "tests/extdata/mid-stage/rt_exp.rda")
