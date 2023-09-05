##############################################################################80
## Multi-phrase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
pkg <- c('data.table', 'plotrix', 'car','afex', 'effsize', 'fBasics', 'infer',
         'emmeans', 'BayesFactor','ggbeeswarm', 'tibble')
sapply(pkg, require, character.only = TRUE)
rm (list = ls())
wk <- ifelse(.Platform$OS.type == "windows", shortPathName("E:/Documents/pedxing/"),
             '/media/yslin/Tui/Documents/pedxing/')

setwd(wk)
options(digits = 4)
load("data/pedxing.RData")
source("R/0_functions/utils_v1.R")
source("R/0_functions/clearup_data_v2.R")
load("R/2_preprocess/data/dpcar.RData"); 

## cb <- Manu::get_pal("Kaka")
## dmean_pcar <- dmean
## dd_pcar <- ddata

rvar <- 's'
wvar <- 'TTANme'
dv <- 'y'
dO <- dat[C==TRUE]

## Which cell resulted in 0 responses
tblsafe <- as_tibble(dO[, c('TTA', 'Side', 'G', 'R', 'RT', 's')])
tmp <- tblsafe %>% dplyr::count(TTA, G, s, .drop = FALSE)
dbtmp <- data.table(tmp[which(tmp$n==0), ])
dbtmp[, .N, .(TTA, G)]
#    TTA      G  N
# 1: 2.5 before 12  ## 12 participants did not respond at 2.5 s
# 2: 3.0 before  1
# 3: 4.0  after  1


## RM-ANOVA on safe proportions  ---------------------------------------------
dP0 <- dat[, .(y = mean(C)), .(TTANme, G, s)]   ## Safe proportion
aez0 <- aov_ez(rvar, dv, data = dP0[G == "before"], within = wvar)
aez1 <- aov_ez(rvar, dv, data = dP0[G == "after"],  within = wvar)
bf0 <- anovaBF(y ~ TTANme + s, data = dP0[G == "before"], whichRandom = rvar) 
bf1 <- anovaBF(y ~ TTANme + s, data = dP0[G == "after"], whichRandom = rvar) 
# [1] TTANme + s : 4.428e+23 ±0.76%
# [1] TTANme + s : 1.182 ±0.32%

dP030 <- dP0[G== "before" & TTANme == 3]
setorder(dP030, y)
dP030

dP025 <- dP0[G== "before" & TTANme == 2.5]
setorder(dP025, y)
dP025
## 2.5
##     TTANme      G  s       y
##  1:    2.5 before  4 0.00000
##  2:    2.5 before  7 0.00000
##  3:    2.5 before  9 0.00000
##  4:    2.5 before 11 0.00000
##  5:    2.5 before 12 0.00000
##  6:    2.5 before 15 0.00000
##  7:    2.5 before 17 0.00000
##  8:    2.5 before 18 0.00000
##  9:    2.5 before 19 0.00000
## 10:    2.5 before 20 0.00000
## 11:    2.5 before 21 0.00000
## 12:    2.5 before 22 0.00000
## 13:    2.5 before 13 0.02326
## 14:    2.5 before 14 0.05263
## 3.0
##     TTANme      G  s      y
##  1:      3 before 22 0.0000
##  2:      3 before 11 0.5000
##  3:      3 before 20 0.5385
##  4:      3 before 21 0.5833
##  5:      3 before 17 0.6250
##  6:      3 before  9 0.8696
##  7:      3 before 14 0.8750
##  8:      3 before 18 0.9200
##  9:      3 before 15 0.9216
## 10:      3 before 19 0.9412
## 11:      3 before 12 0.9643
## 12:      3 before  4 0.9667
## 13:      3 before 13 0.9773
## 14:      3 before  7 1.0000
