##############################################################################80
## Multi-phase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
pkg <- c('data.table', 'plotrix', 'car','afex', 'emmeans', 'BayesFactor','pedxing')
sapply(pkg, require, character.only = TRUE)
rm (list = ls())
wk <- ifelse(.Platform$OS.type =="windows", shortPathName("E:/Documents/pedxing/"),
    '/media/yslin/Tui/Documents/pedxing/')
# install.packages("effsize")
# install.packages("gamlss.dist")
setwd(wk)
options(digits = 4)
dat[E=="exp" & G == "after" & R == "hit"]
dat[E=="eeg" & G == "after" & R == "hit"]

dat[G == "after" & R == "hit" & TTA != 2.5]$RT - 
  (dat[G == "after" & R == "hit" & TTA != 2.5]$TTA + dat[G == "after" & R == "hit" & TTA != 2.5]$Jitter)

dat[G == "after" & R == "hit" & TTA == 2.5]$RT - 
  (dat[G == "after" & R == "hit" & TTA == 2.5]$TTA + dat[G == "after" & R == "hit" & TTA == 2.5]$Jitter)

dat[G == "after" & R == "hit" & TTA == 2.5]

(4.2-1.9)/1.6
1.9/1.6
100*(7/nrow(dat))

# require(gamlss.dist)
# x <- seq(0, 5, .1) + .1
# y <- dexGAUS(x, mu = 1, nu = .5)
# y <- dexGAUS(x, mu = 1, nu = 1)
# plot(x, y)
# psych::describe(y)
dunsafe <- dat[R == "hit" & E == "exp" & TTA == 2.5 & G=="after"]

dunsafe <- dat[R == "safe" & E == "exp" & TTA == 2.5 & G=="after"]

dunsafe <- dat[R == "hit" & E == "exp"]
idx <- which(dat$R == "hit" & dat$E == "exp" & dat$G == "after")
dat[idx]
(2.5-0.02849 + 4.96/16) + (4.2/1.6)

3.973 - (2.5-0.02849 + 4.96/16)
dat[354,]

(4.2-1.9)/1.6

dunsafe[G=="after", .N, .(TTA, G)]

dunsafe <- dat[R == "hit" & E == "eeg"]
dunsafe[G=="after", .N, .(TTA, G)]

## E1 --------------------------------------------------------
# d0[s %in% c(9, 13, 17, 18)]
# load("data/pedxing.RData")
# source("R/0_functions/utils_v1.R")
## 3561 rows
data(ped)
dsafe <- dat[R == "safe" & E == "exp"]
psych::describeBy(dsafe[G=="before"]$RT, dsafe[G=="before"]$TTA)
dat[TTA == 2.5 & R == "safe" & G == "before"]
dat[TTA == 3 & R == "safe" & G == "before"]
dtmp <- dat[, .N, .(TTA, R, G)]
setorder(dtmp, TTA, R, G)


dtmp <- dat[R == "safe" & E == "exp", .(y = psych::skew(RT),
                N = length(RT)), .(TTA, R, G, s)]
dd <- dtmp[, .(yy = mean(y)), .(TTA, R, G)]
setorder(dd, TTA)
dd[R == "safe" & G == "before"]
#    TTA    R      G     yy
# 1: 2.5 safe before     NA
# 2: 3.0 safe before 0.2908
# 3: 3.5 safe before 0.8763
# 4: 4.0 safe before 2.0691

# TTA  emmean    SE df lower.CL upper.CL
# X3.5  0.774 0.267 21    0.219     1.33
# X4    1.851 0.372 21    1.077     2.62
# X3    0.590 0.234 21    0.104     1.08
idx <- which( is.na(dtmp$y) )
dtmp2 <- dtmp[!idx, ]
dtmp3 <- dtmp2[TTA != 2.5 & G == "before" & R == "safe"]
dtmp3$TTA <- factor(dtmp3$TTA, levels = c(3, 3.5, 4))

rvar <- "s"
wvar <- "TTA"
aez0 <- aov_ez(rvar, "y", data = dtmp3,  within = wvar, 
               anova_table = list(correction = "none", es = "none"))
aez0 <- aov_ez(rvar, "y", data = dtmp3,  within = wvar, 
               anova_table = list(es = "pes"))
aez0 <- aov_ez(rvar, "y", data = dtmp3,  within = wvar, 
               anova_table = list(es = "ges"))

# Anova Table (Type 3 tests)
# 
# Response: y
# Effect          df   MSE      F       p.value
#    TTA    2,    26  1.26  10.36***       <.001
#        1.61, 20.98  1.56  10.36**         .001
# pes = .444
# ges = 320

dtmp3$s <- factor(dtmp3$s)
bf <- BayesFactor::anovaBF(y ~ TTA + s, data = dtmp3, whichRandom = rvar)
# Bayes factor analysis
#   [1] TTA + s : 162.2 ±0.53%

r1 <- emmeans(aez0, ~TTA)
contrast(r1, interaction = "pairwise", adjust = "bonferroni")
# TTA_pairwise estimate    SE df t.ratio p.value
# X3 - X3.5      -0.701 0.307 13  -2.283  0.1197
# X3 - X4        -1.907 0.452 13  -4.222  0.0030
# X3.5 - X4      -1.206 0.490 13  -2.459  0.0862

col1 <- dtmp3[s != 22 & TTA %in% c(3)]$y
col2 <- dtmp3[s != 22 & TTA %in% c(3.5)]$y
bftmp1 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)

col1 <- dtmp3[s != 22 & TTA %in% c(3)]$y
col2 <- dtmp3[s != 22 & TTA %in% c(4)]$y
bftmp2 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)

col1 <- dtmp3[s != 22 & TTA %in% c(3.5)]$y
col2 <- dtmp3[s != 22 & TTA %in% c(4)]$y
bftmp3 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)
# [1] Alt., r=0.707 : 1.875 ±0% # 3 v 3.5
# [1] Alt., r=0.707 : 39.21 ±0% # 3 v 4
# [1] Alt., r=0.707 : 2.426 ±0% # 3.5 v 4
extractBF(bftmp1)
extractBF(bftmp2)
extractBF(bftmp3)


# TTA  emmean    SE df lower.CL upper.CL
# X3    0.291 0.108 13   0.0576    0.524
# X3.5  0.992 0.311 13   0.3205    1.663
# X4    2.197 0.442 13   1.2420    3.153


# dSkewness <- afex_plot(aez0, x = wvar, error = 'within', return = "data")
dSkewness <- afex_plot(aez0, x = wvar, error = 'within')
plot(dSkewness)


## E2----------------
dtmp <- dat[R == "safe" & E == "eeg", .(y = psych::skew(RT),
                                        N = length(RT)), .(TTA, R, G, s)]
dd <- dtmp[, .(yy = mean(y)), .(TTA, R, G)]
setorder(dd, TTA)
dd[R == "safe" & G == "before"]
# TTA    R      G     yy
# 1: 3.0 safe before 0.2514
# 2: 3.5 safe before 0.5984
# 3: 4.0 safe before 1.3000

idx <- which( is.na(dtmp$y) )
dtmp2 <- dtmp[!idx, ]
dtmp3 <- dtmp2[TTA != 2.5 & G == "before" & R == "safe"]
dtmp3$TTA <- factor(dtmp3$TTA, levels = c(3, 3.5, 4))

rvar <- "s"
wvar <- "TTA"
aez0 <- aov_ez(rvar, "y", data = dtmp3,  within = wvar, 
               anova_table = list(correction = "none", es = "none"))
aez0 <- aov_ez(rvar, "y", data = dtmp3,  within = wvar, 
               anova_table = list(es = "pes"))
aez0 <- aov_ez(rvar, "y", data = dtmp3,  within = wvar, 
               anova_table = list(es = "ges"))
# Effect       df  MSE         F  ges p.value
# TTA       2, 46 0.44 15.47 ***        <.001
# TTA 1.78, 41.04 0.50 15.47 *** .402   <.001 # pes
# TTA 1.78, 41.04 0.50 15.47 *** .276   <.001 # ges

dtmp3$s <- factor(dtmp3$s)
bf <- BayesFactor::anovaBF(y ~ TTA + s, data = dtmp3, whichRandom = rvar)
# [1] TTA + s : 11530 ±3.7%
r1 <- emmeans(aez0, ~TTA)
# TTA  emmean     SE df lower.CL upper.CL
# X3    0.251 0.1536 23  -0.0663    0.569
# X3.5  0.598 0.0783 23   0.4365    0.760
# X4    1.300 0.1885 23   0.9100    1.690
contrast(r1, interaction = "pairwise", adjust = "bonferroni")
# TTA_pairwise estimate    SE df t.ratio p.value
# X3 - X3.5      -0.347 0.176 23  -1.974  0.1816
# X3 - X4        -1.048 0.223 23  -4.704  0.0003
# X3.5 - X4      -0.702 0.173 23  -4.047  0.0015
# dtmp3[s != 22 & TTA %in% c(3)]$y

col1 <- dtmp3[TTA %in% c(3)]$y
col2 <- dtmp3[TTA %in% c(3.5)]$y
bftmp1 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)

col1 <- dtmp3[TTA %in% c(3)]$y
col2 <- dtmp3[TTA %in% c(4)]$y
bftmp2 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)

col1 <- dtmp3[TTA %in% c(3.5)]$y
col2 <- dtmp3[TTA %in% c(4)]$y
bftmp3 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)
extractBF(bftmp1)
extractBF(bftmp2)
extractBF(bftmp3)
#                 bf     error  
# Alt., r=0.707 1.12   0.0002311 
# Alt., r=0.707 277   3.234e-10 
# Alt., r=0.707 64.49 1.127e-06 
