##############################################################################80
## Multi-phrase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
pkg <- c('data.table', 'plotrix', 'car','afex', 'effsize', 'fBasics', 'infer',
         'emmeans', 'BayesFactor','ggbeeswarm', 'tibble', 'bayestestR')
sapply(pkg, require, character.only = TRUE)
rm (list = ls())
wk <- ifelse(.Platform$OS.type == "windows", shortPathName("E:/Documents/pedxing/"),
             '/media/yslin/Tui/Documents/pedxing/')

setwd(wk)
options(digits = 4)
load("data/pedxing.RData")
source("R/0_functions/clearup_data_v2.R")
d <- dat
as_tibble(dat)

sapply(d[, .(s, G, Side, R, CNme, TTANme)], unique)
sapply(d[, .(s, G, Side, R, CNme, TTANme)], table)
rvar <- 's'
dv <- 'y'

## LR factor has no influence ------------------------------------------
wvar0 <- c('Side')
wvar1 <- c('Side', 'G')

dRT <- d[, .(y = median(RT)), .(TTANme, G, Side, s)]
dP0 <- d[, .(y = mean(GBool)), .(TTANme, Side, s)] ## Pre-car proportion
dP1 <- d[, .(y = mean(C)), .(TTANme, G, Side, s)] ## Safe proportion

dRT0 <- d[, .(y = median(RT)), .(Side, s)]
dRT1 <- d[, .(y = median(RT)), .(Side, G, s)]
dPPreC <- d[, .(y = mean(GBool)), .(Side, s)] ## Pre-car proportion
dPSafe0 <- d[, .(y = mean(C)), .(Side, s)]    ## Safe proportion
dPSafe1 <- d[, .(y = mean(C)), .(Side, G, s)] ## Safe proportion

bf0 <- anovaBF(y ~ Side + G + s, data = dRT1, whichRandom = "s") 
# Bayes factor analysis
# --------------
# [1] Side + s              : 0.2776    ±4.87%
# [2] G + s                 : 1.265e+46 ±1.18%
# [3] Side + G + s          : 3.511e+45 ±1.81%
# [4] Side + G + Side:G + s : 1.291e+45 ±3.12%
#
# [1] Side + s : 1.244 ±4.69%
# [1] Side + s : 0.4529 ±0.92%
bf1 <- anovaBF(y ~ Side + s, data = dRT1[G=="before"], whichRandom = "s") 
bf2 <- anovaBF(y ~ Side + s, data = dRT1[G=="after"], whichRandom = "s") 

bf3 <- anovaBF(y ~ Side + G + s, data = dPSafe1, whichRandom = "s") 
bf4 <- anovaBF(y ~ Side + s, data = dPSafe1[G=="before"], whichRandom = "s") 
bf5 <- anovaBF(y ~ Side + s, data = dPSafe1[G=="after"], whichRandom = "s") 
# Bayes factor analysis
# --------------
# [1] Side + s              : 0.276     ±0.85%
# [2] G + s                 : 1.02e+15  ±1.35%
# [3] Side + G + s          : 3.357e+14 ±1.71%
# [4] Side + G + Side:G + s : 2.037e+14 ±4.41%
# [1] Side + s : 0.5241 ±0.9%
# [1] Side + s : 1.091 ±2.35%
# Against denominator:
#   y ~ s 

bf6 <- anovaBF(y ~ Side + s, data = dPPreC, whichRandom = "s") 
# [1] Side + s : 1.038 ±0.79%


insight::get_parameters(bf1)
hdi(bf1, ci = c(.80, .90, .95))
## Traditional ANOVA ---------------------------------------------------
aez0 <- aov_ez(rvar, dv, data = dRT0,  within = wvar0)
aez1 <- aov_ez(rvar, dv, data = dRT1,  within = wvar1)
aez2 <- aov_ez(rvar, dv, data = dPPreC,  within = wvar0)
aez3 <- aov_ez(rvar, dv, data = dPSafe0,  within = wvar0)
aez4 <- aov_ez(rvar, dv, data = dPSafe1,  within = wvar1)
# nice(aez0, correction = "GG", es = "pes")
# nice(aez1, correction = "GG", es = "pes")
# nice(aez2, correction = "GG", es = "pes")
# nice(aez3, correction = "GG", es = "pes")
# nice(aez4, correction = "GG", es = "pes")

## Effect       df    MSE           F   ges p.value
##    Side   1, 13   0.15        0.59  .043    .458  across pre-and-post car RT 
##
##     RTs separating pre and post car responses
##  Effect       df    MSE           F   ges p.value
##    Side   1, 13   0.00        2.16  .142    .166
##       G   1, 13   0.09     1320.51  .990   <.001*** ## shorter RT in pre car responses
##  Side:G   1, 13   0.00        0.09  .007    .770
## 
##    Side    1, 13   0.00        3.40  .207    .088+   ## Prop. of pre-car response 
##    Side    1, 13   0.00        0.32  .024    .581    ## Prop. of safe  response 
##
## Prop. of safe  response separating pre and post car responses
##  Effect       df    MSE           F   ges p.value
##    Side    1, 13   0.00       0.54   .040    .475
##       G    1, 13   0.00     133.13   .911   <.001 *** more safe post-car passing
##  Side:G    1, 13   0.00       1.33   .093    .270
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘+’ 0.1 ‘ ’ 1

## Right-side scene resulted in more pre-car response 
r2 <- emmeans(aez2, ~ Side)
contrast(r2, interaction = "pairwise", adjust = "bonferroni")
# Side  emmean     SE df lower.CL upper.CL
# left   0.510 0.0484 13    0.406    0.615
# right  0.541 0.0558 13    0.420    0.661
# Confidence level used: 0.95 
# Side_pairwise estimate     SE df t.ratio p.value
# left - right   -0.0302 0.0164 13  -1.843  0.0883

## BF disputed the ANOVA finding.
dPPreC <- d[, .(y = mean(GBool)), .(Side, s)] ## Pre-car proportion
bf2 <- anovaBF(y ~ Side + s, data = dPPreC, whichRandom = "s") 
## [1] Side + s : 1.056 ±0.87% the marginal effect of pre-car response prop. is
## not strong
# help("get_parameters", package = "insight")
# insight::get_parameters(bf2)
hdi(bf2, ci = c(.80, .90, .95))
# Highest Density Interval
# 
# Parameter     |        80% HDI |        90% HDI |        95% HDI
# -------------------------------------------------------------
#   mu         | [ 0.46,  0.59] | [ 0.44,  0.61] | [ 0.42,  0.63]
