##############################################################################80
## Multi-phase decision making 
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
             '/media/yslin/Tui/projects/pedxing/')
setwd(wk)
options(digits = 4)
require(pedxing)
cb <- Manu::get_pal("Kaka")
data(ped)

# source("R/0_functions/clearup_data_v2.R")
# source("R/0_functions/utils_v1.R")
# load("R/5_analyseBayes/data/3_subjects_v0.RData")

ddat <- dat
# dsim <- x1[type == "Simulation"]
# dsim$C <- dsim$R
# dsim$R <- ifelse(dsim$C == TRUE, "safe", "hit")
# dsim$TTANme <- factor(dsim$TTA)
# dsim$CNme <- ifelse(dsim$C == TRUE, "O", "X")
# dsim$G <- ifelse(dsim$GBool == TRUE, "before", "after")



ddat$type <- "Data"
dsim$type <- "Simulation"
d <- rbind(ddat[, c("s", "RT", "R", "TTA", "G", "C", "GBool", "TTANme", "CNme", "type")],
           dsim[, c("s", "RT", "R", "TTA", "G", "C", "GBool", "TTANme", "CNme", "type")])
d$type <- factor(d$type)

# sapply(d[, .(s, G, R, CNme, TTANme, type)], unique)
# sapply(d[, .(s, G, R, CNme, TTANme, type)], table)

rvar <- 's'
wvar <- 'TTANme'
dv <- 'y'

dO <- ddat[C==TRUE & E == "eeg"]
dO$GBool <- ifelse(dO$G == "before", 1, 0)
## TTA affects pre-car prop. -------------------------------------------
#### Pre-car proportion was affected by the TTA ----------
dP0 <- ddat[, .(y = mean(GBool)), .(TTANme, s)] ## Pre-car proportion
dP1 <- dO[, .(y = mean(GBool)), .(TTANme, s)]   ## Pre-car proportion
aez1 <- aov_ez(rvar, dv, data = dP1, within = wvar)
aez0 <- aov_ez(rvar, dv, data = dP0, within = wvar)
nice(aez0, correction = "GG", es = "pes")
nice(aez1, correction = "GG", es = "pes")
# Effect          df  MSE         F  pes p.value
# TTANme 1.27, 16.53 0.03 33.54 *** .721   <.001 both correct and error RT
# TTANme 2.04, 26.55 0.03 78.14 *** .857   <.001 only correct RT

d[, .N, .(TTANme, G, C)]
d[, .(n = .N,
      M = mean(RT)), .(TTANme, G, C)]

r0 <- emmeans(aez0, ~ TTANme)
r1 <- emmeans(aez1, ~ TTANme)
# dP0[, .N, .(TTANme)]
# dP1[, .N, .(TTANme)]
# TTA   emmean     SE df lower.CL upper.CL
# 2.5    0.316 0.0613 13    0.183    0.448
# 3      0.445 0.0550 13    0.326    0.563
# 3.5    0.602 0.0584 13    0.476    0.728
# 4      0.740 0.0591 13    0.613    0.868
# contrast(r0, interaction = "pairwise", adjust = "bonferroni")
# TTA_pairwis   e estimate     SE df t.ratio p.value
# X2.5 - X3         -0.129 0.0286 13  -4.503  0.0036
# X2.5 - X3.5       -0.286 0.0564 13  -5.073  0.0013
# X2.5 - X4         -0.425 0.0666 13  -6.375  0.0001
# X3 - X3.5         -0.157 0.0359 13  -4.383  0.0044
# X3 - X4           -0.296 0.0443 13  -6.679  0.0001
# X3.5 - X4         -0.139 0.0233 13  -5.938  0.0003

# lim <- c(0.2242, 0.4072 )
p0 <- afex_plot(aez0, x = wvar, error = 'within')
p1 <- afex_plot(aez1, x = wvar, error = 'within')

p0 <- p0 + ylab('Prop. of pre-car responses') + xlab('TTA (s)') 
p1 <- p1 + ylab('Prop. of pre-car responses') + xlab('TTA (s)') 
# png(filename = "R/2_preprocess/figs/tmp1.png" , 800, 600)
gridExtra::grid.arrange(p0, p1, ncol = 2)
dev.off()


# bf0 <- anovaBF(y ~ TTANme + s, data = dP1, whichRandom = "s") 
# [1] TTA + s : 1.05 x  10^8 
# require(bayestestR)
# hdi(bf0, ci = c(.80, .90, .95))
# TTANme-2.5 | [-0.24, -0.17] | [-0.25, -0.15] | [-0.25, -0.14]
# TTANme-3   | [-0.11, -0.04] | [-0.12, -0.03] | [-0.13, -0.02]
# TTANme-3.5 | [ 0.04,  0.11] | [ 0.03,  0.12] | [ 0.02,  0.13]
# TTANme-4   | [ 0.17,  0.24] | [ 0.16,  0.25] | [ 0.15,  0.26]

dpre_car_tmp <- afex_plot(aez0, x = 'TTANme', error = 'within', return = "data")
dpre_carO_tmp <- afex_plot(aez1, x = 'TTANme', error = 'within', return = "data")
dpre_car <- replace_TTA(dpre_car_tmp)
dpre_carO <- replace_TTA(dpre_carO_tmp)

dpre_car$mean$type <- "All"
dpre_carO$mean$type <- "Safe"
dpre_car$data$type <- "All"
dpre_carO$data$type <- "Safe"
dmean <- rbind(dpre_car$mean, dpre_carO$mean)
dd <- rbind(dpre_car$data, dpre_carO$data)

