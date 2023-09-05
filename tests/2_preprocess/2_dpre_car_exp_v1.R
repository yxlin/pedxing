##############################################################################80
## Multi-phaase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
# install.packages('car')
# install.packages('afex')
# install.packages('effsize')
# install.packages('fBasic')
# install.packages('infer')
# install.packages('emmeans')
pkg <- c('data.table', 'plotrix', 'car','afex', 'effsize', 'fBasics', 'infer',
         'emmeans', 'BayesFactor','ggbeeswarm', 'tibble')
sapply(pkg, require, character.only = TRUE)
rm (list = ls())
win_path <- "D:/03_Projects/pedxing"
unix_path <- "/media/yslin/kea/03_Projects/pedxing"
wk <- ifelse(.Platform$OS.type == "windows",
    shortPathName(win_path), unix_path
)

setwd(wk)
options(digits = 4)

cb <- Manu::get_pal("Kaka")
# source("R/0_functions/clearup_data_v2.R")

# Pre-car | safe -----------------------------
load("tests/5_analyseBayes/data/3_subjects_v0.RData")
source("tests/0_functions/utils_v1.R")
data(ped)

ddat <- dat[E=="exp"]
ddat <- ddat[s != 2]

length(unique(ddat$s))
length(unique(x1$s))

dsim <- x1[type == "Simulation"]
dtmp <- ddat[G == "before" & C == TRUE, .N, .(s, G, TTA)]

setorder(dtmp, TTA)
range(dtmp[TTA %in% c(3) & s != 2]$N)
range(dtmp[TTA %in% c(3.5) & s != 2]$N)
range(dtmp[TTA %in% c(4) & s != 2]$N)
dtmp[TTA == 3]
dtmp[TTA == 3.5]
dtmp[TTA == 4]

table(ddat$R)
table(dsim$R)

dsim$C <- dsim$R
dsim$R <- ifelse(dsim$C == TRUE, "safe", "hit")
dsim$TTANme <- factor(dsim$TTA)
dsim$CNme <- ifelse(dsim$C == TRUE, "O", "X")
dsim$G <- ifelse(dsim$GBool == TRUE, "before", "after")
# table(dsim$C)

ddat$type <- "data"
dsim$type <- "Simulation"
d <- rbind(ddat[, c("s", "RT", "R", "TTA", "G", "C", "TTANme",
                    "type")],
           dsim[, c("s", "RT", "R", "TTA", "G", "C", "TTANme",
                    "type")])
d$type <- factor(d$type)

rvar <- 's'
dv <- 'y'
wvar0 <- 'TTANme'
wvar1 <- c('TTANme', 'G')

ddat$GBool <- ifelse(ddat$G == "before", 1, 0)
dsim$GBool <- ifelse(dsim$G == "before", 1, 0)

## Pre-car proportion given safe----------
dP1 <- ddat[R == "safe", .(y = mean(GBool)), .(TTANme, s)] ## Pre-car proportion
dP2 <- dsim[R == "safe", .(y = mean(GBool)), .(TTANme, s)] ## Pre-car proportion
aez1 <- aov_ez(rvar, dv, data = dP1, within = wvar0)
aez2 <- aov_ez(rvar, dv, data = dP2, within = wvar0)
dpcar_dat <- afex_plot(aez1, x = wvar0, error = 'within', return = "data")
dpcar_sim <- afex_plot(aez2, x = wvar0, error = 'within', return = "data")

ddat2 <- replace_TTA(dpcar_dat, type = "Data")
dsim2 <- replace_TTA(dpcar_sim, type = "Simulation")
dmean <- rbind(ddat2$means, dsim2$means)
ddata <- rbind(ddat2$data, dsim2$data)

dmean$x <- factor(dmean$x, levels = c("2.5", "2.5\nSimulation", 
                                      "3", "3\nSimulation", 
                                      "3.5", "3.5\nSimulation",
                                      "4", "4\nSimulation"))
ddata$x <- factor(ddata$x, levels = c("2.5", "2.5\nSimulation", 
                                      "3", "3\nSimulation", 
                                      "3.5", "3.5\nSimulation",
                                      "4", "4\nSimulation"))


# p3 <- ggplot() +
#   geom_point(data = dmean, aes(x = x, y = y, colour = type), 
#              size = 5) +
#   geom_point(data = ddata, aes(x = x, y = y, colour = type), 
#              size = 1, shape = 1) +
#   geom_errorbar(data = dmean, width = 0,
#                 aes(x = x, ymin = lower.CL, ymax = upper.CL, colour = type)) +
#   scale_colour_discrete(type = cb) +
#   theme_minimal(base_size = 20)
# p3


# save(dmean, ddata, dpcar_dat, dpcar_sim, 
#   file = "tests/2_preprocess/data/BP_givenSafe_exp.RData")

save(dmean, ddata, dpcar_dat, dpcar_sim, 
  file = "tests/extdata/anova/precar_prop_given_safe_exp.RData")


# Safe | made before -----------------------------
dP3 <- ddat[G == "before", .(y = mean(C)), .(TTANme, G, s)]   ## Safe proportion
dP4 <- dsim[G == "before", .(y = mean(C)), .(TTANme, G, s)]   ## Safe proportion
aez7 <- aov_ez(rvar, dv, data = dP3, within = wvar0)
aez8 <- aov_ez(rvar, dv, data = dP4, within = wvar0)
dSafe_dat <- afex_plot(aez7, x = wvar0, error = 'within', return = "data")
dSafe_sim <- afex_plot(aez8, x = wvar0, error = 'within', return = "data")

ddat2 <- replace_TTA(dSafe_dat, type = "Data")
dsim2 <- replace_TTA(dSafe_sim, type = "Simulation")
dmean <- rbind(ddat2$means, dsim2$means)
ddata <- rbind(ddat2$data, dsim2$data)

dmean$x <- factor(dmean$x, levels = c("2.5", "2.5\nSimulation", 
                                      "3", "3\nSimulation", 
                                      "3.5", "3.5\nSimulation",
                                      "4", "4\nSimulation"))
ddata$x <- factor(ddata$x, levels = c("2.5", "2.5\nSimulation", 
                                      "3", "3\nSimulation", 
                                      "3.5", "3.5\nSimulation",
                                      "4", "4\nSimulation"))

cb <- Manu::get_pal("Kaka")

p3 <- ggplot() +
  geom_point(data = dmean, aes(x = x, y = y, colour = type), 
             size = 5) +
  geom_point(data = ddata, aes(x = x, y = y, colour = type), 
             size = 1, shape = 1) +
  geom_errorbar(data = dmean, width = 0,
                aes(x = x, ymin = lower.CL, ymax = upper.CL, colour = type)) +
  scale_colour_discrete(type = cb) +
  # facet_wrap(.~G, scale = "free_y") +
  theme_minimal(base_size = 20)
p3

p4 <- ggplot() +
  geom_point(data = dSafe_dat$means, aes(x = TTANme, y = y), colour = cb[1], 
             size = 5, fill = cb[1]) +
  geom_point(data = dSafe_sim$means, aes(x = TTANme, y = y), colour = cb[2], 
              size = 5, fill = cb[2]) +
  # geom_point(data = dRTO_dat$data, aes(x = TTANme, y = y), colour = cb[1], 
  #             size = 1, shape = 1) +
  # geom_point(data = dRTO_sim$data, aes(x = TTANme, y = y), colour = cb[2], 
  #            size = 1, shape = 1) +
  geom_errorbar(data = dSafe_dat$means, colour = cb[1], width = 0,
                 aes(x = TTANme, ymin = lower.CL, ymax = upper.CL)) +
  geom_errorbar(data = dSafe_sim$means, colour = cb[2], width = 0,
                 aes(x = TTANme, ymin = lower.CL, ymax = upper.CL)) 
  # facet_wrap(.~G, scale = "free_y")
p4  

save(dmean, ddata, dSafe_dat, dSafe_sim, 
file = "tests/extdata/anova/safe_prop_given_precar_exp.RData")

# Safe | made after -----------------------------
dP3 <- ddat[G == "after", .(y = mean(C)), .(TTANme, G, s)]   ## Safe proportion
dP4 <- dsim[G == "after", .(y = mean(C)), .(TTANme, G, s)]   ## Safe proportion
aez7 <- aov_ez(rvar, dv, data = dP3, within = wvar0)
aez8 <- aov_ez(rvar, dv, data = dP4, within = wvar0)
dSafe_dat <- afex_plot(aez7, x = wvar0, error = 'within', return = "data")
dSafe_sim <- afex_plot(aez8, x = wvar0, error = 'within', return = "data")

ddat2 <- replace_TTA(dSafe_dat, type = "Data")
dsim2 <- replace_TTA(dSafe_sim, type = "Simulation")
dmean <- rbind(ddat2$means, dsim2$means)
ddata <- rbind(ddat2$data, dsim2$data)

dmean$x <- factor(dmean$x, levels = c("2.5", "2.5\nSimulation", 
                                      "3", "3\nSimulation", 
                                      "3.5", "3.5\nSimulation",
                                      "4", "4\nSimulation"))
ddata$x <- factor(ddata$x, levels = c("2.5", "2.5\nSimulation", 
                                      "3", "3\nSimulation", 
                                      "3.5", "3.5\nSimulation",
                                      "4", "4\nSimulation"))
save(dmean, ddata, dSafe_dat, dSafe_sim, 
file = "tests/extdata/anova/safe_prop_given_postcar_exp.RData")
