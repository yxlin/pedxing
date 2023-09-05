##############################################################################80
## Multi-phase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
# install.packages("BayesFactor")
pkg <- c('data.table', 'plotrix', 'car','afex', 'effsize', 'fBasics', 'infer',
         'emmeans', 'BayesFactor',  'tibble')
# 'ggbeeswarm',
sapply(pkg, require, character.only = TRUE)
win_path <- 'D:/03_Projects/pedxing'
rm(list = ls())
wk <- ifelse(.Platform$OS.type == "windows", 
             shortPathName("E:/Documents/pedxing/"),
             '/media/yslin/Tui/Documents/pedxing/')

setwd(wk)
# options(digits = 4)
load("data/pedxing.RData")
source("R/0_functions/clearup_data_v2.R")
source("R/0_functions/utils_v1.R")
load("R/5_analyseBayes/data/3_subjects_v0.RData")

ddat <- dat
dsim <- x1[type == "Simulation"]

dsim$C <- dsim$R
dsim$R <- ifelse(dsim$C == TRUE, "safe", "hit")
dsim$TTANme <- factor(dsim$TTA)
dsim$CNme <- ifelse(dsim$C == TRUE, "O", "X")
dsim$G <- ifelse(dsim$GBool == TRUE, "before", "after")

ddat$type <- "data"
dsim$type <- "Simulation"

d <- rbind(ddat[, c("s", "RT", "R", "TTA", "G", "C", "GBool", "TTANme", "CNme", "type")],
           dsim[, c("s", "RT", "R", "TTA", "G", "C", "GBool", "TTANme", "CNme", "type")])
d$type <- factor(d$type)

rvar <- 's'
dv <- 'y'
wvar0 <- 'TTANme'
wvar1 <- c('TTANme', 'G')
## RT all-----------------------------
dRT0  <- ddat[, .(y = median(RT), n = length(RT)), .(TTANme, G, s)] 
dsRT0 <- dsim[, .(y = median(RT), n = length(RT)), .(TTANme, G, s)] 

dRT0  <- ddat[, .(y = median(RT)), .(TTANme, G, s)] ## RT safe and hit
dsRT0 <- dsim[, .(y = median(RT)), .(TTANme, G, s)] ## RT safe and hit
aez3 <- aov_ez(rvar, dv, data = dRT0,  within = wvar1)
aez4 <- aov_ez(rvar, dv, data = dsRT0,  within = wvar1)

dRT0_dat <- afex_plot(aez3, x = wvar1, error = 'within', return = "data")
dRT0_sim <- afex_plot(aez4, x = wvar1, error = 'within', return = "data")

ddat2 <- replace_TTA(dRT0_dat, type = "Data")
dsim2 <- replace_TTA(dRT0_sim, type = "Simulation")
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
p0 <- ggplot() +
  geom_point(data = dmean, aes(x = x, y = y, colour = type), 
             size = 5) +
  geom_point(data = ddata, aes(x = x, y = y, colour = type), 
             size = 1, shape = 1) +
  geom_errorbar(data = dmean, width = 0,
                aes(x = x, ymin = lower.CL, ymax = upper.CL, colour = type)) +
  scale_colour_discrete(type = cb) +
  facet_wrap(.~G, scale = "free_y") +
  theme_minimal(base_size = 20)
p0

p1 <- ggplot() +
  geom_point(data = dRT0_dat$means, aes(x = TTANme, y = y), colour = cb[1], 
             size = 5, fill = cb[1]) +
  geom_point(data = dRT0_sim$means, aes(x = TTANme, y = y), colour = cb[2], 
             size = 5, fill = cb[2]) +
  geom_point(data = dRT0_dat$data, aes(x = TTANme, y = y), colour = cb[1], 
             size = 1, shape = 1) +
  geom_point(data = dRT0_sim$data, aes(x = TTANme, y = y), colour = cb[2], 
             size = 1, shape = 1) +
  geom_errorbar(data = dRT0_dat$means, colour = cb[1], width = 0,
                aes(x = TTANme, ymin = lower.CL, ymax = upper.CL)) +
  geom_errorbar(data = dRT0_sim$means, colour = cb[2], width = 0,
                aes(x = TTANme, ymin = lower.CL, ymax = upper.CL)) +
  facet_wrap(.~G, scale = "free_y")
p1  

save(dmean, ddata, dRT0_dat, dRT0_sim, file = "R/2_preprocess/data/dRT0.RData")

