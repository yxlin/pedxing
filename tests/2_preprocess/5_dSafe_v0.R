##############################################################################80
## Multi-phrase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
pkg <- c('data.table', 'plotrix', 'car','afex', 'effsize', 'fBasics', 'infer',
         'emmeans', 'BayesFactor','ggbeeswarm', 'tibble')
install.packages("ggbeeswarm")
sapply(pkg, require, character.only = TRUE)
rm (list = ls())
wk <- ifelse(.Platform$OS.type == "windows", shortPathName("E:/Documents/pedxing/"),
             '/media/yslin/Tui/Documents/pedxing/')

setwd(wk)
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
# tibble::as_tibble(d)

rvar <- 's'
dv <- 'y'
wvar0 <- 'TTANme'
wvar1 <- c('TTANme', 'G')
cb <- Manu::get_pal("Kaka")
## Proportion of safe -----------------------------
# dP3 <- ddat[, .(y = mean(C),
#                 n = length(C)), .(TTANme, G)]   ## Safe proportion
# dP4 <- dsim[, .(y = mean(C),
#                 n = length(C)), .(TTANme, G)]   ## Safe proportion
# length(unique(ddat$s))
# length(unique(dsim$s))

dP3 <- ddat[, .(y = mean(C)), .(TTANme, G, s)]   ## Safe proportion
dP4 <- dsim[, .(y = mean(C)), .(TTANme, G, s)]   ## Safe proportion
aez7 <- aov_ez(rvar, dv, data = dP3, within = wvar1)
aez8 <- aov_ez(rvar, dv, data = dP4, within = wvar1)
dSafe_dat <- afex_plot(aez7, x = wvar1, error = 'within', return = "data")
dSafe_sim <- afex_plot(aez8, x = wvar1, error = 'within', return = "data")

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
  facet_wrap(.~G, scale = "free_y") +
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
                 aes(x = TTANme, ymin = lower.CL, ymax = upper.CL)) +
  facet_wrap(.~G, scale = "free_y")
p4  

# save(dmean, ddata, dSafe_dat, dSafe_sim, file = "R/2_preprocess/data/dSafe.RData")

## Plot all three ---------------------
rm(list = ls())
cb <- Manu::get_pal("Kaka")
# save(dmean, ddata, dpcar_dat, dpcar_sim, file = "R/2_preprocess/data/dpcar.RData")

# save(dmean, ddata, dSafe_dat, dSafe_sim, file = "R/2_preprocess/data/dSafe.RData")
# save(dmean, ddata, dRT0_dat, dRT0_sim, file = "R/2_preprocess/data/dRT0.RData")
# save(dmean, ddata, dRTO_dat, dRTO_sim, file = "R/2_preprocess/data/dRTO.RData")

load("R/2_preprocess/data/dpcar.RData"); dmcar <- dmean; dcar <- ddata
load("R/2_preprocess/data/dSafe.RData"); dmsafe <- dmean; dsafe <- ddata
load("R/2_preprocess/data/dRT0.RData"); dmrt <- dmean; drt <- ddata
# load("R/2_preprocess/data/dRTO.RData"); dmrt <- dmean; drt <- ddata
# tibble::as_tibble(dmcar)
# tibble::as_tibble(dmrt)
# tibble::as_tibble(dmsafe)

dmrt$panel   <- "RT (s)"
dmsafe$panel <- "Safe prop." 
dmcar$panel <- "Pre-car prop."
dmcar$G <- ""

drt$panel   <- "RT (s)"
dsafe$panel <- "Safe prop." 
dcar$panel <- "Pre-car prop."
dcar$G <- ""
dmcols <- c("G", "y", "SE", "df", "lower.CL", "upper.CL", "x", "TTA", "type", "panel")
dcols <- c("G", "s", "x", "y", "TTA", "type", "panel")
dmean <- rbind(dmrt[, dmcols], dmsafe[, dmcols], dmcar[, dmcols])
d <- rbind(drt[, dcols], dsafe[, dcols], dcar[, dcols])

tibble::as_tibble(dmean)
tibble::as_tibble(d)

dmean$panelG <- factor( paste0(dmean$panel, "  ", dmean$G) ) 
d$panelG <- factor( paste0(d$panel, "  ", d$G) ) 

lab <- c("2.5", "", "3", "", "3.5", "", "4", "")
p0 <- ggplot() +
  geom_point(data = dmean, aes(x = x, y = y, colour = type), size = 2) +
  geom_point(data = d, aes(x = x, y = y, colour = type), alpha = .5, size = 1, shape = 1) +
  geom_errorbar(data = dmean, width = 0,
                aes(x = x, ymin = lower.CL, ymax = upper.CL, colour = type)) +
  scale_colour_discrete(type = cb) +
  scale_x_discrete(name = "TTA", labels = lab) +
  ylab("") +
  facet_wrap(panelG~., scale = "free_y", ncol = 3) +
  theme_minimal(base_size = 24) +
  theme(legend.position = c(.85, .2), ## "none",
        legend.key.size =  unit(1, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        # strip.text.x = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank())

# png(filename = "R/2_preprocess/figs/5_dSafe_facet_v0.png" , 800, 600)
p0
dev.off()
  
  
