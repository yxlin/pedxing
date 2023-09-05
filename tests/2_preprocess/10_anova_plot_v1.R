##############################################################################80
## Multi-phase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
rm(list = ls())
pkg <- c('data.table', 'plotrix', 'car','afex', 'effsize', 'fBasics', 'infer',
         'emmeans', 'BayesFactor','ggbeeswarm', 'tibble')
sapply(pkg, require, character.only = TRUE)
win_path <- "D:/03_Projects/pedxing"
unix_path <- "/media/yslin/kea/03_Projects/pedxing"
wk <- ifelse(.Platform$OS.type == "windows",
    shortPathName(win_path), unix_path
)

setwd(wk)
options(digits = 4)
aov_theme <- theme_minimal(base_size = 24) +
  theme(legend.position = leg_pos,
        legend.key.size =  unit(1, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        legend.title = element_blank(),
        axis.title.y = element_blank())

## Plot all three ---------------------
cb <- Manu::get_pal("Kaka")
# save(dmean, ddata, dpcar_dat, dpcar_sim, file = "R/2_preprocess/data/dpcar.RData")
# save(dmean, ddata, dSafe_dat, dSafe_sim, file = "R/2_preprocess/data/dSafe.RData")
# save(dmean, ddata, dRT0_dat, dRT0_sim, file = "R/2_preprocess/data/dRT0.RData")
# save(dmean, ddata, dRTO_dat, dRTO_sim, file = "R/2_preprocess/data/dRTO.RData")

# load("tests/2_preprocess/data/archive/dpcar.RData"); dmcar <- dmean; dcar <- ddata
# load("tests/2_preprocess/data/archive/dSafe.RData"); dmsafe <- dmean; dsafe <- ddata

load('tests/extdata/anova/precar_prop_given_safe_exp.RData'); dmcar <- dmean; dcar <- ddata
load('tests/extdata/anova/safe_prop_given_precar_exp.RData'); dmsafe <- dmean; dsafe <- ddata
load('tests/extdata/anova/safe_prop_given_postcar_exp.RData'); dmsafe_post <- dmean; dsafe_post <- ddata
load('tests/extdata/anova/safe_RT_given_precar_exp.RData'); dmrt <- dmean; drt <- ddata
load('tests/extdata/anova/safe_RT_given_postcar_exp.RData'); dmrt_post <- dmean; drt_post <- ddata
load('tests/extdata/anova/unsafe_RT_given_precar_exp.RData'); dm_ert_post <- dmean; d_ert_post <- ddata

# load("tests/2_preprocess/data/archive/dRT0.RData"); dmrt <- dmean; drt <- ddata
# load("R/2_preprocess/data/archive/dRTOX.RData"); dmrt <- dmean; drt <- ddata
# tibble::as_tibble(dmcar)
# tibble::as_tibble(dmrt)
# tibble::as_tibble(dmsafe)

# dmrt$panel   <- "CRT (s)"
dmcar$panel <- "Pre-car\nproportion | safe"
dmsafe$panel <- "Safe proportion |\nmade before" 
dmsafe_post$panel <- "Safe proportion |\nmade after"

p1 <- ggplot() +
  geom_line(data = dmean, aes(x = TTA, y = y, group = gp, colour = type),
            size = line_size) +
  # geom_errorbar(data = dmean, width = 0.1, position = pd, size = line_size,
  #                 aes(x = TTA, ymin = lower.CL, ymax = upper.CL, colour = type)) +
  
  # geom_point(data = d, aes(x = TTA, y = y, colour = type), alpha = .5, 
  #            size = dpoint_size, shape = 1) +
  # geom_point(data = dmean, aes(x = TTA, y = y, colour = type, shape = type), 
  #            size = mpoint_size) +
  # scale_colour_discrete(type = cb) +
  # scale_x_discrete(name = "TTA", labels = lab) +
  # ylab("") +
  # facet_wrap(panelG~., scale = "free_y", ncol = 3) +
  aov_theme

# png(filename = "R/2_preprocess/figs/6_anovaplot_v2.png" , 800, 600)
# png(filename = "docs/manu/figs/6_anovaplot_v2.png" , 800, 600)
fig_path <- paste0(win_path, '/tests/figs/6_anovaplot_v2.pdf')
pdf(fig_path)
p1


dmrt$panel   <- "RT (s)"
dmcar$G <- ""

drt$panel   <- "RT (s)"
dsafe$panel <- "Safe prop." 
dcar$panel <- "Pre-car prop."
dcar$G <- ""
dmcols <- c("G", "y", "SE", "df", "lower.CL", "upper.CL", "x", "TTA", "type", "panel")
dcols <- c("G", "s", "x", "y", "TTA", "type", "panel")

dmean <- rbind(dmrt[, dmcols], dmsafe[, dmcols], dmcar[, dmcols])
d <- rbind(drt[, dcols], dsafe[, dcols], dcar[, dcols])

dmean$panelG <- factor( paste0(dmean$panel, "  ", dmean$G) ) 
d$panelG <- factor( paste0(d$panel, "  ", d$G) ) 

lab <- c("2.5", "", "3", "", "3.5", "", "4", "")

tibble(dmean)
dmean$gp <- paste0(dmean$type, dmean$panelG)
d$gp <- paste0(d$type, d$panelG)
pd <- position_dodge(0.2) # move them .1 to the left and right

mpoint_size <- 4
dpoint_size <- 2
line_size <- .5
leg_pos <- c(.85, .2)

p1 <- ggplot() +
  geom_line(data = dmean, aes(x = TTA, y = y, group = gp, colour = type),
            size = line_size) +
  geom_errorbar(data = dmean, width = 0.1, position = pd, size = line_size,
                  aes(x = TTA, ymin = lower.CL, ymax = upper.CL, colour = type)) +
  
  geom_point(data = d, aes(x = TTA, y = y, colour = type), alpha = .5, 
             size = dpoint_size, shape = 1) +
  geom_point(data = dmean, aes(x = TTA, y = y, colour = type, shape = type), 
             size = mpoint_size) +
  scale_colour_discrete(type = cb) +
  # scale_x_discrete(name = "TTA", labels = lab) +
  ylab("") +
  facet_wrap(panelG~., scale = "free_y", ncol = 3) +
  aov_theme

# png(filename = "R/2_preprocess/figs/6_anovaplot_v2.png" , 800, 600)
# png(filename = "docs/manu/figs/6_anovaplot_v2.png" , 800, 600)
fig_path <- paste0(win_path, '/tests/figs/6_anovaplot_v2.pdf')
pdf(fig_path)
p1
dev.off()

