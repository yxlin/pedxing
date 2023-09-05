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
data(ped)
# source("R/0_functions/clearup_data_v2.R")
# source("R/0_functions/utils_v1.R")
# load("R/2_preprocess/data/dsafeP_v1.RData"); 
cb <- Manu::get_pal("Kaka")

rvar <- 's'
wvar <- 'TTANme'
dv <- 'y'


ddat <- dat
## safe RTs EEG ------------
dO <- ddat[C==TRUE & E == "exp"]
daov <- dO[, .(y = median(RT)), .(TTANme, G, s)]
aez0 <- aov_ez(rvar, dv, data = daov[G=='after'],  within = wvar) 
aez1 <- aov_ez(rvar, dv, data = daov[G=='before'], within = wvar) 
nice(aez0, correction = "GG", es = "pes")
nice(aez1, correction = "GG", es = "pes")

r0 <- emmeans(aez0, ~ TTANme)
r1 <- emmeans(aez1, ~ TTANme)
contrast(r1, interaction = "pairwise", adjust = "bonferroni")

dtmp <- daov[G=='before']
tibble::as_tibble(dtmp)
dtmp$s <- factor(dtmp$s)
bf <- BayesFactor::anovaBF(y ~ TTANme + s, data = dtmp, whichRandom = rvar)
extractBF(bf)

col1 <- dtmp[TTANme %in% c(3)]$y
col2 <- dtmp[TTANme %in% c(3.5)]$y
bftmp1 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)
col1 <- dtmp[TTANme %in% c(3)]$y
col2 <- dtmp[TTANme %in% c(4)]$y
bftmp2 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)
col1 <- dtmp[TTANme %in% c(3.5)]$y
col2 <- dtmp[TTANme %in% c(4)]$y
bftmp3 <- BayesFactor::ttestBF(col1, col2, paired = TRUE)
extractBF(bftmp1)
extractBF(bftmp2)
extractBF(bftmp3)

# Effect          df  MSE           F  pes p.value
# TTANme 1.56, 35.93 0.00 26.05 *** .531   <.001
# TTANme_pairwise estimate      SE df t.ratio p.value
# X3 - X3.5       -0.05825 0.01122 23  -5.192  0.0001
# X3 - X4         -0.06745 0.01157 23  -5.832  <.0001
# X3.5 - X4       -0.00921 0.00697 23  -1.322  0.5979


## unsafe RTs exp ------------
dX <- dat[R == "hit" & E == "exp" & G == "before"]

daov <- dX[s != 2, .(y = median(RT)), .(TTANme, s)]
tibble::as_tibble(daov)
daov$s <- factor(daov$s)
levels(daov$s)
unique(daov$s)
# aez0 <- aov_ez(rvar, dv, data = daov, within = wvar)
# dtmp <- daov[TTANme %in% c(2.5, 3)]
# sort(unique(dtmp$s))
# sort(unique(dat[E == "exp"]$s))

aez0 <- aov_ez(rvar, dv, data = daov[TTANme %in% c(2.5, 3)], within = wvar) 
nice(aez0, correction = "GG", es = "pes")
r0 <- emmeans(aez0, ~ TTANme)
contrast(r0, interaction = "pairwise", adjust = "bonferroni")
# Effect   df  MSE         F  pes p.value
# TTANme_pairwise estimate     SE df t.ratio p.value
# X2.5 - X3         -0.292 0.0304 12  -9.617  <.0001

daov[, .N, .(TTANme)]

tibble::as_tibble(daov)
dtmp <- daov[TTANme %in% c(2.5, 3) & s != 7]
dtmp$s <- factor(dtmp$s)
bf <- BayesFactor::anovaBF(y ~ TTANme + s, data = dtmp, whichRandom = rvar)
tmp <- extractBF(bf)
tmp$bf
# TTANme emmean     SE df lower.CL upper.CL
# X2.5    0.600 0.0297 12    0.536    0.665
# X3      0.892 0.0214 12    0.846    0.939


## unsafe RTs eeg ------------
dX <- dat[R == "hit" & E == "eeg" & G == "before"]
daov <- dX[, .(y = median(RT)), .(TTANme, s)]
tibble::as_tibble(daov)
daov$s <- factor(daov$s)
aez0 <- aov_ez(rvar, dv, data = daov[TTANme %in% c(2.5, 3)], within = wvar) 
nice(aez0, correction = "GG", es = "pes")
r0 <- emmeans(aez0, ~ TTANme)
contrast(r0, interaction = "pairwise", adjust = "bonferroni")
# Effect    df  MSE      F  pes p.value
# TTANme 1, 21 0.13 7.94 * .274    .010
# TTANme emmean     SE df lower.CL upper.CL
# X2.5    0.652 0.0213 21    0.608    0.697
# X3      0.953 0.1001 21    0.745    1.161

# 17, 24
tibble::as_tibble(daov)
dtmp <- daov[TTANme %in% c(2.5, 3) & s != 17 & s !=24]
dtmp$s <- factor(dtmp$s)
bf <- BayesFactor::anovaBF(y ~ TTANme + s, data = dtmp, whichRandom = rvar)
tmp <- extractBF(bf)
tmp$bf

## RM-ANOVA collapsing over correct and error RTs  ------------
daov <- dat[, .(y = median(RT)), .(TTANme, G, s)]
aez0 <- aov_ez(rvar, dv, data = daov[G=='after'], within = wvar) 
aez1 <- aov_ez(rvar, dv, data = daov[G=='before'], within = wvar) 
## after ---------------------------------45
nice(aez0, correction = "GG", es = "pes")
nice(aez1, correction = "GG", es = "pes")
# Anova Table (Type 3 tests): median RTs
# Effect          df  MSE           F  pes p.value
# af-TTA 2.18, 26.18 0.01 1949.57 *** .994   <.001
# be-TTA 1.60, 20.83 0.01 10.47 ** .446    .001

r0 <- emmeans(aez0, ~ TTANme)
r1 <- emmeans(aez1, ~ TTANme)
bf0 <- anovaBF(y ~ TTANme + s, data = daov[G=='after'], whichRandom = rvar)
bf1 <- anovaBF(y ~ TTANme + s, data = daov[G=='before'], whichRandom = rvar)
# [1] TTANme + s : 1.399e+42 ±0.5%
# [1] TTANme + s : 639 ±0.48%

daov_safe <- dO[, .(y = median(RT)), .(TTANme, G, s)]
aez2 <- aov_ez(rvar, dv, data = daov_safe[G=='after'], within = wvar) 
aez3 <- aov_ez(rvar, dv, data = daov_safe[G=='before'], within = wvar) 
bf2 <- anovaBF(y ~ TTANme + s, data = daov_safe[G=='after'], whichRandom = rvar)
bf3 <- anovaBF(y ~ TTANme + s, data = daov_safe[G=='before'], whichRandom = rvar)
# [1] TTANme + s : 1.488e+42 ±0.58%
# [1] TTANme + s : 2.43 ±1.11%


tmp0 <- afex_plot(aez0, x = wvar, error = 'within', return = "data")
tmp1 <- afex_plot(aez1, x = wvar, error = 'within', return = "data")
dOaf <- replace_TTA(tmp0)
dObe <- replace_TTA(tmp1)

dOaf$mean$G <- "After"
dObe$mean$G <- "Before"
dOaf$data$G <- "After"
dObe$data$G <- "Before"
dmean <- rbind(dObe$mean, dOaf$mean)
dd <- rbind(dObe$data, dOaf$data)

dmean_RT <- dmean
dd_RT <- dd

dmean_pcar$G <- "before"
dmean_safeP$type <- "Safe"
dmean_RT$type <- "All"

dmean_pcar$DV <- "Pre-car prop."
dmean_safeP$DV <- "Safe prop."
dmean_RT$DV <- "RT (s)"

dmean <- rbind(dmean_pcar, dmean_safeP, dmean_RT)

mpoint_size <- 4
dpoint_size <- 2
line_size <- 1.25
leg_pos <- c(.1, .9)
dpctheme <- theme_minimal(base_size = 20) +
  theme(legend.position = leg_pos,
        legend.title = element_blank(),
        strip.text.x = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.title.x = element_blank())

p0 <- ggplot() +
    geom_errorbar(data = dmean_pcar, size = line_size,
                  aes(x = TTA, ymin = lower, ymax = upper, colour = type), width = 0) +
    # geom_line(data = dpre_car$data, aes(x = TTA, y = y, group = s),
    #           colour = "grey80") +
    geom_line(data = dmean_pcar, aes(x = TTA, y = y, colour = type), size = line_size) +
    geom_point(data = dd_pcar, aes(x = TTA, y = y, colour = type), shape = 2) +
    geom_point(data = dmean_pcar, aes(x = TTA, y = y, colour = type), size = 3.5,
             shape = 21, fill = "white") + xlab("TTA") + ylab("") + dpctheme
p0

p1 <- ggplot() +
  geom_errorbar(data = dmean_safeP, size = line_size, 
                aes(x = TTA, ymin = lower, ymax = upper), width = 0) +
  geom_point(data = dd_safeP, aes(x = TTA, y = y), size = dpoint_size, 
              shape = 2) +
  geom_point(data = dmean_safeP, aes(x = TTA, y = y), size = mpoint_size, 
             shape = 21, fill = "white") + xlab("TTA") +
  ylab("") +
  facet_wrap(.~G, scales ="free_y") +
  dpctheme

p1

gridExtra::grid.arrange(p0, p1, ncol = 2)

p1 <- ggplot() +
  # geom_point(data=dsp, aes(x = TTA, y = y, shape=G), colour = "grey", size = 2) +
  geom_errorbar(data = dmean, aes(x = TTA, ymin = lower, ymax = upper),
                 width = 0)  +
  # geom_errorbar(data = dmean, aes(x = TTA, ymin = lower.CL, ymax = upper.CL), 
  #               width = 0)  +
  
  # geom_line(data = dsp, aes(x = TTA, y = y, group = Gs, colour = G), linetype = "dotted") +
  # geom_line(data = ds, aes(x = TTA, y = y, group = G), size = 1) +
  geom_point(data = dmean, aes(x = TTA, y = y), size = 3.5, 
             shape = 21, fill = "white") + xlab("TTA") +
  ylab("RT (s)") +
  facet_wrap(DV~type, scales ="free_y") 
# + dpctheme
p1

save(dmean_pcar, dd_pcar, dmean_safeP, dd_safeP,
     ddat, dsim, file = "R/2_preprocess/data/dsafeP_v1.RData")

