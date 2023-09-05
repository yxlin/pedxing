##############################################################################80
## Multi-phase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
# install.packages("gamlss.dist")
# install.packages("ggbeeswarm")
pkg <- c('data.table',  "facetscales", 'ggbeeswarm')
sapply(pkg, require, character.only = TRUE)
rm (list = ls())
win_path <- "D:/03_Projects/pedxing"
unix_path <- "/media/yslin/kea/03_Projects/pedxing"
wk <- ifelse(.Platform$OS.type == "windows",
    shortPathName(win_path), unix_path
)

setwd(wk)
options(digits = 4)
rm(list = ls())

## RT --------------------------------------------------------
# d0[s %in% c(9, 13, 17, 18)]
# load("tests/inst/extdata/exp/pedxing.RData")
# source("tests/0_functions/clearup_data_v2.R")
# source("tests/0_functions/utils_v1.R")
data("ped")

dtmp0 <- d0[d0$sNme == "Bhix" & d0$Trial == 28 & d0$TTA == 2.5,]
dtmp1 <- d0[d0$sNme == "VFFt" & d0$Trial == 17 & d0$TTA == 2.5 & d0$Delay == 3, ]
dtmp2 <- d0[d0$sNme == "Bhix" & d0$Trial == 27 & d0$TTA == 2.5,]
dtmp3 <- d0[d0$sNme == "Bhix" & d0$Trial == 24 & d0$TTA == 2.5,]
outlier1 <- which(d0$sNme == "Bhix" & d0$Trial == 28 & d0$TTA == 2.5)
outlier2 <- which(d0$sNme == "VFFt" & d0$Trial == 17 & d0$TTA == 2.5 & d0$Delay == 3)
outlier3 <- which(d0$sNme == "Bhix" & d0$Trial == 27 & d0$TTA == 2.5)
outlier4 <- which(d0$sNme == "Bhix" & d0$Trial == 24 & d0$TTA == 2.5)
dexp0 <- d0[E=="exp"]
dexp1 <- dexp0[!c(outlier1,outlier2,outlier3,outlier4)]
dexp <- dexp1[E=="exp"]
tibble::as_tibble(dexp)

deeg <- d0[E=="eeg"]
length(unique(deeg$s))


# Separate proportion data
dunsafeExp <- dexp[, .(Avg_TTA = mean(TTA + Jitter),
                    Sd_TTA = sd(TTA + Jitter)), .(TTA)]
dunsafeEeg <- deeg[, .(Avg_TTA = mean(TTA + Jitter),
                       Sd_TTA = sd(TTA + Jitter)), .(TTA)]

getRTdf <- function(x, cols = c("x","y","G","TTA", "CNme"), binwidth = 8/85, rng = c(0, 6)) {
  # x <- sim
  # binwidth = 8/85
  # rng = c(0, 6)
  bk <- seq(rng[1], rng[2], by = binwidth)
  dA <- x[G=="after"]
  dB <- x[G=="before"]
  rt1 <- dA$RT
  rt0 <- dB$RT

  y0 <- cut(rt0, breaks = bk)
  y1 <- cut(rt1, breaks = bk)
  
  tmp0 <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", y0) ),
                upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", y0) ))
  tmp1 <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", y1) ),
                upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", y1) ))
  dB$mid <- .5*rowSums(tmp0)
  dA$mid <- .5*rowSums(tmp1)
  
  d1 <- dB[, .(y = .N / nrow(dB)), .(C, TTA, mid)]
  d2 <- dA[, .(y = .N / nrow(dA)), .(C, TTA, mid)]
  
  names(d1) <- c("C", "TTA", "x", "y")
  names(d2) <- c("C", "TTA", "x", "y")
  d1$CNme <- ifelse(d1$C, "safe", "unsafe")
  d2$CNme <- ifelse(d2$C, "safe", "unsafe")
  d1$G <- "before"
  d2$G <- "after"
  
  out <- rbind(d1[,..cols], d2[,..cols])
  
  return(out)
}

d12 <- getRTdf(dexp)
d12e <- getRTdf(deeg)


dvlineExp <- dexp[, .(x = median(RT),
                      upper = median(RT) + plotrix::std.error(RT),
                      lower = median(RT) - plotrix::std.error(RT)
                      ), .(TTA, G, C)]

dvlineEeg <- deeg[, .(x = median(RT),
                      upper = median(RT) + plotrix::std.error(RT),
                      lower = median(RT) - plotrix::std.error(RT)
), .(TTA, G, C)]

dvlineExp$CNme <- ifelse(dvlineExp$C, "safe", "unsafe")
dvlineEeg$CNme <- ifelse(dvlineEeg$C, "safe", "unsafe")

dvlineExp$CE <- paste0(dvlineExp$CNme)
dvlineEeg$CE <- paste0(dvlineEeg$CNme)

# dvlineExp$CE <- paste0(dvlineExp$CNme, "-E1")
# dvlineEeg$CE <- paste0(dvlineEeg$CNme, "-E2")
dvlineExp$E <- "E1"
dvlineEeg$E <- "E2"
dvline_all <- rbind(dvlineExp, dvlineEeg)

# mean(dvline[G == "before" & C == TRUE & TTA >= 3.0]$x)

# require(scico)

cb1 <- Manu::get_pal("Kotare")
cb2 <- Manu::get_pal("Kaka")
# cb <- scico(6, palette = "devon")

dgood <- d12[which(!is.na(d12$x))]
dgoode <- d12e[which(!is.na(d12e$x))]

dgood$CE <- paste0(dgood$CNme, "-E1")
dgoode$CE <- paste0(dgoode$CNme, "-E2")
dgood$E <- "E1"
dgoode$E <- "E2"
dall <- rbind(dgood, dgoode)

p0 <- ggplot(data = dall[E=='E1']) +
  geom_line( aes(x = x, y = y, group = CNme, colour = CNme), size = .7) +
  geom_point( aes(x = x, y = y, colour = CNme), shape = 1, size = 2) +
  
  geom_point(data = dvline_all[E=='E1'], aes(x = x, y = 0.075, colour = CNme)) +
  geom_segment(data = dunsafeExp, aes(x = Avg_TTA - Sd_TTA, xend = Avg_TTA + Sd_TTA,
                                     y = 0, yend = 0), size = 1) +
  scale_colour_manual( values = cb1 ) +
  xlab("RT (s)") +  ylab("Proportion of trials") +
  facet_grid(TTA~.) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.9, .8), 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y  = element_text(size = 10))

# cairo_pdf("/media/yslin/Tui/projects/pedxing/tests/docs/2nd_submit/figs/2_RTDist_v7.pdf")
fig_path <- paste0(win_path, '/tests/figs/2_RTDist_exp.pdf')
pdf(fig_path)
p0
dev.off()

p1 <- ggplot(data = dall[E=='E2']) +
  geom_line( aes(x = x, y = y, group = CNme, colour = CNme), size = .7) +
  geom_point( aes(x = x, y = y, colour = CNme), shape = 1, size = 2) +
  
  geom_point(data = dvline_all[E=='E2'], aes(x = x, y = 0.075, colour = CNme)) +
  geom_segment(data = dunsafeExp, aes(x = Avg_TTA - Sd_TTA, xend = Avg_TTA + Sd_TTA,
                                     y = 0, yend = 0), size = 1) +
  scale_colour_manual( values = cb1 ) +
  xlab("RT (s)") +  ylab("Proportion of trials") +
  facet_grid(TTA~.) +
  theme_classic(base_size = 20) +
  theme(legend.position = c(.9, .8), 
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        axis.text.y  = element_text(size = 10))

# cairo_pdf("/media/yslin/Tui/projects/pedxing/tests/docs/2nd_submit/figs/2_RTDist_v7.pdf")
fig_path <- paste0(win_path, '/tests/figs/2_RTDist_eeg.pdf')
pdf(fig_path)
p1
dev.off()


# dall[G == "after" & E == "E1"]
# d25 <- dall[G == "after" & E == "E1" & TTA == 2.5]
# d30 <- dall[G == "after" & E == "E1" & TTA == 3.0]
# d35 <- dall[G == "after" & E == "E1" & TTA == 3.5]
# d40 <- dall[G == "after" & E == "E1" & TTA == 4]
# setorder(d25, x)
# setorder(d30, x)
# setorder(d35, x)
# setorder(d40, x)
# head(d25)
# head(d30)
# head(d35)
# head(d40)

dtmp <- dexp[G=="after" & RT <= 8, .(qRT = quantile(RT, prob = seq(.1))), .(s, TTA)]
cor.test(dtmp$qRT, dtmp$TTA)
psycho::bayes_cor.test(dtmp$qRT, dtmp$TTA)
# Results of the Bayesian correlation indicate extreme evidence (BF > 100) in favour of 
# the existence of a positive association between dtmp$qRT and dtmp$TTA (r = 0.48, MAD = 0.11, 90% CI [0.27, 0.67]). The correlation can be considered as large, moderate, small or very small with respective probabilities of 44.55%, 48.98%, 6.31% and 0.08%.
# res1 <- BayesianFirstAid::bayes.cor.test( ~ x + y, data = dtmp)
BayesFactor::correlationBF(dtmp$qRT, dtmp$TTA)


tibble::as_tibble(dtmp)
setorder(dtmp, TTA)
dtmp

dtmp <- deeg[G=="after", .(qRT = quantile(RT, prob = seq(.1))), .(s, TTA)]
cor.test(dtmp$qRT, dtmp$TTA)
BayesFactor::correlationBF(dtmp$qRT, dtmp$TTA)
psycho::bayes_cor.test(dtmp$qRT, dtmp$TTA)
setorder(dtmp, TTA)
dtmp


# res0 <- psycho::bayes_cor.test(dtmp$x, dtmp$y)
# res1 <- BayesianFirstAid::bayes.cor.test( ~ x + y, data = dtmp)
