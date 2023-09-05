pkg <- c("data.table", "ggplot2", "rstan", "rethinking", "dplyr", "pedxing")
sapply(pkg, require, character.only = TRUE)
win_path <- "D:/03_Projects/pedxing"
unix_path <- "/media/yslin/kea/03_Projects/pedxing"
wk <- ifelse(.Platform$OS.type == "windows",
    shortPathName(win_path), unix_path
)

setwd(wk)
options(digits = 4)
rm(list = ls())
data(ped)

## Analysis ---------------------
load("tests/extdata/mid-stage/precar_given_safe_eeg.RData")
summary(m0@stanfit)$summary

load("tests/extdata/mid-stage/precar_given_safe_eeg.RData")
ms0 <- compare(m0, m1, m2, m3, m4)
load("tests/extdata/mid-stage/safe_given_precar_eeg.RData")
ms1 <- compare(m0, m1, m2, m3, m4)
load("tests/extdata/mid-stage/safe_given_postcar_eeg.RData")
ms2 <- compare(m0, m1, m2, m3, m4)
load("tests/extdata/mid-stage/safe_rt_given_precar_eeg.RData")
ms3 <- compare(m0, m1, m2, m3, m4)
load("tests/extdata/mid-stage/safe_rt_given_postcar_eeg.RData")
ms4 <- compare(m0, m1, m2, m3, m4)
load("tests/extdata/mid-stage/unsafe_rt_given_precar_eeg.RData")
ms5 <- compare(m0, m1, m2, m3, m4)


dms0 <- data.table(ms0@output)
dms1 <- data.table(ms1@output)
dms2 <- data.table(ms2@output)
dms3 <- data.table(ms3@output)
dms4 <- data.table(ms4@output)
dms5 <- data.table(ms5@output)

dms0$M <- "(a) Pre-car \nproportion | safe"
dms1$M <- "(b) Safe proportion | \nmade before"
dms2$M <- "(c) Safe proportion | \nmade after"
dms3$M <- "(d) RT | made before"
dms4$M <- "(e) RT | made after"
dms5$M <- "(f) RT | made before, \ncollision"

dms0$m <- factor(rownames(ms0@output),
    levels = rownames(ms0@output)[order(dms0$WAIC)]
)
dms1$m <- factor(rownames(ms1@output),
    levels = rownames(ms1@output)[order(dms1$WAIC)]
)
dms2$m <- factor(rownames(ms2@output),
    levels = rownames(ms2@output)[order(dms2$WAIC)]
)
dms3$m <- factor(rownames(ms3@output),
    levels = rownames(ms3@output)[order(dms3$WAIC)]
)
dms4$m <- factor(rownames(ms4@output),
    levels = rownames(ms4@output)[order(dms4$WAIC)]
)
dms5$m <- factor(rownames(ms5@output),
    levels = rownames(ms5@output)[order(dms5$WAIC)]
)
dms <- rbind(dms0, dms1, dms2, dms3, dms4, dms5)

Models <- c("(a) Pre-car \nproportion | safe", "(b) Safe proportion | \nmade before", "(c) Safe proportion | \nmade after", "(d) RT | made before", "(e) RT | made after", "(f) RT | made before, \ncollision")
dms$M <- factor(dms$M, levels = Models)
p <- ggplot(data = dms) +
    geom_point(aes(x = m, y = WAIC)) +
    geom_errorbar(aes(x = m, y = WAIC, ymin = WAIC + SE, ymax = WAIC - SE),
        width = .5
    ) +
    facet_wrap(M ~ ., scales = "free") +
    xlab("Model") +
    ylab("Deviance") +
    theme_bw(base_size = 18) +
    theme(strip.text.x = element_text(size = 10))


# options(bitmapType="cairo")
cairo_pdf("tests/figs/glm_eeg_waic.pdf")
print(p)
dev.off()

## Simulation ---------------------
# Pre-car Proportion | safe crossing ------------
load("tests/extdata/mid-stage/precar_given_safe_eeg.RData")

ns <- sum(rownames(summary(m0@stanfit)$summary) %in% paste0("a[", 1:24, "]"))
dp0 <- data.table(
    s      = rep(1:ns, each = 8),
    TTAint = rep(c(0, 1, 2, 3, 0, 1, 2, 3), ns),
    TTANme = rep(c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4), ns),
    total  = rep(32, 8 * ns)
)
dp0$TTANme <- factor(dp0$TTANme, levels = c(2.5, 3, 3.5, 4))

m0pred <- ensemble(m0, data = dp0)
m1pred <- ensemble(m1, data = dp0)

# summarize
pred0.p <- apply(m0pred$link, 2, mean)
pred1.p <- apply(m1pred$link, 2, mean)

pred0.p.PI <- apply(m0pred$link, 2, PI)
pred1.p.PI <- apply(m1pred$link, 2, PI)

dp0$m0value <- pred0.p
dp0$m1value <- pred1.p
dp0$m005 <- pred0.p.PI[1, ]
dp0$m094 <- pred0.p.PI[2, ]
dp0$m105 <- pred1.p.PI[1, ]
dp0$m194 <- pred1.p.PI[2, ]

dpred0 <- dp0[, .(
    y = mean(m0value),
    ymin = mean(m005),
    ymax = mean(m094)
), .(TTAint)]
dpred1 <- dp0[, .(
    y = mean(m1value),
    ymin = mean(m105),
    ymax = mean(m194)
), .(TTAint)]

dpred0$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))
dpred1$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))

## Safe prop | before----------------
load("tests/extdata/mid-stage/safe_given_precar_eeg.RData")
d.pred <- data.table(
    s      = rep(1:ns, each = 8),
    TTAint = rep(c(0, 1, 2, 3, 0, 1, 2, 3), ns),
    TTANme = rep(c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4), ns),
    total  = rep(32, 8 * ns)
)
d.pred$TTANme <- factor(d.pred$TTANme, levels = c(2.5, 3, 3.5, 4))

m0pred <- ensemble(m0, data = d.pred)
m1pred <- ensemble(m1, data = d.pred)

# summarize
pred0.p <- apply(m0pred$link, 2, mean)
pred1.p <- apply(m1pred$link, 2, mean)
d.pred$m0value <- pred0.p
d.pred$m1value <- pred1.p

pred0.p.PI <- apply(m0pred$link, 2, PI)
pred1.p.PI <- apply(m1pred$link, 2, PI)
d.pred$m005 <- pred0.p.PI[1, ]
d.pred$m094 <- pred0.p.PI[2, ]
d.pred$m105 <- pred1.p.PI[1, ]
d.pred$m194 <- pred1.p.PI[2, ]

dpred2 <- d.pred[, .(
    y = mean(m0value),
    ymin = mean(m005),
    ymax = mean(m094)
), .(TTAint)]
dpred3 <- d.pred[, .(
    y = mean(m1value),
    ymin = mean(m105),
    ymax = mean(m194)
), .(TTAint)]

dpred2$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))
dpred3$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))

## Safe prop | after ----------------------------
load("tests/extdata/mid-stage/safe_given_postcar_eeg.RData")
d.pred <- data.table(
    s      = rep(1:ns, each = 8),
    TTAint = rep(c(0, 1, 2, 3, 0, 1, 2, 3), ns),
    TTANme = rep(c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4), ns),
    total  = rep(32, 8 * ns)
)
d.pred$TTANme <- factor(d.pred$TTANme, levels = c(2.5, 3, 3.5, 4))

m0pred <- ensemble(m0, data = d.pred)
m1pred <- ensemble(m1, data = d.pred)

# summarize
pred0.p <- apply(m0pred$link, 2, mean)
pred1.p <- apply(m1pred$link, 2, mean)
d.pred$m0value <- pred0.p
d.pred$m1value <- pred1.p

pred0.p.PI <- apply(m0pred$link, 2, PI)
pred1.p.PI <- apply(m1pred$link, 2, PI)
d.pred$m005 <- pred0.p.PI[1, ]
d.pred$m094 <- pred0.p.PI[2, ]
d.pred$m105 <- pred1.p.PI[1, ]
d.pred$m194 <- pred1.p.PI[2, ]

dpred4 <- d.pred[, .(
    y = mean(m0value),
    ymin = mean(m005),
    ymax = mean(m094)
), .(TTAint)]
dpred5 <- d.pred[, .(
    y = mean(m1value),
    ymin = mean(m105),
    ymax = mean(m194)
), .(TTAint)]
dpred4$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))
dpred5$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))

## RT | before ----------------------------
load("tests/extdata/mid-stage/safe_rt_given_precar_eeg.RData")
d.pred <- data.table(
    s      = rep(1:ns, each = 8),
    TTAint = rep(c(0, 1, 2, 3, 0, 1, 2, 3), ns),
    TTANme = rep(c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4), ns),
    total  = rep(32, 8 * ns)
)
d.pred$TTANme <- factor(d.pred$TTANme, levels = c(2.5, 3, 3.5, 4))

m0pred <- ensemble(m0, data = d.pred)
m1pred <- ensemble(m1, data = d.pred)

# summarize
pred0.p <- apply(m0pred$link, 2, mean)
pred1.p <- apply(m1pred$link, 2, mean)
d.pred$m0value <- pred0.p
d.pred$m1value <- pred1.p

pred0.p.PI <- apply(m0pred$link, 2, PI)
pred1.p.PI <- apply(m1pred$link, 2, PI)
d.pred$m005 <- pred0.p.PI[1, ]
d.pred$m094 <- pred0.p.PI[2, ]
d.pred$m105 <- pred1.p.PI[1, ]
d.pred$m194 <- pred1.p.PI[2, ]

dpred6 <- d.pred[, .(
    y = mean(m0value),
    ymin = mean(m005),
    ymax = mean(m094)
), .(TTAint)]
dpred7 <- d.pred[, .(
    y = mean(m1value),
    ymin = mean(m105),
    ymax = mean(m194)
), .(TTAint)]
dpred6$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))
dpred7$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))

## RT | after  ----------------------------
load("tests/extdata/mid-stage/safe_rt_given_postcar_eeg.RData")
d.pred <- data.table(
    s      = rep(1:ns, each = 8),
    TTAint = rep(c(0, 1, 2, 3, 0, 1, 2, 3), ns),
    TTANme = rep(c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4), ns),
    total  = rep(32, 8 * ns)
)
d.pred$TTANme <- factor(d.pred$TTANme, levels = c(2.5, 3, 3.5, 4))


m0pred <- ensemble(m0, data = d.pred)
m1pred <- ensemble(m1, data = d.pred)

# summarize
pred0.p <- apply(m0pred$link, 2, mean)
pred1.p <- apply(m1pred$link, 2, mean)
d.pred$m0value <- pred0.p
d.pred$m1value <- pred1.p

pred0.p.PI <- apply(m0pred$link, 2, PI)
pred1.p.PI <- apply(m1pred$link, 2, PI)
d.pred$m005 <- pred0.p.PI[1, ]
d.pred$m094 <- pred0.p.PI[2, ]
d.pred$m105 <- pred1.p.PI[1, ]
d.pred$m194 <- pred1.p.PI[2, ]

dpred8 <- d.pred[, .(
    y = mean(m0value),
    ymin = mean(m005),
    ymax = mean(m094)
), .(TTAint)]
dpred9 <- d.pred[, .(
    y = mean(m1value),
    ymin = mean(m105),
    ymax = mean(m194)
), .(TTAint)]
dpred8$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))
dpred9$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))

## RT | before, collision  ----------------------------
load("tests/extdata/mid-stage/unsafe_rt_given_precar_eeg.RData")
d.pred <- data.table(
    s      = rep(1:ns, each = 8),
    TTAint = rep(c(0, 1, 2, 3, 0, 1, 2, 3), ns),
    TTANme = rep(c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4), ns),
    total  = rep(32, 8 * ns)
)
d.pred$TTANme <- factor(d.pred$TTANme, levels = c(2.5, 3, 3.5, 4))


m0pred <- ensemble(m0, data = d.pred)
m1pred <- ensemble(m1, data = d.pred)

# summarize
pred0.p <- apply(m0pred$link, 2, mean)
pred1.p <- apply(m1pred$link, 2, mean)
d.pred$m0value <- pred0.p
d.pred$m1value <- pred1.p

pred0.p.PI <- apply(m0pred$link, 2, PI)
pred1.p.PI <- apply(m1pred$link, 2, PI)
d.pred$m005 <- pred0.p.PI[1, ]
d.pred$m094 <- pred0.p.PI[2, ]
d.pred$m105 <- pred1.p.PI[1, ]
d.pred$m194 <- pred1.p.PI[2, ]

dpred10 <- d.pred[, .(
    y = mean(m0value),
    ymin = mean(m005),
    ymax = mean(m094)
), .(TTAint)]
dpred11 <- d.pred[, .(
    y = mean(m1value),
    ymin = mean(m105),
    ymax = mean(m194)
), .(TTAint)]
dpred10$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))
dpred11$TTANme <- factor(c(2.5, 3, 3.5, 4), levels = c(2.5, 3, 3.5, 4))
save(dpred0, dpred1, dpred2, dpred3, dpred4, dpred5, dpred6, dpred7, dpred8, dpred9, dpred10, dpred11,
    file = "tests/extdata/mid-stage/3_glm_eeg_analysis.rda"
)

## Plot data and models ----------
dpred0$M <- "(a) Pre-car \nproportion | safe"
dpred1$M <- "(a) Pre-car \nproportion | safe"
dpred2$M <- "(b) Safe proportion | \nmade before"
dpred3$M <- "(b) Safe proportion | \nmade before"
dpred4$M <- "(c) Safe proportion | \nmade after"
dpred5$M <- "(c) Safe proportion | \nmade after"

dpred6$M <- "(d) RT | made before"
dpred7$M <- "(d) RT | made before"
dpred8$M <- "(e) RT | made after"
dpred9$M <- "(e) RT | made after"

dpred10$M <- "(f) RT | made before, \ncollision"
dpred11$M <- "(f) RT | made before, \ncollision"

Models <- c("(a) Pre-car \nproportion | safe", "(b) Safe proportion | \nmade before", "(c) Safe proportion | \nmade after", "(d) RT | made before", "(e) RT | made after", "(f) RT | made before, \ncollision")

dpredM0 <- rbind(dpred0, dpred2, dpred4, dpred6, dpred8, dpred10)
dpredM1 <- rbind(dpred1, dpred3, dpred5, dpred7, dpred9, dpred11)

dpredM0$M <- factor(dpredM0$M, levels = Models)
dpredM1$M <- factor(dpredM1$M, levels = Models)


## Data
load("tests/extdata/mid-stage/proportion_rt_eeg.rda")
cols <- c("TTANme", "s", "value", "M")
dp_be$M <- "(a) Pre-car \nproportion | safe"
dsafe_be$M <- "(b) Safe proportion | \nmade before"
dsafe_af$M <- "(c) Safe proportion | \nmade after"
dbe$M <- "(d) RT | made before"
daf$M <- "(e) RT | made after"
dbe_hit$M <- "(f) RT | made before, \ncollision"

dsum <- rbind(
    dp_be[, ..cols], dsafe_be[, ..cols], dsafe_af[, ..cols],
    dbe[, ..cols], daf[, ..cols], dbe_hit[, ..cols]
)

dsum_mean <- dsum[, .(
    n = length(value),
    value = mean(value)
), .(TTANme, M)]


dsum$M <- factor(dsum$M, levels = Models)
dsum_mean$M <- factor(dsum_mean$M, levels = Models)

jitter <- position_jitter(width = 0.1, height = 0.01)

# pd <- position_dodge(5) # move them .1 to the left and right

setorder(dsum_mean, TTANme)

p0 <- ggplot() +
    geom_point(
        data = dsum, aes(x = TTANme, y = value), alpha = .2,
        position = jitter
    ) +
    geom_point(
        data = dsum_mean, aes(x = TTANme, y = value),
        colour = "sienna", size = 6
    ) +
    geom_errorbar(
        data = dpredM0,
        aes(x = TTANme, y = y, ymin = ymin, ymax = ymax),
        colour = get_pal("Kaka")[1], width = .3, linewidth = 1
    ) +
    geom_errorbar(
        data = dpredM1, aes(x = TTANme, y = y, ymin = ymin, ymax = ymax),
        colour = get_pal("Kaka")[2], width = .3, linewidth = 1
    ) +
    geom_text(
        data = dsum_mean, size = 4, colour = "white",
        aes(x = TTANme, y = value, label = n)
    ) +
    xlab("TTA (s)") +
    ylab("Prop. or RT (s)") +
    facet_wrap(. ~ M, scales = "free") +
    theme_bw(base_size = 18) +
    theme(strip.text.x = element_text(size = 10))


cairo_pdf("tests/figs/glm_eeg_gof.pdf")
print(p0)
dev.off()
