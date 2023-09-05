#####################################################################80
## Multi-phase decision making 
## Version: 1.0
## Authors: Yi-Shin Lin (y.s.lin@leeds.ac.uk)
## Date: 2022, Jan, 28
## License: GPL
## Description: A hierarchical Bayesian UM model
# install.packages("plotrix")
rm(list = ls())
pkg <- c("data.table", "psych", "plotrix", "ggplot2")
sapply(pkg, require, character.only = TRUE)
win_path <- "D:/03_Projects/pedxing"
unix_path <- "/media/yslin/kea/03_Projects/pedxing"


wk <- ifelse(.Platform$OS.type == "windows",
             shortPathName(win_path), unix_path
)

setwd(wk)
options(digits = 4)


source("tests/0_functions/utils_v1.R")

# Sys.which("make")
# path.expand('~')
# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

## Step over the details -----------------------------------------------------
d0 <- fread("tests/inst/extdata/exp/pilot/ANT4/Commotions_Output_0.csv")
d1 <- fread("tests/inst/extdata/exp/pilot/ANT4/Commotions_Output_1.csv")
# nrow(d0)  ## 256 = 4 * 2 * 8 * 4 rows
# nrow(d1)  ## 256 
## 4 TTAs x 2 (left and right)  x 8 repetitions x 4 blocks
## If the side factor, left vs. right, has no influence, we will have 64 trials
## per cell. This was expected and this was confirmed later in Bayes factor 
## test.
256/4
## list.files("data/exp/pilot/ANT4/")

#### How many trial in each condition? --------------------
res0 <- d0[, .N, .(TTA, Side)]
res1 <- d1[, .N, .(TTA, Side)]
setorder(res0, TTA, Side)
setorder(res1, TTA, Side)
res0
res1

#### Did the jittering work? ------------------------------
DT0 <- data.table(x = rep(1:nrow(d0), 2), 
                 y = c(d0$TTA, d0$TTA + d0$Jitter),
                 group = rep(c("TTA", "J-TTA"), each = nrow(d0)))
envelop <- sort( c( unique(d0$TTA)-.1, unique(d0$TTA) + .1) )
range(d0$Jitter)
p0 <- ggplot(DT0) +
    geom_point(aes(x = x, y = y, color = group)) +
    geom_hline(yintercept = envelop, linetype = "dashed") +
    ylab("TIME (Unity second)") + xlab("Index") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "top")
p0

## Proportion ## 
d0$BeforeCarPassed <- (d0$CrossingkeyPressedTime - d0$TrialStartTime 
                       - d0$signalDelay) < (d0$TTA + d0$Jitter)

prop0 <- d0[, .N, .(TTA, BeforeCarPassed, Hit)]
prop0[, NN := sum(N), .(TTA)]
prop0[, Pcerent := round(N/NN, 2)]
setorder(prop0, TTA, BeforeCarPassed)
cat("Percentage of been hit\n")
as_tibble(prop0)
cat("Only trials passing before the vehicle\n")
prop0[BeforeCarPassed == TRUE]

## RT ## 
RT0 <- d0$CrossingkeyPressedTime - (d0$TrialStartTime + d0$signalDelay)
RT1 <- d0$CrossingkeyPressedTime - (d0$CarAppearingTim)

RT <- d0$CrossingkeyPressedTime - (d0$CarAppearingTim)

cat("Average imprecision in second\n")
all.equal(RT0, RT1)

DT <- data.table(RT = RT, 
                 R  = factor(ifelse(d0$Hit, "was hit", "safely crossed")),
                 TTA = d0$TTA,
                 group = ifelse(d0$BeforeCarPassed, "Crossed before", 
                                "Crossed after"))

DT$TTA <- factor(DT$TTA, levels = sort(unique(DT$TTA)))
binsize <- psych::describe(RT)$se  
bk <- seq(min(RT)-binsize, max(RT)+binsize, by = binsize)
DT$groupTTA <- paste0(DT$group, "-", DT$TTA)

# cb <- Manu::get_pal("Tui")
DT$Rgp <- paste0(DT$group, ", ", DT$R)
p1 <- ggplot(DT) +
    geom_histogram(aes(x = RT, fill = Rgp), breaks = bk) +
    # scale_fill_manual(values = cb) +
    xlab("RT (s)") +
    theme_minimal(base_size = 20) +
    theme(legend.position = c(.60, .90),
          legend.title = element_blank())

# png(file = "RTDist_ANT4_2.png", width = 480, height = 480)
p1
# dev.off()

# postscript(file = "preprocessing.ps", width = 480, height = 480)
# png(file = "preprocessing.png", width = 480, height = 480)
# dev.off()

options(digits = 2)
res <- DT[, .(N = .N,
              Mean = mean(RT),
       Mdn  = median(RT),
       Max  = max(RT),
       Min  = min(RT),
       SD   = sd(RT),
       SE   = plotrix::std.error(RT),
       skew = psych::skew(RT),
       kurtosis = psych::kurtosi(RT)
       ), .(TTA, R)]
res[order(res$TTA), ]

## Batch process ----------------------------------------------
sNme_exp <- data.frame(sq = 1:24,
                   s = c("P4jJ", "mNoy", "Dndk", "Bhix", "69WA", "wLu7",
                         "B0kn", "HLBQ", "wXWs", "QU2L","52ho", "wToD",
                         "YSlV", "VFFt", "amug", "qlNp", "LG2o", "OcQN",
                         "Jv1V", "zFpL", "PjiE", "5NNB", "OEWy", "mkKt"),
                   sex = c("female", "male", "male", "female", "male",
                           "female", "male", "male", "female", "male",
                           "male", "female", "female", "male", "female",
                           "female", "female", "female", "female", "male",
                           "female","male", "female", NA),
                   age = c(22, 41, 47, 38, 27,  48,
                           28, 22, 45, NA, 44, 40,
                           44, 26, 51, 48, 32, 22,
                           37, 64, 40, 38, 51, NA),
                   handedness = c('r', 'r', 'r', 'r', 'r', 'l',
                                  'r', 'r', 'r', NA, 'r', 'r',
                                  'r', 'r', 'r', 'r', 'r', 'l',
                                  'r', 'r', 'r', 'r', 'r', NA))


dropout_exp <- c(1, 8, 23, 24)
discard_exp <- c(3, 5, 6, 10, 16)
valid_exp <- c(2, 4, 7, 9, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22)
sNme_exp[dropout_exp, ]
sNme_exp[discard_exp, ]
sNme_exp[valid_exp, ]
length(valid_exp)   ## 15
length(discard_exp) ## 5
# 24 - 4 - 5
# 5/20
# 'pG2u' = 'wXWs'
# 'Bhix' = 'w3ug'
sNme_eeg <- data.frame(sq = 1:24,
                       s = c('pG2u','snv0','rnd0','NEd0','HBfP','xuzS', 
                             'EtJa','md0M','2uqO','4WWn','pdLD','lDHY',
                             '3maD','5TtL','7pve','12UD','w3ug','bcKj',
                             '55tV','I8BJ','T4Vy','JNvZ','ulp9','dDjo'),
                       sex = c("female", 'female', 'male', 'male', 'female','male',
                               'female', 'male', 'female', 'male', 'female','male',
                               'male', 'male', 'female', 'female', 'female','male',
                               'female', 'female', 'female', 'female', 'male', 'male'),
                       age = c(45, 22, 42, 29, 24, 60,
                               61, 19, 53, 19, 35, NA,
                               47, 28, 26, 28, 38, 19,
                               38, 21, 19, 62, 22, 41),
                       handedness = c('r', 'r', 'r', 'r', 'r', 'r',
                                      'r', 'r', 'r', 'r', 'r', 'r',
                                      'r', 'r', 'r', 'r', 'r', 'r',
                                      'r', 'r', 'r', 'r', 'r', 'r'))
table(sNme_exp$sex)
table(sNme_eeg$sex)

length( unique(sNme_exp$s) )
length( unique(sNme_eeg$s) )



discard_eeg <- NULL
# valid_eeg <- sNme_eeg$sq[-discard_eeg]
valid_eeg <- sNme_eeg$sq
c( length(valid_eeg), length(valid_exp) )
c( nrow(sNme_eeg), nrow(sNme_exp) )

dtmp0 <- sNme_exp[valid_exp, ]
dtmp1 <- sNme_eeg[valid_eeg, ]
dtmp0$exp <- 'exp'
dtmp1$exp <- 'eeg'
dbg <- data.table(rbind(dtmp0, dtmp1))

## Exp experiment
# dage <- tibble::as_tibble(dbg[exp == 'exp'])
# dage$age_gp <- cut(dage$age, breaks = c(18, 30, 40, 50, 60, 70), right = TRUE)
# dage[is.na(dage$age_gp),]
# table(dage$age_gp)
# psych::describe(dage$age)
# (18,30] (30,40] (40,50] (50,60] (60,70]
#       3       6       4       1       1
#
#        vars  n  mean    sd median trimmed  mad min max range skew kurtosis   se
# X1    1 15 39.33 10.33     40   38.77 5.93  22  64    42 0.46     0.12 2.67

# require(data.table)
# table(dbg$sex)

## EEG experiment
dage <- tibble::as_tibble(dbg[exp == 'eeg'])
dage$age_gp <- cut(dage$age, breaks = c(18, 30, 40, 50, 60, 70), right = TRUE)
dage[is.na(dage$age_gp),]

length(unique(dage$s))
table(dage$sex)
range(dage$age, na.rm = TRUE)
table(dage$age_gp)
# One NA did not fill in his age
# (18,30] (30,40] (40,50] (50,60] (60,70] 
#     12       3       4       2       2 


## Two experiments together; all participants
## length(valid_exp)  ## 21 participants in the final analysis
dage <- tibble::as_tibble(dbg)
dage$age_gp <- cut(dage$age, breaks = c(18, 30, 40, 50, 60, 80), right = TRUE)
table(dage$age_gp)
# (20,30] (30,40] (40,50] (50,60] (60,70] 
#     15       9       8       3       3 
# sum(table(dage$age_gp))

# range(dbg$age, na.rm = TRUE) 
# 19 64

# hist(dbg$age)
# table(dbg$sex)
# female   male 
#    22     17 
22+17


#### Load choice RT data ------------------------------------------------------
E <- c('exp', 'eeg')
path <- c("inst/extdata/exp/valid/", "inst/extdata/eeg/valid/") 
ns0 <- length(list.files(path[1])); ns0
ns1 <- length(list.files(path[2])); ns1

fn <- "Commotions_Output_1.csv"

## Participation sequence
##  [1] "mNoy" "Bhix" "B0kn" "wXWs" "52ho" "WtoD" "YSIV" "VFFt" "amug" "LG2o"
## [11] "OcQN" "Jv1V" "zFpL" "PjiE" "5NNB"
tmpNme <- c("Participant", "Trial", "TTA", "Jitter", "Side", "signalDelay",
            "TrialStartTime", "CarAppearingTime", "CrossingkeyPressedTime", 
            "Hit", "OtherKeyPressedTime", "FirstKeyPressed", "AvgFPS")
d0 <- NULL
# i <- 1
# j
for(i in seq_len(2)) {
    path <- paste0('inst/extdata/', E[i], '/valid/')    
    sNme <- dbg[exp == E[i]]$s
    ns <- length(sNme)
    
    for(j in seq_len(ns)) {
        dfile <- paste0(path, sNme[j],"/",fn); dfile
        output <- paste0(path, sNme[j], '/result_bk.ps')
        tmp <- default_result(fn = fn, path = paste0(path, sNme[j]),
                              output = output, verbose = FALSE)    
        
        dtmp <- fread(dfile)
        names(dtmp) <- tmpNme
    
        dtmp$s <- dbg[s == sNme[j]]$sq
        dtmp$age <- dbg[s == sNme[j]]$age
        dtmp$sex <- dbg[s == sNme[j]]$sex
        
        if(!dbg[s == sNme[j]]$s == sNme[j]) stop("Subject name mismatched")
        
        dtmp$sNme <- sNme[j]
        dtmp$exp <- E[i]
        
        d0 <- rbind(d0, dtmp[, 2:18])
    }
}    

## d0
names(d0) <- c("Trial", "TTA", "Jitter",  "Side", "Delay", "TrialStartTime", 
  "CarAppearingTime", "CrossingkeyPressedTime", "Hit", "OtherKeyPressedTime",
  "FirstKeyPressed", "AvgFPS",  "s", "age", "sex", "sNme", "E")
d0

#### Rename columns ---------------------------
d0$BeforeCarPassed <- (d0$CrossingkeyPressedTime - d0$CarAppearingTim) < (d0$TTA + d0$Jitter)
d0$RT <- d0$CrossingkeyPressedTime - d0$CarAppearingTim
d0$TTANme <- factor(d0$TTA)
d0$Side <- factor(d0$Side)

# d0$RT0 <- d0$CrossingkeyPressedTime - (d0$TrialStartTime + d0$signalDelay)
# d0$RT1 <- d0$CrossingkeyPressedTime - (d0$CarAppearingTim)

d0$R <- factor(ifelse(d0$Hit, "hit", "safe"))
d0$C <- as.logical(ifelse(d0$Hit, 0, 1))
# levels(d0$TTANme)

d0$G <- factor(ifelse(d0$BeforeCarPassed, "before", "after"))
d0$D <- 40 

dat <- d0[, c('Trial', 'Jitter',"BeforeCarPassed", "E", 'TTA', 'Side','R', 'C',
              'G', 'TTANme', 'RT', 's')]

# table(dat[, .N, .(s, E)]$N)
# table(dat[, .N, .(s, E, TTA)]$N)

## Add D and DTTA columns for legacy code
dat$D <- 40 
dat$Side <- factor(dat$Side)
dat$DTTA <- factor(paste0(dat$D, "-", dat$TTA), 
                   levels = c('40-2.5', '40-3', '40-3.5', '40-4'))


## Participants in exp and eeg were different groups, so use E and s to 
## differeniate them.
# load("tests/inst/extdata/exp/pedxing.RData")
valid_exp <- c(2, 4, 7, 9, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22)
# valid_eeg <- c(2:7, 9:24)

dexp <- dat[E == "exp" & s %in% valid_exp]
deeg <- dat[E == "eeg" & s %in% valid_eeg]

dexp$s <- factor(dexp$s, levels = valid_exp)
deeg$s <- factor(deeg$s, levels = valid_eeg)
dexp$snew <- paste0(dexp$s, dexp$E)
deeg$snew <- paste0(deeg$s, deeg$E)

d <- rbind(dexp, deeg)
d$TTANme <- factor(d$TTANme, levels = c(2.5, 3, 3.5, 4)) 

d$sold <- d$s
d$s <- as.integer( factor(d$snew) )
setorder(d, s, TTA, Side, G, R)

d$TTAint <- as.integer( as.integer(factor(d$TTA)) - 1 )
d$side   <- as.integer( as.integer(d$Side) - 1 )
d$TTAint <- as.integer( as.integer(factor(d$TTA)) - 1 )
d$side   <- as.integer( as.integer(d$Side) - 1 )



#### Save to pedxing.RData ---------
## dat reduces columns
## d0 has all the columns
## d for stan 
save(d0, dat, d, file = "../data/ped-tmp.rda")
tools::resaveRdaFiles("../data/ped-tmp.rda", compress = "xz",
                           compression_level = 9)
length(unique(deeg$s))
length(unique(dexp$s))

#### Two abnormal responses -------------------------------------------
## Road width: 4.2m
## Initial distance from Kerb: 0.5 m
## Car Length: 4.96 m
## Car Width: 1.9 m
## Car initial distance: 40 m
40/c(2.5, 3, 3.5, 4)

dtmp0 <- d0[d0$sNme == "Bhix" & d0$Trial == 28 & d0$TTA == 2.5,]
dtmp1 <- d0[d0$sNme == "VFFt" & d0$Trial == 17 & d0$TTA == 2.5 & d0$Delay == 3, ]
dtmp2 <- d0[d0$sNme == "Bhix" & d0$Trial == 27 & d0$TTA == 2.5,]
dtmp3 <- d0[d0$sNme == "Bhix" & d0$Trial == 24 & d0$TTA == 2.5,]

time2pass <- (4.2+.5) / 1.6  
car2pass <- 4.96 / 40
dtmp0$CarAppearingTime + dtmp0$TTA + dtmp0$Jitter
dtmp0$CrossingkeyPressedTime 

dtmp1$CarAppearingTime + dtmp1$TTA + dtmp1$Jitter
dtmp1$CrossingkeyPressedTime 

dtmp2$CarAppearingTime + dtmp2$TTA + dtmp2$Jitter
dtmp2$CrossingkeyPressedTime 

dtmp3$CarAppearingTime + dtmp3$TTA + dtmp3$Jitter
dtmp3$CrossingkeyPressedTime 

## Conclusion: According to the Unity programmer, Unity 
## represents the pedestrian and car by using cylinders. Therefore, imprecision 
## representation might be the reason causing these two abnormal responses. 

## RT distributions  ----------------------------------------------------------
# load("data/pedxing.RData")
# dplt <- dat[E == "exp" & RT <= 6 & RT > .1]
dplt <- dat[RT <= 6 & RT > .1]
1 - nrow(dplt) / nrow(dat)
dtmp <- dplt[, .N, .(TTA, G, R)]
setorder(dtmp, R, G, TTA)
dtmp[TTA == 2.5 & G == "before"]
dtmp[TTA == 3& G == "before"]
dtmp[TTA == 3.5& G == "before"]


# Simple get_break function
rt <- dplt$RT
binsize <- psych::describe(rt)$se  
bk <- seq(min(rt)-binsize, max(rt)+binsize, by = binsize)

# DT$groupTTA <- paste0(DT$group, "-", DT$TTA)
dv <- data.frame(TTA = c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4),
                 Etta = c("exp, 2.5", "exp, 3", "exp, 3.5", "exp, 4",
                          "eeg, 2.5", "eeg, 3", "eeg, 3.5", "eeg, 4"))

# sort(unique(DT$TTA))
# [1] 2.5 3.0 3.5 4.0 

dplt$Etta <- factor( paste0(dplt$E, ", ", dplt$TTA), 
                     levels = c("exp, 2.5", "exp, 3", "exp, 3.5", "exp, 4",
                                "eeg, 2.5", "eeg, 3", "eeg, 3.5", "eeg, 4"))

dplt$Rgp <- factor( paste0(dplt$G, ", ", dplt$R), 
                    levels = c("before, safe", "before, hit", "after, safe", "after, hit"))
table(dplt$Rgp)
pal <- Manu::get_pal("Tui")
danno <- data.table()

danno1 <- data.table(x = c(.15, .15), 
                    y = c(120, 120),
                    Etta = c("exp, 2.5", "eeg, 2.5"),
                    lab = c("EXP", "EEG"),
                    E = c("exp", "eeg"))
p2 <- ggplot() +
    geom_histogram(data = dplt[E=="exp"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    geom_histogram(data = dplt[E=="eeg"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    geom_text(data = danno1, aes(x = x, y = y, label = lab), color = "black", 
            size = 8, family = "Times") + 
    scale_fill_manual(values = pal) +
    geom_vline(data = dv, aes(xintercept = TTA), linetype = "dashed") +
    scale_x_continuous(name = "RT (s)", breaks = c(0, 1, 2, 3, 4, 5, 6)) +
    ylab("") +
    facet_grid(E~., switch = "y") +
    theme_minimal(base_size = 12) +
    # theme(legend.position = "none",
    theme(legend.position = c(.80, .85),
          axis.title.x = element_blank(),
          axis.text.x =  element_blank(),
          legend.title = element_blank(),
          strip.text.y = element_blank(),
          strip.placement = "outside")

danno2 <- data.table(x = c(5.5, 5.5), 
                    y = c(40, 40),
                    Etta = c("exp, 2.5", "eeg, 2.5"),
                    lab = c("EXP", "EEG"),
                    E = c("exp", "eeg"))

p3 <- ggplot() +
    geom_histogram(data = dplt[E=="exp"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    geom_histogram(data = dplt[E=="eeg"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    geom_text(data = danno2, aes(x = x, y = y, label = lab), color = "black", 
            size = 4, family = "Times") + 
    scale_fill_manual(values = Manu::get_pal("Tui")) +
    geom_vline(data = dv, aes(xintercept = TTA), linetype = "dashed") +
    scale_x_continuous(name = "RT (s)", breaks = c(0, 1, 2, 3, 4, 5, 6)) +
    ylab("") +
    facet_grid(Etta~., switch = "y") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
    # theme(legend.position = c(.85, .95),
    strip.text.y = element_blank(),
          legend.title = element_blank(),
          strip.placement = "outside")
# p3

cairo_pdf("inst/docs/2nd_submit/figs/RTDist_exp_eeg.pdf")
# png("inst/docs/2nd_submitfigs/1RTDist_exp_eeg.png", 800, 600)
gridExtra::grid.arrange(p2, p3, ncol = 1)
dev.off()

    

