\name{pedxing}
\title{Multi-phase decision data}
\alias{ped}
\usage{data(ped)}
\format{\strong{ped} stores three data frames, \code{d}, \code{d0} and \code{dat}. \strong{snew} and \strong{sold} are two labelling schemes for participant labels. The multiple data frames are to enable different software to process conveniently. \code{d0} has all the variables being aggregated across participants. \code{dat} reduces the columns. \code{d} is to be processed by \strong{RStan} via \strong{rethinking} package.

    \itemize{
      \item Trial. A label recording also the sequence of trial appearing in the task.
      \item Jitter. A randomly selected time between -0.1 to 0.1 s. This added a small variability to the time-to-arrival so as to prevent responses based on anticipating the imperative stimulus. 
	    \item BeforeCarPassed. A Boolean variable, indicating a response is made before the car
	passes the crossing zone.
	    \item E. A nominal variable, indicating whether a trial is from the behavioural
	study, \strong{exp}, or from the EEG study, \strong{eeg}.
	    \item TTA. Time-to-arrival in second. 2.5, 3, 3.5 or 4 seconds.
	    \item Side. A nominal variable, indicating a trial showing the avator, at the
	beginning, stands on the right or left side of scene. 
	    \item R. A nominal variable, indicating a trial resulting in a safe response or a collision. 
	    \item C. A Boolean variable, indicating a trial resulting in a safe response (TRUE) or a collision (FALSE).
	    \item G. A nominal variable, indicating a response is made before the car passes the crossing zone.
	    \item TTANme. A nominal variable for the time-to-arrival.
	    \item RT. A continuous variable, indicating response times in second.
	    \item s. A nominal variable, indicating different participants.
	    \item D. A continuous variable of vehicle initial distance. This column has
	only one value, i.e., 40 m. This is to calculate vehicle kinematics. 
	    \item DTTA. A nominal label of distance and TTA for the purpose of plotting figures.
	    }
}
	    
  \description{Behaviroual data from human participants doing a road-crossing task.}
\source{
  Institute for Transport Studies, University of Leeds. This work was supported 
  by Engineering and Physical Sciences Research Council (EPSRC) under Grant 
  number, EP/S005056/1.
}

\references{
  
}
\keyword{datasets}
\examples{
data(ped)
dplt <- dat[RT <= 6 & RT > .1]
1 - nrow(dplt) / nrow(dat)


# Simple get_break function
rt <- dplt$RT
binsize <- psych::describe(rt)$se  
bk <- seq(min(rt)-binsize, max(rt)+binsize, by = binsize)

dv <- data.frame(TTA = c(2.5, 3, 3.5, 4, 2.5, 3, 3.5, 4),
                 Etta = c("exp, 2.5", "exp, 3", "exp, 3.5", "exp, 4",
                          "eeg, 2.5", "eeg, 3", "eeg, 3.5", "eeg, 4"))

# sort(unique(DT$TTA))
# [1] 2.5 3.0 3.5 4.0 

dplt$Rgp <- factor( paste0(dplt$G, ", ", dplt$R), 
                    levels = c("after, safe", "after, hit", "before, safe", "before, hit"))
table(dplt$Rgp)

p0 <- ggplot() +
    geom_histogram(data = dplt[E=="exp"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    geom_histogram(data = dplt[E=="eeg"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    ## scale_fill_manual(values = Manu::get_pal("Tui")) +
    geom_vline(data = dv, aes(xintercept = TTA), linetype = "dashed") +
    scale_x_continuous(name = "RT (s)", breaks = c(0, 1, 2, 3, 4, 5, 6)) +
    facet_grid(E~., switch = "y") +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none",
    # theme(legend.position = c(.80, .90),
          legend.title = element_blank(),
          strip.placement = "outside")

dplt$Etta <- factor( paste0(dplt$E, ", ", dplt$TTA), 
                    levels = c("exp, 2.5", "exp, 3", "exp, 3.5", "exp, 4",
                               "eeg, 2.5", "eeg, 3", "eeg, 3.5", "eeg, 4"))

p1 <- ggplot() +
    geom_histogram(data = dplt[E=="exp"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    geom_histogram(data = dplt[E=="eeg"], aes(x = RT, fill = Rgp), 
                   breaks = bk) +
    ## scale_fill_manual(values = Manu::get_pal("Tui")) +
    geom_vline(data = dv, aes(xintercept = TTA), linetype = "dashed") +
    scale_x_continuous(name = "RT (s)", breaks = c(0, 1, 2, 3, 4, 5, 6)) +
    ylab("") +
    facet_grid(Etta~., switch = "y") +
    theme_minimal(base_size = 16) +
    theme(legend.position = c(.85, .95),
          legend.title = element_blank(),
          strip.placement = "outside")

gridExtra::grid.arrange(p0, p1, ncol = 2)

}
    
