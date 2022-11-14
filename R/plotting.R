##' Plot choice RT distributions
##' 
##' Plot choice RT distributions
##' 
##' @param x a data frame
##' @param binwidth bin width
##' @param cols key columns 
##' @param print.fig Boolean 
##' 
##' @examples 
##' ## some examples
##' @export
plot_hist <- function(x, binwidth = .1, cols = c("RT", "R"), print.fig = TRUE) {
  # x <- dat
  # x <- data.frame(dat_tmp)
  # binwidth <- .1
  # cols <- c("RT", "R")
  upper <- max(x$RT) + 3 * binwidth
  bk <- seq(0, upper, by = binwidth)
  c1 <- names(table(x$R))[1]
  c2 <- names(table(x$R))[2]
  dA <- x[x$R == c1, cols]
  dB <- x[x$R == c2, cols]
  
  rtA <- dA$RT
  rtB <- dB$RT
  y0 <- cut(rtA, breaks = bk)
  y1 <- cut(rtB, breaks = bk)
  tmp0 <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", y0) ),
                upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", y0) ))
  tmp1 <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", y1) ),
                upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", y1) ))
  dA$mid <- .5*rowSums(tmp0)
  dB$mid <- .5*rowSums(tmp1)
  y1 <- table(dA$mid) / nrow(x)
  y2 <- table(dB$mid) / nrow(x)
  
  d1 <- data.frame(x = as.numeric(names(y1)), y1)
  d2 <- data.frame(x = as.numeric(names(y2)), y2)
  names(d1) <- c("x", "xr", "y")
  names(d2) <- c("x", "xr", "y")
  d1$R <- c1
  d2$R <- c2
  dout <- rbind(d1, d2)
  
  p2 <- ggplot(data = dout) +
    geom_line( aes(x = x, y = y, colour = R), size = 1) +
    geom_point( aes(x = x, y = y, colour = R), size = 2) +
    xlab("RT (s)") +  ylab("Proportion of trials") +
    theme_minimal(base_size = 20) 
  if (print.fig) print(p2)
  return(dout)
}

##' Get trajectory data out of a rModel function
##' 
##' Get trajectory data out of a rModel function
##' 
##' @param x a model output
##' @param model a string identify which model 
##' @import data.table
##' @examples 
##' ## some examples
##' @export
get_traj <- function(x, model = "inhibition") {
  # x <- tmp0
  # require(data.table)
  len <- x@counter[1,1] + 1; 
  len <- ifelse(len >= 10000, 500, len); 
  idx <- 1:len; idx
  if (model == "inhibition") {
    idxm <- c(2:4, 6)
  } else if (model == "leak") {
    idxm <- c(2:4)
  } else if (model == "input") {
    idxm <- c(2:3)
  }
  groups <- slotNames(x)[idxm]
  
  out <- NULL
  for(i in groups) {
    keyvar <- slot(x, i)[,,1]
    dtmp <- data.frame(x = c(idx, idx), 
                       y = c(keyvar[1, idx], keyvar[2, idx]),
                       z = rep(c("1", "2"), each = len))
    dtmp$gp <- i
    out <- rbind(out, dtmp)
  }
  
  return(out)
}
